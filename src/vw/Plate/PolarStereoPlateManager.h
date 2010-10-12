// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_POLAR_STEREO_PLATEMANAGER_H__
#define __VW_PLATE_POLAR_STEREO_PLATEMANAGER_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Plate/PlateManager.h>
#include <vector>
#include <sstream>
#include <boost/foreach.hpp>

namespace vw {
namespace platefile {

  class PlateFile;
  class TileInfo;

  // Figures out the tiles effected by image. It's a we bit smarter
  // than just a transformed bbox.
  void stereo_image_tiles( BBox2i const&,
                           cartography::GeoTransform const&, int32,
                           std::list<TileInfo> & );

  // -------------------------------------------------------------------------
  //                            PLATE MANAGER
  // -------------------------------------------------------------------------

  template <class PixelT>
  class PolarStereoPlateManager : public PlateManager {
    FifoWorkQueue m_queue;
    cartography::Datum m_datum;
  public:

  PolarStereoPlateManager(boost::shared_ptr<PlateFile> platefile,
                          cartography::Datum const& datum ) :
    PlateManager(platefile), m_queue(1), m_datum(datum) {}

    /// Add an image to the plate file. Return used transaction id.
    template <class ViewT>
    void insert( ImageViewBase<ViewT> const& imagebase,
                std::string const& description,
                int transaction_id_override,
                cartography::GeoReference const& input_georef,
                bool tweak_settings_for_terrain, bool /*verbose*/ = false,
                const ProgressCallback &progress = ProgressCallback::dummy_instance()) {

      ViewT const& image = imagebase.impl();

      // Determine if input is North or South Pole from points
      Vector2 test_points[5];
      test_points[0] = Vector2( image.cols()/2, image.rows()/2 );
      test_points[1] = Vector2( image.cols()*3/4, image.rows()/2 );
      test_points[2] = Vector2( image.cols()*1/4, image.rows()/2 );
      test_points[3] = Vector2( image.cols()/2, image.rows()*3/4 );
      test_points[4] = Vector2( image.cols()/2, image.rows()*1/4 );
      uint8 north_count = 0;
      for ( uint8 i = 0; i < 5; i++ )
        if ( input_georef.pixel_to_lonlat( test_points[i] )[1] > 0 )
          north_count++;
      bool is_north = north_count > 2;

      // Work out output resolution from 5 points
      cartography::GeoReference output_georef( m_datum );
      output_georef.set_stereographic( is_north ? 90.0 : -90.0, 0, 1.0, 0, 0 );
      Matrix3x3 tmp_transform = math::identity_matrix<3>();
      tmp_transform(1,1) = -1;
      tmp_transform(0,2) = -m_datum.semi_major_axis();
      tmp_transform(1,2) = m_datum.semi_major_axis();
      output_georef.set_transform( tmp_transform );
      cartography::GeoTransform geotx( input_georef, output_georef );
      // We are seeding pixel_per_meters with the lowest resolution possible
      double requested_pixels_per_meters=256.0/(2*m_datum.semi_major_axis());
      for ( uint i = 0; i < 5; i++ ) {
        Vector2 i_pos = geotx.forward( test_points[i] );
        Vector2 x_res = geotx.forward( test_points[i]+Vector2(1,0) ) - i_pos;
        Vector2 y_res = geotx.forward( test_points[i]+Vector2(0,1) ) - i_pos;
        double i_resolution = 1.0 / std::min( norm_2(x_res), norm_2(y_res) );
        if ( i_resolution > requested_pixels_per_meters )
          requested_pixels_per_meters = i_resolution;
      }

      // Fit requested_pixels_per_meters to the nearest (256*2^n) / (2*major)
      int32 pyramid_level =
        boost::numeric_cast<int32>(ceil(log(requested_pixels_per_meters*2*m_datum.semi_major_axis()/256)/log(2)));
      requested_pixels_per_meters =
        256*pow(2,pyramid_level)/(2*m_datum.semi_major_axis());

      Matrix3x3 aff_transform = math::identity_matrix<3>();
      aff_transform(0,0) = 1/requested_pixels_per_meters;
      aff_transform(1,1) = -1/requested_pixels_per_meters;
      aff_transform(0,2) = -m_datum.semi_major_axis();
      aff_transform(1,2) = m_datum.semi_major_axis();
      output_georef.set_transform( aff_transform );

      // Work out which tiles we affect
      geotx = cartography::GeoTransform( input_georef, output_georef );
      BBox2i output_bbox = geotx.forward_bbox( bounding_box(image) );
      vw_out() << "\t    Placing image at level " << pyramid_level
               << " with bbox " << output_bbox << "\n"
               << "\t    (Total Stereographic resolution at this level =  "
               << requested_pixels_per_meters*2*m_datum.semi_major_axis() << " pixels.)\n";
      if ( is_north )
        vw_out() << "\t    This is a North Pole image.\n";
      else
        vw_out() << "\t    This is a South Pole image.\n";

      ImageViewRef<typename ViewT::pixel_type> stereo_view =
        transform(image,geotx,ZeroEdgeExtension(),BicubicInterpolation());

      // chop up the image into small chunks of what actually needs rendering
      std::list<TileInfo> tiles;
      stereo_image_tiles( bounding_box(image), geotx, 256, tiles );

      // Compute the affected tiles.
      BBox2i affected_tiles_bbox;
      BOOST_FOREACH( TileInfo const& tile, tiles ) {
        affected_tiles_bbox.grow(Vector2i(tile.i,tile.j));
      }
      size_t tiles_size = tiles.size();

      // Obtain a transaction ID for this image.
      //
      // Note: the user may have specified a transaction_id to use,
      // which was passed in with transaction_id_override.  If not,
      // then transaction_id_override == -1, and we get an
      // automatically assigned t_id.
      int32 transaction_id =
        m_platefile->transaction_request(description,
                                         transaction_id_override);

      vw_out(InfoMessage, "platefile")
        << "\t    Rasterizing " << tiles_size << " image tiles.\n"
        << "\t    Platefile ID: " << (m_platefile->index_header().platefile_id()) << "\n"
        << "\t    Transaction ID: " << transaction_id << "\n"
        << "\t    Affected tiles @ root: " << affected_tiles_bbox << "\n";

      // Grab a lock on a blob file to use for writing tiles during
      // the two operations below.
      m_platefile->write_request();

      // Add each tile
      progress.report_progress(0);
      BOOST_FOREACH( TileInfo const& tile, tiles ) {
        typedef WritePlateFileTask<ImageViewRef<typename ViewT::pixel_type> > Job;
        m_queue.add_task(boost::shared_ptr<Task>(
          new Job(m_platefile, transaction_id,
                  tile, pyramid_level,
                  stereo_view, tweak_settings_for_terrain,
                  false, tiles_size, progress)));
      }
      m_queue.join_all();
      progress.report_finished();

      // Sync the index
      m_platefile->sync();

      // Mipmap the tiles.
      if (m_platefile->num_levels() > 1) {
        std::ostringstream mipmap_str;
        mipmap_str << "\t--> Mipmapping from level " << pyramid_level << ": ";
        this->mipmap(pyramid_level, affected_tiles_bbox, transaction_id,
                     (!tweak_settings_for_terrain), // mipmap preblur = !tweak_settings_for_terrain
                     TerminalProgressCallback( "plate", mipmap_str.str()));
      }

      // Release the blob id lock.
      m_platefile->write_complete();

      // Notify the index that this transaction is complete.  Do not
      // update the read cursor (false).
      m_platefile->transaction_complete(transaction_id, false);
    }

    /// This function generates a specific mipmap tile at the given
    /// col, row, and level, and transaction_id.
    void generate_mipmap_tile(int col, int row, int level,
                              int transaction_id, bool preblur) const {
      // Create an image large enough to store all of the child nodes
      int tile_size = m_platefile->default_tile_size();
      ImageView<PixelT> super(2*tile_size, 2*tile_size);

      // Iterate over the children, gathering them and (recursively)
      // regenerating them if necessary.
      for( int j=0; j<2; ++j ) {
        for( int i=0; i<2; ++i ) {
          try {
            int child_col = 2*col+i;
            int child_row = 2*row+j;
            vw_out(VerboseDebugMessage, "platefile") << "Reading tile "
              << child_col << " " << child_row
              << " @  " << (level+1) << "\n";
            ImageView<PixelT> child;
            m_platefile->read(child, child_col, child_row, level+1,
                transaction_id, true); // exact_transaction_match == true
            crop(super,tile_size*i,tile_size*j,tile_size,tile_size) = child;
          } catch (TileNotFoundErr &e) {
            // If that fails, then there is no tile. Do nothing to
            // this quadrant of super.
          }
        }
      }

      // We subsample after blurring with a standard 2x2 box filter.
      std::vector<float> kernel(2);
      kernel[0] = kernel[1] = 0.5;

      ImageView<PixelT> new_tile;
      if (preblur) {
        new_tile = subsample( separable_convolution_filter( super,
                                                            kernel,
                                                            kernel,
                                                            1, 1,
                                                            ConstantEdgeExtension() ), 2);
      } else {
        new_tile = subsample( super, 2 );
      }

      if (!is_transparent(new_tile)) {
        vw_out(VerboseDebugMessage, "platefile") << "Writing " << col << " " << row
          << " @ " << level << "\n";
        m_platefile->write_update(new_tile, col, row, level, transaction_id);
      }
    }
  };

}} // end vw::platefile

#endif//__VW_PLATE_POLAR_STEREO_PLATEMANAGER_H__
