// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_PLATE_CARREE_PLATEMANAGER_H__
#define __VW_PLATE_PLATE_CARREE_PLATEMANAGER_H__

#include <vw/Image.h>
#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Cartography/ToastTransform.h>
#include <vw/Mosaic/ImageComposite.h>

#include <vw/Plate/Index.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/PlateManager.h>

#include <vector>
#include <fstream>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------------
  //                            PLATE MANAGER
  // -------------------------------------------------------------------------

  template <class PixelT>
  class PlateCarreePlateManager : public PlateManager {

    //    PixelT m_fix;
    FifoWorkQueue m_queue;

  public:

    PlateCarreePlateManager(boost::shared_ptr<PlateFile> platefile) :
      PlateManager(platefile), m_queue(1)  {} // Set threads to 1 for now...

    // Create a georeference object for this plate file.  The user
    // supplies the desired level for which they want the
    // georeference.
    cartography::GeoReference georeference(int level) {
      int tile_size = m_platefile->default_tile_size();
      int resolution = (1<<level)*tile_size;

      cartography::GeoReference output_georef =
        cartography::output::kml::get_output_georeference(resolution,resolution);
      return output_georef;
    }

    /// Compute the bounding boxes for [ tile_size x tile_size ] tiles
    /// to generate for an image with image_bbox in a TOAST
    /// projection space that is [ resolution x resolution ] pixels in
    /// size.
    ///
    /// NOTE: Tiles in the TOAST projection overlap with their
    /// neighbors ( the last row and last column of pixels is the same
    /// as the top row and left column of the next images.)  This
    /// means that these bounding boxes are a little funny.  This
    /// function computes those bounding boxes for the tiles at the
    /// bottom of the pyramid so that there is the proper amount of
    /// overlap.
    ///
    /// TODO: This function is slow and not very smart -- it computes
    /// all possible bounding boxes before selecting down to the
    /// possible handful that overlap with the image_bbox.  It could
    /// be made to go MUCH faster by being smart about its search.
    std::vector<TileInfo> kml_image_tiles( BBox2i const& image_bbox,
                                           int32 const resolution,
                                           int32 const tile_size);

    /// Add an image to the plate file.
    template <class ViewT>
    void insert(ImageViewBase<ViewT> const& image, std::string const& description,
                int transaction_id_override, cartography::GeoReference const& input_georef,
                bool /*verbose*/ = false,
                const ProgressCallback &progress = ProgressCallback::dummy_instance()) {

      // We build the ouput georeference from the input image's
      // georeference by copying the datum.
      cartography::GeoReference output_georef;
      output_georef.set_datum(input_georef.datum());
      int resolution = 256;

      // Right now, we only need a WGS84 output geoereference to compute
      // the resolution. The rest of the output info will get set later.
      cartography::GeoTransform geotx( input_georef, output_georef );

      // Calculate the best resolution at 5 different points in the image,
      // as occasionally there's a singularity at the center pixel that
      // makes it extremely tiny (such as in pole-centered images).
      const int cols = image.impl().cols();
      const int rows = image.impl().rows();
      Vector2 res_pixel[5];
      res_pixel[0] = Vector2( cols/2, rows/2 );
      res_pixel[1] = Vector2( cols/2 + cols/4, rows/2 );
      res_pixel[2] = Vector2( cols/2 - cols/4, rows/2 );
      res_pixel[3] = Vector2( cols/2, rows/2 + rows/4 );
      res_pixel[4] = Vector2 (cols/2, rows/2 - rows/4 );
      int res;
      for(int i=0; i < 5; i++) {
        res = cartography::output::kml::compute_resolution(geotx, res_pixel[i]);
        if( res > resolution ) resolution = res;
      }

      // Round the resolution to the nearest power of two.  The
      // base of the pyramid is 2^8 or 256x256 pixels.
      int pyramid_level = (int)ceil(log(resolution) / log(2)) - 8;
      int tile_size = m_platefile->default_tile_size();
      resolution = (1<<pyramid_level)*tile_size;

      output_georef = cartography::output::kml::get_output_georeference(resolution,resolution);
      output_georef.set_datum( input_georef.datum() );

      // Set up the KML transform and compute the bounding box of this
      // image in the KML projection space.
      cartography::GeoTransform kml_tx( input_georef, output_georef );
      BBox2i input_bbox = BBox2i(0,0,image.impl().cols(),image.impl().rows());
      BBox2i output_bbox = kml_tx.forward_bbox(input_bbox);
      output_bbox.crop( BBox2i(0,0,resolution,resolution) );

      vw_out() << "\t    Placing image at level " << pyramid_level
               << " with bbox " << output_bbox << "\n"
               << "\t    (Total KML resolution at this level =  "
               << resolution << " pixels.)\n";

      // Create the output view and crop it to the proper size.
      ImageViewRef<typename ViewT::pixel_type> kml_view =
        transform(image,kml_tx, ZeroEdgeExtension(),BicubicInterpolation());

      if( input_georef.proj4_str()=="+proj=longlat" &&
          fabs(input_georef.lonlat_to_pixel(Vector2(-180,0)).x()) < 1 &&
          fabs(input_georef.lonlat_to_pixel(Vector2(180,0)).x() - image.impl().cols()) < 1 &&
          fabs(input_georef.lonlat_to_pixel(Vector2(0,90)).y()) < 1 &&
          fabs(input_georef.lonlat_to_pixel(Vector2(0,-90)).y() - image.impl().rows()) < 1 ) {
        vw_out() << "\t--> Detected global overlay.  "
                 << "Using cylindrical edge extension to hide the seam.\n";
        kml_view = transform(image,kml_tx,
                             CylindricalEdgeExtension(), BicubicInterpolation());
      }

      // chop up the image into small chunks
      // Keep the "this"! It makes kml_image_tiles dependent on template type.
      // http://www.parashift.com/c++-faq-lite/templates.html#faq-35.19
      std::vector<TileInfo> tiles = this->kml_image_tiles( output_bbox,
                                                           resolution,
                                                           m_platefile->default_tile_size());


      // Compute the affected tiles.
      BBox2i affected_tiles_bbox;
      for (size_t i = 0; i < tiles.size(); ++i)
        affected_tiles_bbox.grow(Vector2i(tiles[i].i,tiles[i].j));

      // Obtain a transaction ID for this image.
      //
      // Note: the user may have specified a transaction_id to use,
      // which was passed in with transaction_id_override.  If not,
      // then transaction_id_override == -1, and we get an
      // automatically assigned t_id.
      int transaction_id = m_platefile->transaction_request(description,
                                                            transaction_id_override);

      std::cout << "\t    Rasterizing " << tiles.size() << " image tiles.\n"
                << "\t    Platefile ID: " << (m_platefile->index_header().platefile_id()) << "\n"
                << "\t    Transaction ID: " << transaction_id << "\n"
                << "\t    Affected tiles @ root: " << affected_tiles_bbox << "\n";

      // Grab a lock on a blob file to use for writing tiles during
      // the two operations below.
      m_platefile->write_request();

      // Add each tile
      progress.report_progress(0);
      for (size_t i = 0; i < tiles.size(); ++i) {
        m_queue.add_task(boost::shared_ptr<Task>(
          new WritePlateFileTask<ImageViewRef<typename ViewT::pixel_type> >(m_platefile,
                                                                            transaction_id,
                                                                            tiles[i],
                                                                            pyramid_level,
                                                                            kml_view,
                                                                            false,
                                                                            tiles.size(),
                                                                            progress)));
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
    virtual void generate_mipmap_tile(int col, int row, int level, int transaction_id) const;

  };


}} // namespace vw::plate

#endif // __VW_PLATE_PLATE_CARREE_PLATEMANAGER_H__
