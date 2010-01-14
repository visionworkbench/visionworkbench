// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_KML_PLATEMANAGER_H__
#define __VW_PLATE_KML_PLATEMANAGER_H__

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

// Protocol Buffer
#include <vw/Plate/ProtoBuffers.pb.h>

namespace fs = boost::filesystem;

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------------
  //                            PLATE MANAGER
  // -------------------------------------------------------------------------
  
  class KmlPlateManager {

    boost::shared_ptr<PlateFile> m_platefile;
    FifoWorkQueue m_queue;

  public:
  
    KmlPlateManager(boost::shared_ptr<PlateFile> platefile, int num_threads) : 
      m_platefile(platefile), m_queue(num_threads)  {}

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
    void insert(ImageViewBase<ViewT> const& image, 
                cartography::GeoReference const& input_georef,
                int read_transaction_id, int write_transaction_id,
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
      std::vector<TileInfo> tiles = kml_image_tiles( output_bbox, resolution,
                                                     m_platefile->default_tile_size());

      // And save each tile to the PlateFile
      std::cout << "\t    Rasterizing " << tiles.size() << " image tiles.\n";
      progress.report_progress(0);
      for (size_t i = 0; i < tiles.size(); ++i) {
        m_queue.add_task(boost::shared_ptr<Task>(
          new WritePlateFileTask<ImageViewRef<typename ViewT::pixel_type> >(m_platefile, 
                                                                            read_transaction_id,
                                                                            write_transaction_id,
                                                                            tiles[i], 
                                                                            pyramid_level, 
                                                                            kml_view,
                                                                            false,
                                                                            tiles.size(),
                                                                            progress)));
      }
      m_queue.join_all();
      progress.report_finished();
    }

    // Read a previously-written tile in from disk.  Cache the most
    // recently accessed tiles, since each will be used roughly four
    // times.
    template <class PixelT>
    void load_tile( vw::ImageView<PixelT> &tile, int32 level, int32 x, int32 y, 
                    int read_transaction_id, int write_transaction_id ) {
    
      // First we try to access the indexrecord for this tile.  If that
      // fails, then we must be trying to access a node in the tree that
      // simply doesn't exist.  In this case, we create a totally empty
      // tile with zero pixels and return it.
      tile.reset();
      IndexRecord rec;
      try {
        rec = m_platefile->read_record(x, y, level, write_transaction_id);
      } catch (TileNotFoundErr &e) {
        return;
      }

      // If the record lookup succeded, we look at the current status of
      // the tile to decide what to do next.
      if (rec.status() == INDEX_RECORD_VALID) {

        // CASE 1 : Valid tiles can be returned without any further processing.
        m_platefile->read(tile, x, y, level, write_transaction_id);
        return;

      } else if (rec.status() == INDEX_RECORD_EMPTY || 
                 rec.status() == INDEX_RECORD_STALE) {
    
        // CASE 2 : Empty tiles need to be regenerated from scratch.

        // Create an image large enough to store all of the child nodes
        int tile_size = m_platefile->default_tile_size();
        ImageView<PixelT> super(2*tile_size, 2*tile_size);
        
        // Iterate over the children, gathering them and (recursively)
        // regenerating them if necessary.
        for( int j=0; j<2; ++j ) {
          for( int i=0; i<2; ++i ) {
            ImageView<PixelT> child;
            load_tile(child,level+1,2*x+i,2*y+j,read_transaction_id,write_transaction_id);
            if( child ) crop(super,tile_size*i,tile_size*j,tile_size,tile_size) = child; 
          }
        }
        
        // We subsample after blurring with a standard 2x2 box filter.
        std::vector<float> kernel(2);
        kernel[0] = kernel[1] = 0.5;
    
        tile = subsample( separable_convolution_filter( super, 
                                                        kernel, 
                                                        kernel, 
                                                        1, 1,
                                                        ConstantEdgeExtension() ), 2);
    
        if (rec.status() == INDEX_RECORD_STALE) {
          std::cout << "\t    [ " << x << " " << y << " @ " << level << " ] -- Regenerating tile.\n";
      
          if( ! is_transparent(tile) ) {
            ImageView<PixelT> old_data(tile.cols(), tile.rows());
            try {
              m_platefile->read(old_data, x, y, level, read_transaction_id);
            } catch (TileNotFoundErr &e) { 
              // Do nothing... we already have a default constructed empty image above! 
            }

            VW_ASSERT(old_data.cols() == tile.cols() && old_data.rows() == tile.rows(),
                      LogicErr() << "WritePlateFileTask::operator() -- new tile dimensions do not " 
                      << "match old tile dimensions.");
        
            vw::mosaic::ImageComposite<PixelT> composite;
            composite.insert(old_data, 0, 0);
            composite.insert(tile, 0, 0);
            composite.set_draft_mode( true );
            composite.prepare();
      
            ImageView<PixelT> composite_tile = composite;
            if( ! is_transparent(composite_tile) ) 
              m_platefile->write(composite_tile, x, y, level, write_transaction_id);
          }
      
        } else {
          std::cout << "\t    [ " << x << " " << y << " @ " << level << " ] -- Creating tile.\n";
          if( ! is_transparent(tile) ) 
            m_platefile->write(tile, x, y, level,write_transaction_id);
        }
      }

      return;
    }

  };


}} // namespace vw::plate

#endif // __VW_PLATE_KML_PLATEMANAGER_H__
