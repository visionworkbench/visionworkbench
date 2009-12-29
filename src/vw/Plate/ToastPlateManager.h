// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_TOAST_PLATEMANAGER_H__
#define __VW_PLATE_TOAST_PLATEMANAGER_H__

#include <vw/Image.h>
#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>
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

  template <class PixelT>
  class ToastPlateManager : public PlateManager {

    // Private variables
    FifoWorkQueue m_queue;

    struct CacheEntry {
      int32 level, x, y, transaction_id;
      ImageView<PixelT> tile;
    };
    typedef std::list<CacheEntry> cache_t;
    cache_t m_cache;

    // Private methods
    
    // Compute the bounding boxes for [ tile_size x tile_size ] tiles
    // to generate for an image with image_bbox in a TOAST
    // projection space that is [ resolution x resolution ] pixels in
    // size.
    
    // NOTE: Tiles in the TOAST projection overlap with their
    // neighbors ( the last row and last column of pixels is the same
    // as the top row and left column of the next images.)  This
    // means that these bounding boxes are a little funny.  This
    // function computes those bounding boxes for the tiles at the
    // bottom of the pyramid so that there is the proper amount of
    // overlap.
    
    // TODO: This function is slow and not very smart -- it computes
    // all possible bounding boxes before selecting down to the
    // possible handful that overlap with the image_bbox.  It could
    // be made to go MUCH faster by being smart about its search.
    std::vector<TileInfo> wwt_image_tiles( BBox2i const& input_bbox,
                                           cartography::ToastTransform const& toast_tx,
                                           BBox2i const& image_bbox, 
                                           int32 const resolution,
                                           int32 const tile_size) {
      std::vector<TileInfo> result;
      
      // There's no point in starting the search before there is good
      // image data, so we adjust the start point here.
      int32 minx = int(floor(image_bbox.min().x() / (tile_size-1)) * (tile_size-1));
      int32 miny = int(floor(image_bbox.min().y() / (tile_size-1)) * (tile_size-1));
      int x = minx / (tile_size-1);
      int y = miny / (tile_size-1);

      // Iterate over the bounding boxes in the entire TOAST space...
      int curx = minx;
      int cury = miny;
      while (cury < image_bbox.max().y() - 1) {
        while (curx < image_bbox.max().x() - 1) {
      
          TileInfo be(x, y, BBox2i(curx, cury, tile_size, tile_size));
      
          // ...but only add bounding boxes that overlap with the image.
          if (image_bbox.intersects(be.bbox)) {

            // Images that cross the edges of the TOAST space have
            // very, very large image_bbox's (sometimes containing the
            // whole space!)  We take each individual tile under
            // consideration, and transform it the OTHER way to make
            // sure it actually does still intersect with the source
            // imagery.
            //
            // approximate == true for reverse_bbox() to speed things
            // up.
            if (input_bbox.intersects(toast_tx.reverse_bbox(be.bbox, true))) {
              result.push_back(be);
            } 
          }

          curx += (tile_size-1);
          ++x;
        }
        curx = minx;
        x = minx / (tile_size-1);
        cury += (tile_size-1);
        ++y;
      }
      return result;
    }
 
    // Read a previously-written tile in from disk.  Cache the most
    // recently accessed tiles, since each will be used roughly four
    // times.
    void load_tile( vw::ImageView<PixelT> &tile, int32 level, int32 x, int32 y, 
                    int transaction_id, int max_level ) {
      tile = this->load_tile_impl(level, x, y, transaction_id, max_level);
    }

    // Ok. This is one of those really annoying and esoteric c++
    // template problems: we can't call load_tile_impl<> directly from
    // the CRTP superclass because the template appears in the return
    // type of this method.  Instead, we add the extra layer of
    // indirection (load_tile<>, above), which has the return value in
    // the function arguments.  
    ImageView<PixelT> load_tile_impl( int32 level, int32 x, int32 y, 
                                      int transaction_id, int max_level );

 public:
  
    ToastPlateManager(boost::shared_ptr<PlateFile> platefile, int num_threads) : 
      PlateManager(platefile), m_queue(num_threads) {}

    /// Add an image to the plate file.
    template <class ViewT>
    void insert(ImageViewBase<ViewT> const& image, std::string const& description,
                cartography::GeoReference const& georef, bool verbose = false,
                const ProgressCallback &progress = ProgressCallback::dummy_instance()) {

      // Compute the pyramid level at which to store this image.  The
      // number of required levels is broken down so that the very top
      // of the pyramid covers the entire globe and has a size of
      // tile_size.
      int tile_size = m_platefile->default_tile_size();
      Vector2 p0 = georef.pixel_to_lonlat(Vector2(image.impl().cols()/2,
                                                  image.impl().rows()/2));
      Vector2 p1 = georef.pixel_to_lonlat(Vector2(image.impl().cols()/2+1,
                                                  image.impl().rows()/2));
      Vector2 p2 = georef.pixel_to_lonlat(Vector2(image.impl().cols()/2,
                                                  image.impl().rows()/2+1));
      double delta = sqrt(pow(p1.y()-p0.y(),2)+pow(p2.y()-p0.y(),2));
      int pyramid_level = (int)round(log(360/delta/(tile_size-1)) /
                                     log(2));
      
      // Compute the resolution of the TOAST output space at the given
      // pyramid_level.  The formula below was carefully chosen to
      // ensure that the output space at each level is sized just
      // slightly shy of an even power of two; therby allowing for the
      // slight overlap of one line of pixels on the right and the
      // bottom of a toast tile.  This overlap is necessary for proper
      // positioning when these tiles are rendered as texture in a 3D
      // graphics environment.
      int32 resolution = (1<<pyramid_level)*(tile_size-1)+1;

      // Set up the toast transform and compute the bounding box of this
      // image in the toast projection space.
      cartography::ToastTransform toast_tx( georef, resolution );
      BBox2i input_bbox = BBox2i(0,0,image.impl().cols(),image.impl().rows());
      BBox2i output_bbox = toast_tx.forward_bbox(input_bbox);

      std::cout << "\t    Placing image at level " << pyramid_level 
                << " with bbox " << output_bbox << "\n"
                << "\t    (Total TOAST resolution at this level =  " 
                << resolution << " pixels.)\n";

      // Create the output view and crop it to the proper size.
      ImageViewRef<typename ViewT::pixel_type> toast_view = 
        transform(image,toast_tx, ZeroEdgeExtension(),BicubicInterpolation());

      if( georef.proj4_str()=="+proj=longlat" &&
          fabs(georef.lonlat_to_pixel(Vector2(-180,0)).x()) < 1 &&
          fabs(georef.lonlat_to_pixel(Vector2(180,0)).x() - image.impl().cols()) < 1 &&
          fabs(georef.lonlat_to_pixel(Vector2(0,90)).y()) < 1 &&
          fabs(georef.lonlat_to_pixel(Vector2(0,-90)).y() - image.impl().rows()) < 1 ) {
        vw_out(0) << "\t--> Detected global overlay.  " 
                  << "Using cylindrical edge extension to hide the seam.\n";
        toast_view = transform(image,toast_tx,
                               CylindricalEdgeExtension(),BicubicInterpolation());
      } 

      // chop up the image into small chunks
      std::vector<TileInfo> tiles = wwt_image_tiles( input_bbox, toast_tx, output_bbox, 
                                                     resolution, 
                                                     m_platefile->default_tile_size());

      // Obtain a transaction ID for this image.  To do so, we must
      // create a list of root tiles that we will be writing so that
      // the index can go ahead and mark those tiles a "being
      // processed."  Marking tiles as part of a transaction_request
      // is an atomic operation on the index, and is critical for
      // ensuring that concurrent image mosaics don't trample on each
      // other by informing the index which tiles will be (eventually)
      // modified under this transaction id.
      std::vector<TileHeader> tile_headers;
      BBox2i effected_tiles_bbox;
      for (size_t i = 0; i < tiles.size(); ++i) {
        TileHeader hdr;
        hdr.set_col(tiles[i].i);
        hdr.set_row(tiles[i].j);
        hdr.set_level(pyramid_level);
        tile_headers.push_back(hdr);
        effected_tiles_bbox.grow(Vector2i(tiles[i].i,tiles[i].j));
      }
      int platefile_id = m_platefile->index_header().platefile_id();
      int transaction_id = m_platefile->transaction_request(description, tile_headers);
      std::cout << "\t    Rasterizing " << tiles.size() << " image tiles.\n" 
                << "\t    Platefile ID: " << platefile_id << "\n"
                << "\t    Transaction ID: " << transaction_id << "\n"
                << "\t    Effected tiles @ root: " << effected_tiles_bbox << "\n";


      // // For debugging: 
      // // 
      // // Test: terminate clients half the time
      // srandom(time(0));
      // float r = float(random()) / (powf(2.0,31)-1.0);
      
      // if (r > 0.2) {
      //   vw_out(0) << "\n\n***********************************************************\n";
      //   vw_out(0) << "                          FAILING...\n";
      //   vw_out(0) << "***********************************************************\n";
      //   vw_throw( IOErr() << "terminating randomly....");
      // }

      // Add each tile.
      progress.report_progress(0);
      for (size_t i = 0; i < tiles.size(); ++i) {
        m_queue.add_task(boost::shared_ptr<Task>(
          new WritePlateFileTask<ImageViewRef<typename ViewT::pixel_type> >(m_platefile, 
                                                                            transaction_id,
                                                                            tiles[i], 
                                                                            pyramid_level, 
                                                                            toast_view, 
                                                                            verbose, // verbose
                                                                            tiles.size(),
                                                                            progress)));
      }
      m_queue.join_all();
      progress.report_finished();

      // Mipmap the tiles.
      if (m_platefile->num_levels() > 1) {
        std::cout << "\t--> Generating mipmap tiles for transaction_id " << transaction_id << "\n";
      
        // Adjust the size of the bbox for this level
        effected_tiles_bbox.min().x() = floor( float(effected_tiles_bbox.min().x()) / 2.0 );
        effected_tiles_bbox.min().y() = floor( float(effected_tiles_bbox.min().y()) / 2.0 );
        effected_tiles_bbox.max().x() = ceil( float(effected_tiles_bbox.max().x()+1) / 2.0 );
        effected_tiles_bbox.max().y() = ceil( float(effected_tiles_bbox.max().y()+1) / 2.0 );        
        
        this->mipmap(m_platefile->num_levels()-2, true, transaction_id, false, effected_tiles_bbox);
      }

      // Notify the index that this transaction is complete.
      m_platefile->transaction_complete(transaction_id);
    }

    virtual void regenerate_tile(int col, int row, 
                                 int level, int read_transaction_id, 
                                 int write_transaction_id, 
                                 bool this_transaction_only) const;

    ImageView<PixelT> composite_child_tile(int &num_composited, 
                                           int col, int row, 
                                           int level, 
                                           int read_transaction_id,
                                           int write_transaction_id,
                                           bool this_transaction_only) const;

  };


}} // namespace vw::plate

#endif // __VW_PLATE_TOAST_PLATEMANAGER_H__
