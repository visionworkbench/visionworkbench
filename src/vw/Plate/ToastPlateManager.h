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

#include <vw/Plate/Index.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/ImageComposite.h>
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
  
  class ToastPlateManager : public PlateManager<ToastPlateManager> {

  public:
  
    ToastPlateManager(boost::shared_ptr<PlateFile> platefile, int num_threads) : 
      PlateManager<ToastPlateManager>(platefile, num_threads) {}

    /// Destructor
    virtual ~ToastPlateManager() {}

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
    std::vector<TileInfo> wwt_image_tiles( BBox2i const& image_bbox, 
                                              int32 const resolution,
                                              int32 const tile_size);

    /// Add an image to the plate file.
    template <class ViewT>
    void insert(ImageViewBase<ViewT> const& image, cartography::GeoReference const& georef,
                std::string description,
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
      std::vector<TileInfo> tiles = wwt_image_tiles( output_bbox, resolution,
                                                     m_platefile->default_tile_size());

      // Determine the read and write transaction ids to use for this image.
      int read_transaction_id = m_platefile->transaction_cursor();
      int write_transaction_id = m_platefile->transaction_request(description);
      
      // And save each tile to the PlateFile
      std::cout << "\t    Rasterizing " << tiles.size() << " image tiles.\n";
      progress.report_progress(0);
      for (int i = 0; i < tiles.size(); ++i) {
        m_queue.add_task(boost::shared_ptr<Task>(
          new WritePlateFileTask<ImageViewRef<typename ViewT::pixel_type> >(m_platefile, 
                                                                            read_transaction_id,
                                                                            write_transaction_id,
                                                                            tiles[i], 
                                                                            pyramid_level, 
                                                                            toast_view, 
                                                                            false,
                                                                            tiles.size(),
                                                                            progress)));
      }
      m_queue.join_all();
      m_platefile->transaction_complete(write_transaction_id);
      progress.report_finished();
    }


    // Read a previously-written tile in from disk.  Cache the most
    // recently accessed tiles, since each will be used roughly four
    // times.
    template <class PixelT>
    void load_tile( vw::ImageView<PixelT> &tile, int32 level, int32 x, int32 y ) {
      tile = this->load_tile_impl<PixelT>(level, x, y);
    }

    // Ok. This is one of those really annoying and esoteric c++
    // template problems: we can't call load_tile_impl<> directly from
    // the CRTP superclass because the template appears in the return
    // type of this method.  Instead, we add the extra layer of
    // indirection (load_tile<>, above), which has the return value in
    // the function arguments.  
    template <class PixelT>
    ImageView<PixelT> load_tile_impl( int32 level, int32 x, int32 y ) {
      int32 num_tiles = 1 << level;
      if( x==-1 ) {
        if( y==-1 ) {
          return load_tile_impl<PixelT>(level, num_tiles-1, num_tiles-1);
        }
        if( y==num_tiles ) {
          return load_tile_impl<PixelT>(level, num_tiles-1, 0);
        }
        ImageView<PixelT> tile = load_tile_impl<PixelT>(level, 0, num_tiles-1-y);
        if( tile ) return rotate_180(tile);
        else return tile;
      }
      if( x==num_tiles ) {
        if( y==-1 ) {
          return load_tile_impl<PixelT>(level, 0, num_tiles-1);
        }
        if( y==num_tiles ) {
          return load_tile_impl<PixelT>(level, 0, 0);
        }
        ImageView<PixelT> tile = load_tile_impl<PixelT>(level, num_tiles-1, num_tiles-1-y);
        if( tile ) return rotate_180(tile);
        else return tile;
      }
      if( y==-1 ) {
        ImageView<PixelT> tile = load_tile_impl<PixelT>(level, num_tiles-1-x, 0);
        if( tile ) return rotate_180(tile);
        else return tile;
      }
      if( y==num_tiles ) {
        ImageView<PixelT> tile = load_tile_impl<PixelT>(level, num_tiles-1-x, num_tiles-1);
        if( tile ) return rotate_180(tile);
        else return tile;
      }
    
      // TODO: Reenable cache
      //
      // // Check the cache
      // for( typename cache_t::iterator i=m_cache.begin(); i!=m_cache.end(); ++i ) {
      //   if( i->level==level && i->x==x && i->y==y ) {
      //     CacheEntry e = *i;
      //     m_cache.erase(i);
      //     m_cache.push_front(e);
      //     return e.tile;
      //   }
      // }


      // First we try to access the indexrecord for this tile.  If that
      // fails, then we must be trying to access a node in the tree that
      // simply doesn't exist.  In this case, we create a totally empty
      // tile and return it.
      ImageView<PixelT> tile;
      IndexRecord rec;
      try {
        rec = m_platefile->read_record(x, y, level);
      } catch (TileNotFoundErr &e) {
        return tile;
      }

      // If the record lookup succeded, we look at the current status of
      // the tile to decide what to do next.


      if (rec.status() == INDEX_RECORD_VALID) {

        // CASE 1 : Valid tiles can be returned without any further processing.
        m_platefile->read(tile, x, y, level);
        return tile;

      } else if (rec.status() == INDEX_RECORD_EMPTY || 
                 rec.status() == INDEX_RECORD_STALE) {
    
        // CASE 2 : Empty tiles need to be regenerated from scratch.

        // Create an image large enough to store all of the child nodes
        int tile_size = m_platefile->default_tile_size();
        ImageView<PixelT> super(4*tile_size-3, 4*tile_size-3);
        
        // Iterate over the children, gathering them and (recursively)
        // regenerating them if necessary.
        for( int j=-1; j<3; ++j ) {
          for( int i=-1; i<3; ++i ) {
            ImageView<PixelT> child = load_tile_impl<PixelT>(level+1,2*x+i,2*y+j);
            if( child ) crop(super,(tile_size-1)*(i+1),(tile_size-1)*(j+1),tile_size,tile_size) = child;	    
          }
        }
        
        // In the WWT implementation of TOAST the pixel centers
        // (rather than the than pixel corners) are grid-aligned, so
        // we need to use an odd-sized antialiasing kernel instead of
        // the usual 2x2 box filter.  The following 5-pixel kernel was
        // optimized to avoid the extra blurring associated with using
        // a kernel wider than 2 pixels.  Math was involved.
        std::vector<float> kernel(5);
        kernel[0] = kernel[4] = -0.0344;
        kernel[1] = kernel[3] = 0.2135;
        kernel[2] = 0.6418;

        tile = subsample( crop( separable_convolution_filter( super, 
                                                              kernel, 
                                                              kernel, 
                                                              NoEdgeExtension() ),
                                tile_size-1, tile_size-1, 2*tile_size, 2*tile_size ), 2 );

        if (rec.status() == INDEX_RECORD_STALE) {
          std::cout << "\t    [ " << x << " " << y << " @ " << level << " ] -- Regenerating tile.\n";

          if( ! is_transparent(tile) ) {
            ImageView<PixelT> old_data(tile.cols(), tile.rows());
            try {
              m_platefile->read(old_data, x, y, level);
            } catch (TileNotFoundErr &e) { 
              // Do nothing... we already have a default constructed empty image above! 
            }

            VW_ASSERT(old_data.cols() == tile.cols() && old_data.rows() == tile.rows(),
                      LogicErr() << "WritePlateFileTask::operator() -- new tile dimensions do not " 
                      << "match old tile dimensions.");
        
            vw::platefile::CompositeView<PixelT> composite;
            composite.insert(old_data, 0, 0);
            composite.insert(tile, 0, 0);
            composite.prepare();

            ImageView<PixelT> composite_tile = composite;
            if( ! is_transparent(composite_tile) ) 
              m_platefile->write(composite_tile, x, y, level);
          }

        } else {
          std::cout << "\t    [ " << x << " " << y << " @ " << level << " ] -- Creating tile.\n";
          if( ! is_transparent(tile) ) 
            m_platefile->write(tile, x, y, level);
        }
      }

      // TODO: Reenable cache
      //
      // // Save it in the cache.  The cache size of 1024 tiles was chosen
      // // somewhat arbitrarily.
      // if( m_cache.size() >= 1024 )
      //   m_cache.pop_back();
      // CacheEntry e;
      // e.level = level;
      // e.x = x;
      // e.y = y;
      // e.tile = tile;
      // m_cache.push_front(e);
      return tile;
    }

  };


}} // namespace vw::plate

#endif // __VW_PLATE_TOAST_PLATEMANAGER_H__
