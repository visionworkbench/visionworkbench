#ifndef __VW_PLATE_TOAST_PLATEMANAGER_H__
#define __VW_PLATE_TOAST_PLATEMANAGER_H__

#include <vw/Math/Vector.h>
#include <vw/Image.h>
#include <vw/Mosaic/ImageComposite.h>

#include <vw/Plate/Index.h>
#include <vw/Plate/Blob.h>

#include <vector>
#include <fstream>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>

// Protocol Buffer
#include <vw/Plate/IndexRecord.pb.h>

namespace fs = boost::filesystem;

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------------
  //                            PLATE MANAGER
  // -------------------------------------------------------------------------
  
  class ToastPlateManager {
    boost::shared_ptr<PlateFile> m_platefile;
    FifoWorkQueue m_queue;

    // -------------------------------------------------------------------------
    //                       PLATE FILE TASKS
    // -------------------------------------------------------------------------
    
    template <class ViewT>
    class WritePlateFileTask : public Task {
      boost::shared_ptr<PlateFile> m_parent;
      ViewT const& m_view;
      int m_i, m_j, m_depth;
      BBox2i m_bbox; 
      
    public:
      WritePlateFileTask(boost::shared_ptr<PlateFile> parent, int i, int j, int depth, 
                         ImageViewBase<ViewT> const& view, BBox2i bbox) : 
        m_parent(parent), m_i(i), m_j(j), m_depth(depth), m_view(view.impl()), m_bbox(bbox) {}
      
      virtual ~WritePlateFileTask() {}
      virtual void operator() () { 
        std::cout << "\t    Generating block: [ " << m_j << " " << m_i 
                  << " @ level " <<  m_depth << "]    BBox: " << m_bbox << "\n";

        m_parent->write(m_i, m_j, m_depth, crop(m_view, m_bbox)); }
    };

    // The Block Entry is used to keep track of the bounding box of
    // tiles and their location in the grid.
    struct BlockEntry {
      int i, j;
      BBox2i bbox;
      BlockEntry(int i, int j, BBox2i const& bbox) : i(i), j(j), bbox(bbox) {}
    };

  public:
  
    ToastPlateManager(boost::shared_ptr<PlateFile> platefile) : m_platefile(platefile) {}

    /// The destructor saves the platefile to disk. 
    ~ToastPlateManager() {}

    /// Tiles in the TOAST projection overlap with their neighbors (
    /// the last row and last column of pixels is the same as the top
    /// row and left column of the next images.)  This means that
    /// these bounding boxes are a little funny.  This function
    /// computes those bounding boxes for the tiles at the bottom of
    /// the pyramid so that there is the proper amount of overlap.
    template <class ViewT>
    std::vector<BlockEntry> wwt_image_blocks( ImageViewBase<ViewT> const& image,
                                              int32 tile_size) {
      std::vector<BlockEntry> result;
      int x = 0, y = 0;
      int32 minx = 0, miny = 0;

      while (miny + tile_size <= image.impl().rows()) {
        while (minx + tile_size <= image.impl().cols()) {
          
          BlockEntry be(x, y, BBox2i(minx, miny, tile_size, tile_size));
          result.push_back(be);

          minx += (tile_size-1);
          ++x;
        }
        minx = 0;
        x = 0;
        miny += (tile_size-1);
        ++y;
      }
      return result;
    }

    /// Add an image to the plate file.
    template <class ViewT>
    void insert(ImageViewBase<ViewT> const& image, int max_level) {

      // Do some quick sanity checks
      VW_ASSERT(image.impl().cols()==image.impl().rows(), ArgumentErr() 
                << "TOAST Mosaic  requires a square source image.");
      VW_ASSERT(image.impl().cols()==((m_platefile->default_block_size()-1)*
                                      (1<<(max_level))+1), ArgumentErr() 
                << "TOAST mosaic requires a source image with dimensions (tilesize-1)*2^numlevels+1.");

      // chop up the image into small chunks
      std::vector<BlockEntry> bboxes = wwt_image_blocks( image.impl(), 
                                                         m_platefile->default_block_size());
      
      // And save each block to the PlateFile
      for (int i = 0; i < bboxes.size(); ++i) {
        m_queue.add_task(boost::shared_ptr<Task>(new WritePlateFileTask<ViewT>(m_platefile, 
                                                                               bboxes[i].i, 
                                                                               bboxes[i].j, 
                                                                               max_level,
                                                                               image, 
                                                                               bboxes[i].bbox)));
      }
      m_queue.join_all();
    }


    // Read a previously-written tile in from disk.  Cache the most
    // recently accessed tiles, since each will be used roughly four
    // times.
    ImageView<PixelRGBA<uint8> > load_tile( int32 level, int32 x, int32 y ) {
      int32 num_tiles = 1 << level;
      if( x==-1 ) {
	if( y==-1 ) {
	  return load_tile(level, num_tiles-1, num_tiles-1);
	}
	if( y==num_tiles ) {
	  return load_tile(level, num_tiles-1, 0);
	}
	ImageView<PixelRGBA<uint8> > tile = load_tile(level, 0, num_tiles-1-y);
	if( tile ) return rotate_180(tile);
	else return tile;
      }
      if( x==num_tiles ) {
	if( y==-1 ) {
	  return load_tile(level, 0, num_tiles-1);
	}
	if( y==num_tiles ) {
	  return load_tile(level, 0, 0);
	}
	ImageView<PixelRGBA<uint8> > tile = load_tile(level, num_tiles-1, num_tiles-1-y);
	if( tile ) return rotate_180(tile);
	else return tile;
      }
      if( y==-1 ) {
	ImageView<PixelRGBA<uint8> > tile = load_tile(level, num_tiles-1-x, 0);
	if( tile ) return rotate_180(tile);
	else return tile;
      }
      if( y==num_tiles ) {
	ImageView<PixelRGBA<uint8> > tile = load_tile(level, num_tiles-1-x, num_tiles-1);
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
      
      // Read it in from the platefile
      ImageView<PixelRGBA<uint8> > tile;
      IndexRecord rec;
      try {
        rec = m_platefile->read(tile, x, y, level);
      } catch (TileNotFoundErr &e) {} // Do nothing... the IndexRecord will be invalid. 

      // If the tile does not exist, then we try to regenerate it by
      // fetching our children and rebuilding this tile.
      if (!rec.valid()) {

        std::cout << "\t    Generating mipmap: [ " << x << " " << y << " @ " << level << " ]\n";

        // If none of the termination conditions are met, then we must
        // be at an invalid record that needs to be regenerated.
        int tile_size = m_platefile->default_block_size();

        // Create an image large enough to store all of the child nodes
        ImageView<PixelRGBA<uint8> > super(4*tile_size-3, 4*tile_size-3);
        
        // Iterate over the children, gathering them and (recursively)
        // regenerating them if necessary.
        for( int j=-1; j<3; ++j ) {
          for( int i=-1; i<3; ++i ) {
            ImageView<PixelRGBA<uint8> > child = load_tile(level+1,2*x+i,2*y+j);
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
        
        if( ! is_transparent(tile) ) 
          m_platefile->write(x, y, level, tile);

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

    void mipmap() { 
      ImageView<PixelRGBA<uint8> > tile = this->load_tile(0,0,0); 
      m_platefile->write(0,0,0,tile);
    }
  };


}} // namespace vw::plate

#endif // __VW_PLATE_TOAST_PLATEMANAGER_H__
