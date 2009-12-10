// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_PLATEMANAGER_H__
#define __VW_PLATE_PLATEMANAGER_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/ToastTransform.h>
#include <vw/Mosaic/ImageComposite.h>

#include <vw/Plate/Index.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/PlateFile.h>

#include <vector>
#include <fstream>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>

// Protocol Buffer
#include <vw/Plate/ProtoBuffers.pb.h>

namespace fs = boost::filesystem;

namespace vw {
namespace platefile {


  // The Tile Entry is used to keep track of the bounding box of
  // tiles and their location in the grid.
  struct TileInfo {
    int i, j;
    BBox2i bbox;
    TileInfo(int i, int j, BBox2i const& bbox) : i(i), j(j), bbox(bbox) {}
  };

  // -------------------------------------------------------------------------
  //                              TILE COMPOSITING
  // -------------------------------------------------------------------------

  template <class PixelT>
  class PlateCompositor {
    
    struct RecordCacheEntry {
      int32 level, x, y, transaction_id;
      IndexRecord record;
    };
    typedef std::list<RecordCacheEntry> record_cache_t;
    record_cache_t m_record_cache;

    // Save the tile in the cache.  The cache size of 10000 records was chosen
    // somewhat arbitrarily.
    void save_record(IndexRecord record, int32 x, int32 y, int32 level, int32 transaction_id) {
      if( m_record_cache.size() >= 10000 )
        m_record_cache.pop_back();
      RecordCacheEntry e;
      e.level = level;
      e.x = x;
      e.y = y;
      e.transaction_id = transaction_id;
      e.record = record;
      m_record_cache.push_front(e);
    }

    bool restore_record(IndexRecord &record, int32 x, int32 y, int32 level, int32 transaction_id) {
      for( typename record_cache_t::iterator i=m_record_cache.begin(); i!=m_record_cache.end(); ++i ) {
        if( i->level==level && i->x==x && i->y==y && i->transaction_id == transaction_id ) {
          RecordCacheEntry e = *i;
          m_record_cache.erase(i);
          record = e.record;
          return true;
        }
      }
      return false;
    }

  public:

    // Composite into the mosaic. The composite_mosaic_tile() function
    // looks for any tiles at equal or lower resolution in the mosaic,
    // and composites this tile on top of those tiles, supersampling the
    // low-res tile if necessary.
    ImageView<PixelT> composite_mosaic_tile(boost::shared_ptr<PlateFile> platefile, 
                                            ImageView<PixelT> tile,
                                            int col, int row, int level,
                                            int max_depth, int transaction_id,
                                            const ProgressCallback &progress_callback = ProgressCallback::dummy_instance());

  };

  // -------------------------------------------------------------------------
  //                            WRITE PLATEFILE TASK
  // -------------------------------------------------------------------------

  template <class ViewT>
  class WritePlateFileTask : public Task {
    boost::shared_ptr<PlateFile> m_platefile;
    int m_transaction_id;
    TileInfo m_tile_info;
    int m_depth;
    ViewT const& m_view;
    bool m_verbose;
    SubProgressCallback m_progress;
    PlateCompositor<typename ViewT::pixel_type> m_compositor;
      
  public:
    WritePlateFileTask(boost::shared_ptr<PlateFile> platefile, 
                       int transaction_id,
                       TileInfo const& tile_info, 
                       int depth, ImageViewBase<ViewT> const& view,
                       bool verbose, int total_num_blocks, 
                       const ProgressCallback &progress_callback = ProgressCallback::dummy_instance()) :
      m_platefile(platefile), m_transaction_id(transaction_id),
      m_tile_info(tile_info), m_depth(depth), m_view(view.impl()), 
      m_verbose(verbose), m_progress(progress_callback,0.0,1.0/float(total_num_blocks)) {}
      
    virtual ~WritePlateFileTask() {}
    virtual void operator() () {
      if (m_verbose) 
        std::cout << "\t    Generating tile: [ " << m_tile_info.j << " " << m_tile_info.i 
                  << " @ level " <<  m_depth << "]    BBox: " << m_tile_info.bbox << "\n";

      // Generate the tile from the image data
      ImageView<typename ViewT::pixel_type> tile = crop(m_view, m_tile_info.bbox);
      
      // Composite into the mosaic. The composite_mosaic_tile() function
      // looks for any tiles at equal or lower resolution in the mosaic,
      // and composites this tile on top of those tiles, supersampling the
      // low-res tile if necessary.
      m_compositor.composite_mosaic_tile(m_platefile, tile, m_tile_info.i, m_tile_info.j, 
                                         m_depth, m_depth, m_transaction_id, m_progress);
    }
  };

}} // namespace vw::plate

#endif // __VW_PLATE_PLATEMANAGER_H__
