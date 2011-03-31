// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_SNAPSHOT_MANAGER_H__
#define __VW_PLATE_SNAPSHOT_MANAGER_H__

#include <vw/Plate/IndexData.pb.h>
#include <vw/Plate/Datastore.h>
#include <vw/Plate/FundamentalTypes.h>
#include <vw/Image/ImageView.h>

namespace vw {
namespace platefile {

  class PlateFile;

  template <class PixelT>
  class SnapshotManager {

    // Tile access is highly localized during snapshotting, so it
    // speeds things up considerably to cache them here.
    struct TileCacheEntry {
      int32 level, x, y;
      TransactionOrNeg transaction_id;
      ImageView<PixelT> tile;
    };
    typedef std::list<TileCacheEntry> tile_cache_t;
    mutable tile_cache_t m_tile_cache;

    // Save the tile in the cache.  The cache size of 1000 records was chosen
    // somewhat arbitrarily.
    void save_tile(ImageView<PixelT> &tile, int32 x, int32 y,
                     int32 level, Transaction transaction_id) const {
      //      std::cout << "\t    Saving tile: " << x << " " << y << " @ "
      //                << level << " for " << transaction_id << "\n";
      if( m_tile_cache.size() >= 1000 )
        m_tile_cache.pop_back();
      TileCacheEntry e;
      e.level = level;
      e.x = x;
      e.y = y;
      e.transaction_id = transaction_id;
      e.tile = tile;
      m_tile_cache.push_front(e);
    }

    bool restore_tile(ImageView<PixelT> &tile, int32 x, int32 y,
                      int32 level, Transaction transaction_id) const {
      for( typename tile_cache_t::iterator i=m_tile_cache.begin();
           i!=m_tile_cache.end(); ++i ) {
        if( i->level==level && i->x==x && i->y==y && i->transaction_id == transaction_id ) {
          TileCacheEntry e = *i;
          m_tile_cache.erase(i);
          tile = e.tile;
          //          std::cout << "\t--> Found tile in cache: " << x << " " << y << " @ " << level << " for " << transaction_id << "\n";
          return true;
        }
      }
      // std::cout << "\t    Tile not found: " << x << " " << y << " @ "
      //           << level << " for " << transaction_id << "\n";
      return false;
    }

    int snapshot_helper(uint32 current_col,
                        uint32 current_row,
                        uint32 current_level,
                        vw::BBox2i const& target_region,
                        uint32 target_level,
                        TransactionRange read_transaction_range) const;

  protected:
    boost::shared_ptr<PlateFile> m_platefile;
    bool m_tweak_settings_for_terrain;

  public:

    SnapshotManager(boost::shared_ptr<PlateFile> platefile, bool tweak_settings_for_terrain)
      : m_platefile(platefile), m_tweak_settings_for_terrain(tweak_settings_for_terrain) {}

    // ---------------------------- SNAPSHOTTING --------------------------------

    // snapshot() creates a complete, composited view of the mosaic.
    //
    //   level -- select the pyramid level on which to carry out mipmapping
    //   bbox -- bounding box (in terms of tiles) containing the tiles that need
    //           to be snapshotted at starting_level.  Use to specify affected tiles.
    //   start_transaction_id -- select a transaction_id to use when accessing tiles.
    //   end_transaction_id -- select a transaction_id to use when accessing tiles.
    //
    void snapshot(uint32 level, BBox2i const& bbox, TransactionRange read_transaction_range) const;

    // Create a full snapshot of every level and every region in the mosaic.
    void full_snapshot(TransactionRange read_transaction_range) const;

  };

}} // namespace vw::plate

#endif // __VW_PLATE_SNAPSHOT_MANAGER_H__
