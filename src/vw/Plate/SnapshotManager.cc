// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/PlateManager.h>    // for bbox_tiles()...
#include <vw/Plate/SnapshotManager.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Image/ImageView.h>
#include <set>

namespace vw {
namespace platefile {

  bool is_leaf(boost::shared_ptr<PlateFile> platefile, TileHeader const& tile_header) {

    int current_col = tile_header.col();
    int current_row = tile_header.row();

    for (int new_row = current_row*2; new_row < current_row*2+2; ++new_row) {
      for (int new_col = current_col*2; new_col < current_col*2+2; ++new_col) {

        try {
          std::list<TileHeader> tile_records;
          tile_records = platefile->search_by_location(new_col, new_row, 
                                                       tile_header.level() + 1,
                                                       tile_header.transaction_id(),
                                                       tile_header.transaction_id(),
                                                       false); // fetch_one_additional_entry
          
          // If this node has child tiles, then it is not a leaf node.
          if (tile_records.size() > 0) 
            return false;
        } catch (TileNotFoundErr &e) { /* do nothing */ }
      }
    }

    // If not child tiles were found, then this is a leaf node.
    return true;
  }

  template <class PixelT> 
  int vw::platefile::SnapshotManager<PixelT>::snapshot_helper(int current_col, 
                       int current_row, 
                       int current_level, 
                       std::map<int32, TileHeader> composite_tiles,
                       vw::BBox2i const& target_region,
                       int target_level, 
                       int start_transaction_id, 
                       int end_transaction_id, 
                       int write_transaction_id) const {

    // Check to see if there are any tiles at this level that need to be
    // snapshotted.  We fetch one additional tile outside of the
    // specified range so that we get any tiles that were not actually
    // composited during the last snapshot.
    std::list<TileHeader> tile_records;
    try {
      tile_records = m_platefile->search_by_location(current_col, current_row, current_level,
                                                   start_transaction_id, end_transaction_id,
                                                   true); // fetch_one_additional_entry
    } catch (TileNotFoundErr &e) { /* do nothing */ }

    // If there are no valid tiles at this level, then there is nothing
    // further for us to do here on this branch of the recursion.  We
    // return immediately.
    if (tile_records.size() == 0)
      return 0;

    // Insert the TileHeaders for this level into the composite_tiles
    // map<>.  Elements will be inserted in order, and any duplicate
    // transaction_ids will be overwritten by new tile records.  In
    // this fashion, we build up a list of the highest resolution
    // valid tiles in order of transaction id.
    for (std::list<TileHeader>::iterator iter = tile_records.begin();
         iter != tile_records.end(); ++iter) {

      if (current_level == target_level || is_leaf(m_platefile, *iter) ) {

        // We always add any valid tiles at the target_level, since
        // they should always be included in the composite.
        //
        // When adding tiles that may need to be supersampled, we are
        // only interested in leaf nodes.  All others will have higher
        // resolution tiles available as their children that we should
        // use instead.
        composite_tiles[iter->transaction_id()] = *iter;

      }
      
    }

    // If we have reached the target level, then we use the
    // accumulated tiles to generate a snapshot for this tile.
    if (current_level == target_level) {

      // Check to see if the second entry is the start_transaction_id.
      // If this is the case, then we can safely skip the first tile
      // (which is then assumed to be the "additional_entry" that
      // falls earlier than the given transaction range) since it
      // will have already been incorporated into the previous
      // snapshot.
      if (composite_tiles.size() >= 2) {
        std::map<int32, TileHeader>::iterator lowest_tid_iter = composite_tiles.begin();
        std::map<int32, TileHeader>::iterator second_tid_iter = lowest_tid_iter;
        ++second_tid_iter;
        if (int(second_tid_iter->first) == start_transaction_id) {
          composite_tiles.erase(lowest_tid_iter);
        } 
      }
      
      // If we arrive at the target level and find that we have fewer
      // than two tiles to composite, then there's no actual
      // compositing that needs to be done.  We can safely bail on
      // this tile.
      if (composite_tiles.size() < 2)
        return 0;

      // Create a tile into which we will accumulate composited data.
      ImageView<PixelT> composite_tile(m_platefile->default_tile_size(),
                                       m_platefile->default_tile_size());

      int num_composited = 0;
      for (std::map<int32, TileHeader>::reverse_iterator iter = composite_tiles.rbegin(); 
           iter != composite_tiles.rend(); ++iter) {

        TileHeader &current_hdr = iter->second;

        // std::cout << "\t--> [ " << current_hdr.transaction_id() << " ]  " 
        //           << current_hdr.col() << " " << current_hdr.row()
        //           << " @ " << current_hdr.level() << "\n";
        
        // Tile access is highly localized during snapshotting, so it
        // speeds things up considerably to cache them here.
        ImageView<PixelT> new_tile;
        if ( !(this->restore_tile(new_tile, current_hdr.col(), current_hdr.row(), 
                                  current_hdr.level(), current_hdr.transaction_id())) ) {
          try {
            // If the tile isn't in our cache, then we take the
            // performance hit and read the tile from the platefile.
            m_platefile->read(new_tile, 
                            current_hdr.col(), current_hdr.row(), 
                            current_hdr.level(), current_hdr.transaction_id(),
                            true); // exact_transaction_match
            this->save_tile(new_tile, current_hdr.col(), current_hdr.row(), 
                            current_hdr.level(), current_hdr.transaction_id());
          } catch (BlobIoErr &e) {
            // If we get a BlobIO error, that's bad news, but not worth
            // killing the snapshot for.  Instead we log the error here and move
            // onto the next location.
            std::ostringstream ostr;
            ostr << "WARNING: error reading tile from blob: " << e.what();
            m_platefile->log(ostr.str());
            continue;
          }

        }
        
        // If this is the first tile in this location, it's already at
        // the target_level, and it's opaque (which is a very common
        // case), then we don't need to save it as part of the
        // snapshot since its already the top tile in the snapshot.
        // This optimization saves us both space and time!
        if ( current_hdr.level() == target_level && 
             iter == composite_tiles.rbegin() &&
             is_opaque(new_tile) ) {
          return 0;
        }

        // If we are compositing a tile from a lower resolution in the
        // pyramid, it needs to be supersampled and cropped before we
        // can use it here.
        if (current_hdr.level() != target_level) {
          
          int level_diff = target_level - current_hdr.level();
          int scaling_factor = pow(2,level_diff);
          int subtile_u = current_col - current_hdr.col() * scaling_factor;
          int subtile_v = current_row - current_hdr.row() * scaling_factor;
          BBox2i subtile_bbox(new_tile.cols() * subtile_u,
                              new_tile.rows() * subtile_v,
                              new_tile.cols(), new_tile.rows());
        
          // Scale up and interpolate the old_tile, then crop out the
          // subtile that we need.
          ImageView<PixelT> subtile = 
            crop(transform(new_tile, ResampleTransform(scaling_factor, scaling_factor), 
                           ConstantEdgeExtension(), BilinearInterpolation()),
                 subtile_bbox);
          
          // Replace the old data
          new_tile = subtile;
        }
        
        // Add the new tile UNDERNEATH the tiles we have already composited.
        vw::mosaic::ImageComposite<PixelT> composite;
        composite.insert(new_tile, 0, 0);
        composite.insert(composite_tile, 0, 0);
        composite.set_draft_mode( true );
        composite.prepare();

        // Overwrite the composite_tile
        composite_tile = composite;
        ++num_composited;

        // Check to see if the composite tile is opaque.  If it is,
        // then there's no point in compositing any more tiles because
        // they will end up hidden anyway.
        if ( is_opaque(composite_tile) ) {
          break;
        }
      }

      if (num_composited > 0) {
        // -- debug
        vw_out(DebugMessage, "plate::snapshot")  << "\t    Compositing " << current_col << " " 
                                                 << current_row << " @ " << current_level 
                                                 << "  [ " << num_composited << " ]\n";
        //---

        m_platefile->write_update(composite_tile, 
                                current_col, current_row, current_level,
                                write_transaction_id);
        return 1;
      } else {
        return 0;
      }

    } else {
      
      // If we haven't reached the target level, then are four
      // possible children that we can explore for more snapshotting
      // information.
      int num_tiles_updated = 0;
      for (int new_row = current_row*2; new_row < current_row*2+2; ++new_row) {
        for (int new_col = current_col*2; new_col < current_col*2+2; ++new_col) {

          // We should only explore those child tiles that overlap
          // with the target_region at the target level.  To test
          // this, we create a bounding box for the child tile we are
          // testing, and scale that bounding box up to see how big it
          // is when projected into the target_level of the pyramid.
          BBox2i parent_region(new_col, new_row, 1, 1);
          parent_region.min() *= pow(2,target_level-(current_level+1));
          parent_region.max() *= pow(2,target_level-(current_level+1));


          // With the scaled up parent region and the target region in
          // hand, we can do a simple intersection test to assess
          // whether this child tile is worth exploring.
          if (parent_region.intersects(target_region)) {
            num_tiles_updated += snapshot_helper(new_col, new_row, current_level+1,
                                                 composite_tiles, target_region, target_level, 
                                                 start_transaction_id, end_transaction_id, 
                                                 write_transaction_id);
          }
        }
      }
      return num_tiles_updated;
    }
  }
}} // namespace vw::platefile


template <class PixelT>
void vw::platefile::SnapshotManager<PixelT>::snapshot(int level, BBox2i const& tile_region, 
                                                      int start_transaction_id, 
                                                      int end_transaction_id, 
                                                      int write_transaction_id) const {
  
  // Subdivide the bbox into smaller workunits if necessary.
  // This helps to keep operations efficient.
  std::list<BBox2i> tile_workunits = bbox_tiles(tile_region, 1024, 1024);
  for ( std::list<BBox2i>::iterator region_iter = tile_workunits.begin(); 
        region_iter != tile_workunits.end(); ++region_iter) {

    // Create an empty list of composite tiles and then kick off the
    // recursive snapshotting process.
    std::map<int32, TileHeader> composite_tiles;
    int num_tiles_updated = snapshot_helper(0, 0, 0, composite_tiles, *region_iter, level, 
                                            start_transaction_id, end_transaction_id, 
                                            write_transaction_id);

    if (num_tiles_updated > 0)
      vw_out() << "\t--> Snapshot " << *region_iter << " @ level " << level 
               << " (" << num_tiles_updated << " tiles updated).\n";

  //   // It will save us time (and result in fewer RPC calls to
  //   // valid_tiles()) if we start the search by popping up several
  //   // levels in the pyramid to first confirm that there is indeed new
  //   // data that needs to be snapshotted.  If there isn't any
  //   // snapshotting to be done at these higher levels, then there
  //   // won't be any snapshotting to do at this level, so we can bail
  //   // early.
  //   int search_level = level-5;
  //   bool worth_continuing = true;
  //   while (search_level != level && search_level >= 0) {

  //     // Subsample this region, returning a region covering the same
  //     // area at a lower level of the pyramid.  
  //     BBox2i level_region = *region_iter;
  //     level_region.min() /= pow(2,level-search_level);
  //     level_region.max().x() = ceil(float(level_region.max().x()) / 
  //                                   powf(2.0,level-search_level));
  //     level_region.max().y() = ceil(float(level_region.max().y()) / 
  //                                   powf(2.0,level-search_level));
      
  //     // Check to see if there are any tiles at this level.  
  //     std::list<TileHeader> tile_records = m_platefile->search_by_region(search_level, 
  //                                                                   level_region,
  //                                                                   start_transaction_id,
  //                                                                   end_transaction_id, 2,
  //                                                                   true); // fetch_one_additional_entry

  //     // If there are, then we continue the search.  If not, then we
  //     // bail early on the search in this particular region.
  //     if (tile_records.size() == 0) {
  //       worth_continuing = false;
  //       break;
  //     }

  //     // Move down the pyramid
  //     search_level++;
  //   }

  //   // It is only worth processing this region if there were tiles at
  //   // higher levels of the pyramid that also need to be mipmapped.
  //   // If our search at higher levels of the pyramid failed, then we
  //   // skip this region.  This optimization prevents us from wasting
  //   // time looking for valid tiles in an area of the mosaic that is
  //   // devoid of any valid tiles.
  //   if (!worth_continuing)
  //     continue;

  //   // Fetch the list of valid tiles in this particular workunit.  
  //   std::list<TileHeader> tiles_in_region = m_platefile->search_by_region(level, *region_iter,
  //                                                                         start_transaction_id,
  //                                                                         end_transaction_id, 2,
  //                                                                         true); // fetch_one_additional_entry

  //   // If there were no valid tiles at *this* level, then we can
  //   // continue without further processing.
  //   if (tiles_in_region.size() == 0) 
  //     continue;

  //   vw_out() << "\t--> Snapshotting " << *region_iter << " @ level " << level 
  //            << ". [ " << tiles_in_region.size() << " snapshottable tiles. ]\n";

  //   // For debugging:
  //   //    if (tiles_in_region.size() != 0)
  //   // std::cout << "\t    Processing Workunit: " << *region_iter 
  //   //           << "    Found " << tiles_in_region.size() << " tile records.\n";    

  //   for ( std::list<TileHeader>::iterator header_iter = tiles_in_region.begin(); 
  //         header_iter != tiles_in_region.end(); ++header_iter) {

  //     // For each location with valid tiles, we query the index for
  //     // the complete list of tiles that match our transaction range
  //     // query.
  //     std::list<TileHeader> tiles_at_location = m_platefile->search_by_location(header_iter->col(),
  //                                                                  header_iter->row(),
  //                                                                  header_iter->level(),
  //                                                                  start_transaction_id,
  //                                                                  end_transaction_id, 
  //                                                                  true); // fetch_one_additional_entry

  //     // Create a tile into which we will accumulate composited data.
  //     ImageView<PixelT> composite_tile(m_platefile->default_tile_size(),
  //                                      m_platefile->default_tile_size());

  //     // Check to see if the second entry is the start_transaction_id.
  //     // If this is the case, then we can safely skip the first tile
  //     // (which is then assumed to be the "additional_entry" that
  //     // falls earlier than the given transaction range."  since it
  //     // will have already been incorporated into the previous
  //     // snapshot.
  //     bool ignore_additional_entry = false;
  //     if (tiles_at_location.size() >= 2) {
  //       std::list<TileHeader>::reverse_iterator lowest_tid_iter = tiles_at_location.rbegin();
  //       std::list<TileHeader>::reverse_iterator second_tid_iter = lowest_tid_iter;
  //       ++second_tid_iter;
  //       if (second_tid_iter->transaction_id() == start_transaction_id) {
  //         ignore_additional_entry = true;
  //       } 
  //     }

  //     //      std::cout << "There are " << tiles_at_location.size() << " tiles.\n";

  //     // Iterate over the tiles at this location.  Start at the top.  
  //     int num_composited = 0;
  //     for ( std::list<TileHeader>::iterator location_iter = tiles_at_location.begin(); 
  //           location_iter != tiles_at_location.end(); ++location_iter) {
  //       //        std::cout << "\tCompositing " << location_iter->transaction_id() << "\n";

  //       // Check the transaction_id and the ignore_additional_entry
  //       // flag to make sure we don't duplicate any effort from the last snapshot.
  //       if (ignore_additional_entry && location_iter->transaction_id() < start_transaction_id) 
  //         continue;

  //       std::cout << "\t--> [ " << location_iter->transaction_id() << " ]  " 
  //                 << location_iter->col() << " " << location_iter->row()
  //                 << " @ " << location_iter->level() << "\n";


  //       // Read the tile from the platefile.
  //       ImageView<PixelT> new_tile;
  //       try {
  //         m_platefile->read(new_tile, 
  //                           header_iter->col(), header_iter->row(), 
  //                           header_iter->level(), location_iter->transaction_id(),
  //                           true); // exact_transaction_match
  //       } catch (BlobIoErr &e) {

  //         // If we get a BlobIO error, that's bad news, but not worth
  //         // killing the snapshot for.  Instead we log the error here and move
  //         // onto the next location.
  //         std::ostringstream ostr;
  //         ostr << "WARNING: error reading tile from blob: " << e.what();
  //         m_platefile->log(ostr.str());
  //         continue;
  //       }
      
  //       // If this is the first tile in the location, and it is opaque
  //       // (which is a very common case), then we don't need to save
  //       // it as part of the snapshot since its already on top.  This
  //       // optimization saves us both space and time!
  //       if ( is_opaque(new_tile) && location_iter == tiles_at_location.begin()) 
  //         break;

  //       // if (debug_me)
  //       //   std::cout << "\t--> Adding " << location_iter->col() << " " << location_iter->row() << "   --   " << location_iter->transaction_id() << "\n";

  //       // Add the new tile UNDERNEATH the tiles we have already composited.
  //       vw::mosaic::ImageComposite<PixelT> composite;
  //       composite.insert(new_tile, 0, 0);
  //       composite.insert(composite_tile, 0, 0);
  //       composite.set_draft_mode( true );
  //       composite.prepare();

  //       // Overwrite the composite_tile
  //       composite_tile = composite;
  //       ++num_composited;

  //       // Check to see if the composite tile is opaque.  If so, then
  //       // there's no point in compositing any more tiles because they
  //       // will end up hidden anyway.
  //       if ( is_opaque(composite_tile) )
  //         break;
  //     }
        
  //     if (num_composited > 1) {
  //       // -- debug
  //       std::cout << "Compositing " << header_iter->col() << " " 
  //                 << header_iter->row() << " @ " << header_iter->level() << "  [ "
  //                 << tiles_at_location.size() << "]\n";
  //       //---

  //       m_platefile->write_update(composite_tile, 
  //                                 header_iter->col(),
  //                                 header_iter->row(),
  //                                 header_iter->level(),
  //                                 write_transaction_id);
  //     }
  //   }
    //    vw_out() << "\t    Region complete.\n";
  }
}

// Create a full snapshot of every level and every region in the mosaic.
template <class PixelT>
void vw::platefile::SnapshotManager<PixelT>::full_snapshot(int start_transaction_id, 
                                                           int end_transaction_id, 
                                                           int write_transaction_id) const {

  for (int level = 0; level < m_platefile->num_levels(); ++level) {    

    // For debugging:
    //    for (int level = 11; level < 12; ++level) {    

    // Snapshot the entire region at each level.  These region will be
    // broken down into smaller work units in snapshot().
    int region_size = pow(2,level);
    int subdivided_region_size = region_size / 16;
    if (subdivided_region_size < 1024) subdivided_region_size = 1024;
    BBox2i full_region(0,0,region_size,region_size);
    std::list<BBox2i> workunits = bbox_tiles(full_region, 
                                             subdivided_region_size, 
                                             subdivided_region_size);
    for ( std::list<BBox2i>::iterator region_iter = workunits.begin(); 
          region_iter != workunits.end(); ++region_iter) {
      snapshot(level, *region_iter, start_transaction_id, 
               end_transaction_id, write_transaction_id);
    }
  }
}

// Explicit template instatiation
namespace vw { 
namespace platefile {
  
  template 
  void vw::platefile::SnapshotManager<PixelGrayA<uint8> >::snapshot(int level, 
                                                                    BBox2i const& bbox, 
                                                                    int start_transaction_id, 
                                                                    int end_transaction_id, 
                                                                    int write_transaction_id) const;
  template 
  void vw::platefile::SnapshotManager<PixelGrayA<int16> >::snapshot(int level, 
                                                                    BBox2i const& bbox, 
                                                                    int start_transaction_id, 
                                                                    int end_transaction_id, 
                                                                    int write_transaction_id) const;
  template
  void vw::platefile::SnapshotManager<PixelRGBA<uint8> >::snapshot(int level, 
                                                                   BBox2i const& bbox, 
                                                                   int start_transaction_id, 
                                                                   int end_transaction_id, 
                                                                   int write_transaction_id) const;

  template
  void vw::platefile::SnapshotManager<PixelGrayA<uint8> >::full_snapshot(int start_transaction_id, 
                                                                         int end_transaction_id, 
                                                                         int write_transaction_id) const;
  template
  void vw::platefile::SnapshotManager<PixelGrayA<int16> >::full_snapshot(int start_transaction_id, 
                                                                         int end_transaction_id, 
                                                                         int write_transaction_id) const;
  template
  void vw::platefile::SnapshotManager<PixelRGBA<uint8> >::full_snapshot(int start_transaction_id, 
                                                                         int end_transaction_id, 
                                                                         int write_transaction_id) const;



}}
