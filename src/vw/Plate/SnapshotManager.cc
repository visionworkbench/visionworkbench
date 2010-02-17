// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/PlateManager.h>    // for bbox_tiles()...
#include <vw/Plate/SnapshotManager.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Image/ImageView.h>

// mipmap() generates mipmapped (i.e. low resolution) tiles in the mosaic.
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

    // It will save us time (and result in fewer RPC calls to
    // valid_tiles()) if we start the search by popping up several
    // levels in the pyramid to first confirm that there is indeed new
    // data that needs to be snapshotted.  If there isn't any
    // snapshotting to be done at these higher levels, then there
    // won't be any snapshotting to do at this level, so we can bail
    // early.
    int search_level = level-5;
    bool worth_continuing = true;
    while (search_level != level && search_level >= 0) {

      // Subsample this region, returning a region covering the same
      // area at a lower level of the pyramid.  
      BBox2i level_region = *region_iter;
      level_region.min() /= pow(2,level-search_level);
      level_region.max().x() = ceil(float(level_region.max().x()) / 
                                    powf(2.0,level-search_level));
      level_region.max().y() = ceil(float(level_region.max().y()) / 
                                    powf(2.0,level-search_level));
      
      // Check to see if there are any tiles at this level.  
      std::list<TileHeader> tile_records = m_platefile->search_by_region(search_level, 
                                                                    level_region,
                                                                    start_transaction_id,
                                                                    end_transaction_id, 2,
                                                                    true); // fetch_one_additional_entry

      // If there are, then we continue the search.  If not, then we
      // bail early on the search in this particular region.
      if (tile_records.size() == 0) {
        worth_continuing = false;
        break;
      }

      // Move down the pyramid
      search_level++;
    }

    // If our search at higher levels of the pyramid failed, then we
    // skip this region.  This optimization prevents us from wasting
    // time looking for valid tiles in an area of the mosaic that is
    // devoid of any valid tiles.
    if (!worth_continuing)
      continue;

    // Fetch the list of valid tiles in this particular workunit.  
    std::list<TileHeader> tiles_in_region = m_platefile->search_by_region(level, *region_iter,
                                                                          start_transaction_id,
                                                                          end_transaction_id, 2,
                                                                          true); // fetch_one_additional_entry

    // If there were no valid tiles at *this* level, then we can
    // continue also.
    if (tiles_in_region.size() == 0) 
      continue;

    vw_out() << "\t--> Snapshotting " << *region_iter << " @ level " << level 
             << ". [ " << tiles_in_region.size() << " snapshottable tiles. ]\n";

    // For debugging:
    //    if (tiles_in_region.size() != 0)
    // std::cout << "\t    Processing Workunit: " << *region_iter 
    //           << "    Found " << tiles_in_region.size() << " tile records.\n";    

    for ( std::list<TileHeader>::iterator header_iter = tiles_in_region.begin(); 
          header_iter != tiles_in_region.end(); ++header_iter) {

      // For each location with valid tiles, we query the index for
      // the complete list of tiles that match our transaction range
      // query.
      std::list<TileHeader> tiles_at_location = m_platefile->search_by_location(header_iter->col(),
                                                                   header_iter->row(),
                                                                   header_iter->level(),
                                                                   start_transaction_id,
                                                                   end_transaction_id, 
                                                                   true); // fetch_one_additional_entry


      // Create a tile into which we will accumulate composited data.
      ImageView<PixelT> composite_tile(m_platefile->default_tile_size(),
                                       m_platefile->default_tile_size());

      // Check to see if the second entry is the start_transaction_id.
      // If this is the case, then we can safely skip the first tile
      // (which is then assumed to be the "additional_entry" that
      // falls earlier than the given transaction range."  since it
      // will have already been incorporated into the previous
      // snapshot.
      bool ignore_additional_entry = false;
      if (tiles_at_location.size() >= 2) {
        std::list<TileHeader>::reverse_iterator lowest_tid_iter = tiles_at_location.rbegin();
        std::list<TileHeader>::reverse_iterator second_tid_iter = lowest_tid_iter;
        ++second_tid_iter;
        if (second_tid_iter->transaction_id() == start_transaction_id) {
          ignore_additional_entry = true;
        } 
      }

      //      std::cout << "There are " << tiles_at_location.size() << " tiles.\n";

      // Iterate over the tiles at this location.  Start at the top.  
      int num_composited = 0;
      for ( std::list<TileHeader>::iterator location_iter = tiles_at_location.begin(); 
            location_iter != tiles_at_location.end(); ++location_iter) {
        //        std::cout << "\tCompositing " << location_iter->transaction_id() << "\n";

        // Check the transaction_id and the ignore_additional_entry
        // flag to make sure we don't duplicate any effort from the last snapshot.
        if (ignore_additional_entry && location_iter->transaction_id() < start_transaction_id) 
          continue;

        // Read the tile from the platefile.
        ImageView<PixelT> new_tile;
        try {
          m_platefile->read(new_tile, 
                            header_iter->col(), header_iter->row(), 
                            header_iter->level(), location_iter->transaction_id(),
                            true); // exact_transaction_match
        } catch (BlobIoErr &e) {

          // If we get a BlobIO error, that's bad news, but not worth
          // killing the snapshot for.  Instead we log the error here and move
          // onto the next location.
          std::ostringstream ostr;
          ostr << "WARNING: error reading tile from blob: " << e.what();
          m_platefile->log(ostr.str());
          continue;
        }
      
        // If this is the first tile in the location, and it is opaque
        // (which is a very common case), then we don't need to save
        // it as part of the snapshot since its already on top.  This
        // optimization saves us both space and time!
        if ( is_opaque(new_tile) && location_iter == tiles_at_location.begin()) 
          break;

        // if (debug_me)
        //   std::cout << "\t--> Adding " << location_iter->col() << " " << location_iter->row() << "   --   " << location_iter->transaction_id() << "\n";

        // Add the new tile UNDERNEATH the tiles we have already composited.
        vw::mosaic::ImageComposite<PixelT> composite;
        composite.insert(new_tile, 0, 0);
        composite.insert(composite_tile, 0, 0);
        composite.set_draft_mode( true );
        composite.prepare();

        // Overwrite the composite_tile
        composite_tile = composite;
        ++num_composited;

        // Check to see if the composite tile is opaque.  If so, then
        // there's no point in compositing any more tiles because they
        // will end up hidden anyway.
        if ( is_opaque(composite_tile) )
          break;
      }
        
      if (num_composited > 1) {
        // -- debug
        // std::cout << "Compositing " << header_iter->col() << " " 
        //           << header_iter->row() << " @ " << header_iter->level() << "  [ "
        //           << tiles_at_location.size() << "]\n";
        //---

        m_platefile->write_update(composite_tile, 
                                  header_iter->col(),
                                  header_iter->row(),
                                  header_iter->level(),
                                  write_transaction_id);
      }
    }
    vw_out() << "\t    Region complete.\n";
  }
}

// Create a full snapshot of every level and every region in the mosaic.
template <class PixelT>
void vw::platefile::SnapshotManager<PixelT>::full_snapshot(int start_transaction_id, 
                                                           int end_transaction_id, 
                                                           int write_transaction_id) const {

  //  for (int level = 0; level < m_platefile->num_levels(); ++level) {    
  for (int level = 0; level < 1; ++level) {    

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
