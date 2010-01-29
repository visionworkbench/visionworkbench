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
  
  std::cout << "\t--> Snapshotting " << tile_region << " @ level " << level << ".\n";

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
      std::list<TileHeader> tile_records = m_platefile->valid_tiles(search_level, level_region,
                                                                    start_transaction_id,
                                                                    end_transaction_id, 2);

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

    std::cout << "Found some tiles worth considering: " << *region_iter << ".\n";


    // Fetch the list of valid tiles in this particular workunit.  
    std::cout << "Calling valid_tiles with region " << *region_iter << "\n";
    std::list<TileHeader> tile_records = m_platefile->valid_tiles(level, *region_iter,
                                                                  start_transaction_id,
                                                                  end_transaction_id, 2);

    //    if (tile_records.size() != 0)
      // std::cout << "\t    Processing Workunit: " << *region_iter 
      //           << "    Found " << tile_records.size() << " tile records.\n";    

    for ( std::list<TileHeader>::iterator header_iter = tile_records.begin(); 
          header_iter != tile_records.end(); ++header_iter) {

      std::cout << "Calling multi_read.\n";
      
      std::list<ImageView<PixelT> > tiles;
      std::list<TileHeader> headers = m_platefile->multi_read(tiles, 
                                                              header_iter->col(),
                                                              header_iter->row(),
                                                              header_iter->level(),
                                                              start_transaction_id,
                                                              end_transaction_id);

      // If there is only one tile at this location, then there is no
      // need to replace it with a new version in the snapshot.  We
      // only create a new version if there are 2 or more tiles that
      // need to be composited together.
      if (tiles.size() > 1) {
        vw::mosaic::ImageComposite<PixelT> composite;

        // The list of tiles comes sorted from newest to oldest, but
        // we actually want to composite the images in the opposite
        // order.  We use reverse iterators here instead.
        for (typename std::list<ImageView<PixelT> >::reverse_iterator tile_iter = tiles.rbegin();
             tile_iter != tiles.rend();  ++tile_iter) {
          composite.insert(*tile_iter, 0, 0);
        }
        composite.set_draft_mode( true );
        composite.prepare();
        ImageView<PixelT> composited_tile = composite;
        std::cout << "Calling write_update with header\n";

        m_platefile->write_update(composited_tile, 
                                  header_iter->col(),
                                  header_iter->row(),
                                  header_iter->level(),
                                  write_transaction_id);
      }
    }
  }
}

// Create a full snapshot of every level and every region in the mosaic.
template <class PixelT>
void vw::platefile::SnapshotManager<PixelT>::full_snapshot(int start_transaction_id, 
                                                           int end_transaction_id, 
                                                           int write_transaction_id) const {

  for (int level = 0; level < m_platefile->num_levels(); ++level) {    

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
