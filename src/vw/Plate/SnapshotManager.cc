#include <vw/Plate/PlateManager.h>    // for bbox_tiles()...
#include <vw/Plate/SnapshotManager.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Image/ImageView.h>

// mipmap() generates mipmapped (i.e. low resolution) tiles in the mosaic.
template <class PixelT>
void vw::platefile::SnapshotManager<PixelT>::snapshot(int blob_id, 
                                                      int level, BBox2i const& tile_region, 
                                                      int start_transaction_id, 
                                                      int end_transaction_id, 
                                                      int write_transaction_id) const {
  
  std::cout << "\t--> Snapshotting " << tile_region << " @ level " << level << ".\n";

  // Subdivide the bbox into smaller workunits if necessary.
  // This helps to keep operations efficient.
  std::list<BBox2i> tile_workunits = bbox_tiles(tile_region, 1024, 1024);
  for ( std::list<BBox2i>::iterator region_iter = tile_workunits.begin(); 
        region_iter != tile_workunits.end(); ++region_iter) {

    // Fetch the list of valid tiles in this particular workunit.  
    std::list<TileHeader> tile_records = m_platefile->valid_tiles(level, *region_iter,
                                                                  start_transaction_id,
                                                                  end_transaction_id, 2);

    //    if (tile_records.size() != 0)
      // std::cout << "\t    Processing Workunit: " << *region_iter 
      //           << "    Found " << tile_records.size() << " tile records.\n";    

    for ( std::list<TileHeader>::iterator header_iter = tile_records.begin(); 
          header_iter != tile_records.end(); ++header_iter) {
      
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
        m_platefile->write_update(blob_id,
                                  composited_tile, 
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
void vw::platefile::SnapshotManager<PixelT>::full_snapshot(int blob_id, 
                                                           int start_transaction_id, 
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
      snapshot(blob_id, level, *region_iter, start_transaction_id, 
               end_transaction_id, write_transaction_id);
    }
  }
}

// Explicit template instatiation
namespace vw { 
namespace platefile {
  
  template 
  void vw::platefile::SnapshotManager<PixelGrayA<uint8> >::snapshot(int blob_id, int level, 
                                                                    BBox2i const& bbox, 
                                                                    int start_transaction_id, 
                                                                    int end_transaction_id, 
                                                                    int write_transaction_id) const;
  template
  void vw::platefile::SnapshotManager<PixelRGBA<uint8> >::snapshot(int blob_id, int level, 
                                                                   BBox2i const& bbox, 
                                                                   int start_transaction_id, 
                                                                   int end_transaction_id, 
                                                                   int write_transaction_id) const;

  template
  void vw::platefile::SnapshotManager<PixelGrayA<uint8> >::full_snapshot(int blob_id, 
                                                                         int start_transaction_id, 
                                                                         int end_transaction_id, 
                                                                         int write_transaction_id) const;
  template
  void vw::platefile::SnapshotManager<PixelRGBA<uint8> >::full_snapshot(int blob_id, 
                                                                        int start_transaction_id, 
                                                                         int end_transaction_id, 
                                                                         int write_transaction_id) const;



}}
