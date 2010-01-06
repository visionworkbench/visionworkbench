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
  
  std::cout << "\t--> Snapshotting @ level " << level << "\n";

  // Subdivide the bbox into smaller workunits if necessary.
  // This helps to keep operations efficient.
  std::list<BBox2i> tile_workunits = bbox_tiles(tile_region, 1024, 1024);
  for ( std::list<BBox2i>::iterator region_iter = tile_workunits.begin(); 
        region_iter != tile_workunits.end(); ++region_iter) {
    std::cout << "\t    Processing Workunit: " << *region_iter << "\n";

    // Fetch the list of valid tiles in this particular workunit.  
    std::list<TileHeader> tile_records = m_platefile->valid_tiles(level, *region_iter,
                                                                  start_transaction_id,
                                                                  end_transaction_id);

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
        for (typename std::list<ImageView<PixelT> >::iterator tile_iter = tiles.begin();
             tile_iter != tiles.end();  ++tile_iter) {
          composite.insert(*tile_iter, 0, 0);
        }
        composite.set_draft_mode( true );
        composite.prepare();
        ImageView<PixelT> composited_tile = composite;
        m_platefile->write(composited_tile, 
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
    BBox2i region(0,0,pow(2,level),pow(2,level));
    
    snapshot(level, region, start_transaction_id, end_transaction_id, write_transaction_id);
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
  void vw::platefile::SnapshotManager<PixelRGBA<uint8> >::full_snapshot(int start_transaction_id, 
                                                                         int end_transaction_id, 
                                                                         int write_transaction_id) const;



}}
