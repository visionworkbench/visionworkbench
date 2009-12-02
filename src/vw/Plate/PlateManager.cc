#include <vw/Plate/PlateManager.h>

template <class ViewT>
void vw::platefile::WritePlateFileTask<ViewT>::operator() () { 
  if (m_verbose) 
    std::cout << "\t    Generating tile: [ " << m_tile_info.j << " " << m_tile_info.i 
              << " @ level " <<  m_depth << "]    BBox: " << m_tile_info.bbox << "\n";
  ImageView<typename ViewT::pixel_type> new_data = crop(m_view, m_tile_info.bbox);

  // If this tile contains no data at all, then we bail early without
  // doing anything.
  if (is_transparent(new_data)) {
    m_progress.report_incremental_progress(1.0);
    return;
  }
  
  // If this tile is opaque, then we can plate it directly into the
  // mosaic.
  //
  // XXX TODO: Take care of cases where an opaque tile masks other
  // tiles at higher resolutions in the pyramid.
  if(is_opaque(new_data) ) {

    m_platefile->write(new_data, m_tile_info.i, m_tile_info.j, 
                       m_depth, m_write_transaction_id);


  // If the new_data is not transparent, then we need to go
  // fetch the existing data so that we can composite the two
  // tiles.  This search begins at the current level, but if no
  // tile exists at that level, then we search up the tree for a
  // lower resolution tile that we can safely supersample and
  // crop.
  } else {

    std::cout << "\t--> This tile is NOT opaque!\n";
    bool found = false;
    IndexRecord closest_record;
    int search_depth = m_depth;
    int search_col = m_tile_info.i;
    int search_row = m_tile_info.j;

    // Search up the tree, looking for the nearest VALID tile that
    // exists on top of which we can composite the new data.
    while (!found && search_depth >= 0) {
      try {
        closest_record = m_platefile->read_record(search_col, search_row, 
                                                  search_depth, m_write_transaction_id-1);

        // Mosaicking must happen strictly in order of transaction
        // ID.  Here we check to see if the underlying data tile is
        // available yet.  It may be "locked" by another mosaicking
        // session running on a different instance of image2plate.
        // If that's the case, we sleep for a short time and then
        // try again to obtain the lock.
        while (closest_record.status() == INDEX_RECORD_LOCKED) {
          vw_out(0) << "WAITING for tile [ " << search_col << " " << search_row 
                    << " @ " << m_depth << "]\n";
          sleep(5.0);
          closest_record = m_platefile->read_record(search_col, search_row, 
                                                    search_depth, m_write_transaction_id-1);
        }

        std::cout << "\t    Found potential match: " << closest_record.DebugString() << "\n";

        // If we find a valid tile, the search is over.
        if (closest_record.status() == INDEX_RECORD_VALID) {
          found = true;
        } else {
          --search_depth;
          search_col /= 2;
          search_row /= 2;
        }
      } catch (TileNotFoundErr &e) {
        // If the tile is not found, the search continues at the next level
        --search_depth;
        search_col /= 2;
        search_row /= 2;
      }
    }
      
    // If no tile was found, then we can safely place the raw data
    // into the mosaic directly.
    if (!found) {
      m_platefile->write(new_data, m_tile_info.i, m_tile_info.j, 
                         m_depth, m_write_transaction_id);

    } else {

      ImageView<typename ViewT::pixel_type> old_data(new_data.cols(), new_data.rows());
      m_platefile->read(old_data, search_col, search_row, 
                        search_depth, m_write_transaction_id-1);

      // If we found a valid tile at a lower resolution in the
      // pyramid, it needs to be supersampled and cropped before we
      // can use it here.
      if (search_depth != m_depth) {
        std::cout << "Didn't find tile: [ " << m_tile_info.i << " " 
                  << m_tile_info.j << " @ " << m_depth << "]         but did at [" 
                  << search_col << " " << search_row << " @ " << search_depth << "]\n";
        std::cout << "Tile dimensions : " <<old_data.cols() << " " << old_data.rows() << "\n";
      }

      // Create the image composite and render to an image view
      if (search_depth != m_depth) 
        std::cout << "creating mosaic\n";
      vw::mosaic::ImageComposite<typename ViewT::pixel_type> composite;
      composite.insert(old_data, 0, 0);
      composite.insert(new_data, 0, 0);
      composite.set_draft_mode( true );
      composite.prepare();

      if (search_depth != m_depth) 
        std::cout << "rendering mosaic " << composite.cols() << " " << composite.rows() << "\n";
      ImageView<typename ViewT::pixel_type> composite_tile = composite;
      if (search_depth != m_depth) 
        std::cout << "rendered " << composite_tile.cols() << " " << composite_tile.rows() << "\n";
      
      // Write the result to the platefile
      if (search_depth != m_depth) 
        std::cout << "writing data\n";
      m_platefile->write(composite_tile, m_tile_info.i, m_tile_info.j, 
                         m_depth, m_write_transaction_id);

      if (search_depth != m_depth)
        std::cout << "done.\n";
    }
  }


  // Report progress
  m_progress.report_incremental_progress(1.0);
}


// Explicit template instatiation
namespace vw { 
namespace platefile {

  template
  void vw::platefile::WritePlateFileTask<ImageViewRef<PixelGrayA<uint8> > >::operator() ();

  template
  void vw::platefile::WritePlateFileTask<ImageViewRef<PixelGray<uint8> > >::operator() ();

  template
  void vw::platefile::WritePlateFileTask<ImageViewRef<PixelRGB<uint8> > >::operator() ();

}}
