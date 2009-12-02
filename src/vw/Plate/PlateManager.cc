#include <vw/Plate/PlateManager.h>

template <class PixelT>
vw::ImageView<PixelT> vw::platefile::composite_mosaic_tile(boost::shared_ptr<PlateFile> platefile, 
                                          ImageView<PixelT> tile,
                                          int col, int row, int level,
                                          int transaction_id,
                                          const ProgressCallback &progress_callback) {
  
  // If this tile contains no data at all, then we bail early without
  // doing anything.
  if (is_transparent(tile)) {
    progress_callback.report_incremental_progress(1.0);
    return tile;
  }
  
  // Store the result.  This first assignment is a shallow copy.
  ImageView<PixelT> result_tile = tile;

  // If this tile is opaque, then we can add it directly into the
  // mosaic without any compositing.
  //
  // XXX TODO: Take care of cases where an opaque tile masks other
  // tiles at higher resolutions in the pyramid.
  if(is_opaque(tile) ) {

    platefile->write(tile, col, row, level, transaction_id);

  // If the tile is not transparent, then we need to go
  // fetch the existing data so that we can composite the two
  // tiles.  This search begins at the current level, but if no
  // tile exists at that level, then we search up the tree for a
  // lower resolution tile that we can safely supersample and
  // crop.
  } else {

    bool found = false;
    IndexRecord closest_record;
    int search_level = level;
    int search_col = col;
    int search_row = row;

    // Search up the tree, looking for the nearest VALID tile that
    // exists on top of which we can composite the new data.
    while (!found && search_level >= 0) {
      try {
        closest_record = platefile->read_record(search_col, search_row, 
                                                  search_level, transaction_id-1);

        // Mosaicking must happen strictly in order of transaction
        // ID.  Here we check to see if the underlying data tile is
        // available yet.  It may be "locked" by another mosaicking
        // session running on a different instance of image2plate.
        // If that's the case, we sleep for a short time and then
        // try again to obtain the lock.
        while (closest_record.status() == INDEX_RECORD_LOCKED || 
               closest_record.status() == INDEX_RECORD_STALE) {
          vw_out(0) << "WAITING for tile [ " << search_col << " " << search_row 
                    << " @ " << search_level << "]\n";
          sleep(5.0);
          closest_record = platefile->read_record(search_col, search_row, 
                                                    search_level, transaction_id-1);
        }

        // If we find a valid tile, the search is over.
        if (closest_record.status() == INDEX_RECORD_VALID) {
          found = true;
        } else {
          --search_level;
          search_col /= 2;
          search_row /= 2;
        }
      } catch (TileNotFoundErr &e) {
        // If the tile is not found, the search continues at the next level
        --search_level;
        search_col /= 2;
        search_row /= 2;
      }
    }
      
    // If no tile was found, then we can safely place the raw data
    // into the mosaic directly.
    if (!found) {

      platefile->write(tile, col, row, level, transaction_id);

    } else {

      ImageView<PixelT> old_tile(tile.cols(), tile.rows());
      platefile->read(old_tile, search_col, search_row, search_level, transaction_id-1);

      // If we found a valid tile at a lower resolution in the
      // pyramid, it needs to be supersampled and cropped before we
      // can use it here.
      if (search_level != level) {

        int level_diff = level - search_level;
        int scaling_factor = pow(2,level_diff);
        int subtile_u = col - search_col * scaling_factor;
        int subtile_v = row - search_row * scaling_factor;
        BBox2i subtile_bbox(old_tile.cols() * subtile_u,
                            old_tile.rows() * subtile_v,
                            old_tile.cols(), old_tile.rows());
        
        // std::cout << "Didn't find tile: [ " << col << " " 
        //           << row << " @ " << level << "]         but did at [" 
        //           << search_col << " " << search_row << " @ " << search_level << "]\n";
        // std::cout << "Tile dimensions : " <<old_tile.cols() << " " << old_tile.rows() << "\n";
        // std::cout << "Subtile : " << subtile_u << " " << subtile_v << "   " << subtile_bbox << "\n";
        // Scale up and interpolate the old_tile, then crop out the
        // subtile that we need.
        ImageView<PixelT> subtile = 
          crop(transform(old_tile, ResampleTransform(scaling_factor, scaling_factor), 
                         ConstantEdgeExtension(), BicubicInterpolation()),
               subtile_bbox);
        
        // Replace the old data
        old_tile = subtile;
      }

      // Create the image composite and render to an image view
      vw::mosaic::ImageComposite<PixelT> composite;
      composite.insert(old_tile, 0, 0);
      composite.insert(tile, 0, 0);
      composite.set_draft_mode( true );
      composite.prepare();
      result_tile = composite;
      
      // Write the result to the platefile
      platefile->write(result_tile, col, row, level, transaction_id);
    }
  }

  // Report progress
  progress_callback.report_incremental_progress(1.0);
  return result_tile;
}


// Explicit template instatiation
namespace vw { 
namespace platefile {

  template ImageView<PixelGrayA<uint8> > 
  composite_mosaic_tile(boost::shared_ptr<PlateFile> platefile, 
                        ImageView<PixelGrayA<uint8> > tile,
                        int col, int row, int level,
                        int transaction_id,
                        const ProgressCallback &progress_callback);

  template ImageView<PixelRGBA<uint8> > 
  composite_mosaic_tile(boost::shared_ptr<PlateFile> platefile, 
                        ImageView<PixelRGBA<uint8> > tile,
                        int col, int row, int level,
                        int transaction_id,
                        const ProgressCallback &progress_callback);

}}
