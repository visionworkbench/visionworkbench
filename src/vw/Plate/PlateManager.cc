#include <vw/Plate/PlateManager.h>

// Given a bbox, returns a list of smaller bboxes that perfectly
// tile the space of the larger bbox.
static std::list<vw::BBox2i> bbox_tiles(vw::BBox2i const& bbox, int width, int height) {
  std::list<vw::BBox2i> bboxes;
  
  vw::int32 j_offset = bbox.min().y();
  while ( j_offset < bbox.max().y() ) {
    vw::int32 j_dim = (bbox.max().y() - j_offset) < height ? (bbox.max().y() - j_offset) : height;
    vw::int32 i_offset = bbox.min().x();
    while ( i_offset < bbox.max().x() ) {
      vw::int32 i_dim = (bbox.max().x() - i_offset) < width ? (bbox.max().x() - i_offset) : width;
      bboxes.push_back(vw::BBox2i(i_offset,j_offset,i_dim,j_dim));
      i_offset += i_dim;
    }
    j_offset += j_dim;
  }
  return bboxes;
}


// mipmap() generates mipmapped (i.e. low resolution) tiles in the mosaic.
//
//   starting_level -- select the pyramid level on which to carry out mipmapping
//   ascend_pyramid -- choose whether to build tiles at all pyramid levels (true), or just this one (false).
//   transaction_id -- select a transaction_id to use when accessing tiles.
//   this_transaction_only -- select whether to read tiles for mipmapping using ONLY this t_id (true), or 
//                            mipmap all tiles >= to transaction_id (false).
//   starting_level_bbox -- bounding box (in terms of tiles) containing the tiles that need 
//                          to be mipmapped at starting_level.  Use to specify effected tiles.
//
void vw::platefile::PlateManager::mipmap(int starting_level, bool ascend_pyramid, 
                                         int transaction_id, bool this_transaction_only, 
                                         vw::BBox2i const& bbox) const {
      
      
      // Set the ending level depending on whether or not the user has
      // chosen to build tiles at all pyramid levels, or just this
      // one.
      int ending_level = starting_level;
      if (ascend_pyramid) 
        ending_level = 0;

      BBox2i level_bbox = bbox;
      for ( int level = starting_level; level >= ending_level; --level) {
        std::cout << "\t--> Mipmapping @ level " << level << "\n";

        // Subdivide the bbox into smaller workunits if necessary.
        // This helps to keep operations efficient.
        std::list<BBox2i> tile_workunits = bbox_tiles(level_bbox, 1024, 1024);
        for ( std::list<BBox2i>::iterator iter = tile_workunits.begin();
              iter != tile_workunits.end(); ++iter) {
          std::cout << "\t    Workunit: " << *iter << "\n";

          for (int j = iter->min().y(); j < iter->max().y(); ++j) {
            for (int i = iter->min().x(); i < iter->max().x(); ++i) {
              this->regenerate_tile(i,j,level,transaction_id,transaction_id,true);
            }
          }
        }

        // Adjust the size of the bbox for this level
        level_bbox.min().x() = floor( float(level_bbox.min().x()) / 2 );
        level_bbox.min().y() = floor( float(level_bbox.min().y()) / 2 );
        level_bbox.max().x() = ceil( float(level_bbox.max().x()) / 2 );
        level_bbox.max().y() = ceil( float(level_bbox.max().y()) / 2 );        
      }
    }



template <class PixelT>
vw::ImageView<PixelT> 
vw::platefile::PlateCompositor<PixelT>::composite_mosaic_tile(boost::shared_ptr<PlateFile> platefile, 
                                                           ImageView<PixelT> tile,
                                                           int col, int row, int level,
                                                           int max_level, int transaction_id,
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
    
    // TODO: This is where we could strip the tile of its alpha
    // channel to save space in the placefile.  This will require a
    // view that strips off the alpha channel.
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
        
        bool was_in_cache = this->restore_record(closest_record, search_col, search_row, 
                                                 search_level, transaction_id-1);
        
        if (!was_in_cache)
          closest_record = platefile->read_record(search_col, search_row, 
                                                  search_level, transaction_id-1);

        // Mosaicking must happen strictly in order of transaction
        // ID.  Here we check to see if the underlying data tile is
        // available yet.  It may be "locked" by another mosaicking
        // session running on a different instance of image2plate.
        // If that's the case, we sleep for a short time and then
        // try again to obtain the lock.
        while (search_level + 10 > max_level &&
               (closest_record.status() == INDEX_RECORD_LOCKED || 
                closest_record.status() == INDEX_RECORD_STALE)) {
          vw_out(0) << "\nWAITING for tile [ " << search_col << " " << search_row 
                    << " @ " << search_level << "]\n";
          sleep(5.0);
          closest_record = platefile->read_record(search_col, search_row, 
                                                  search_level, transaction_id-1);
        }

        // Save the tile in the cache if it VALID or EMPTY
        if (closest_record.status() == INDEX_RECORD_VALID ||
            closest_record.status() == INDEX_RECORD_EMPTY) {

          this->save_record(closest_record, search_col, search_row, 
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

        // Save the invalid index record to the cache.  It will become
        // an EMPTY record, but that will result in the same behavior
        // of found = false.
        IndexRecord not_found;
        this->save_record(not_found, search_col, search_row, 
                          search_level, transaction_id-1);

        // If the tile is not found, the search continues at the next level
        --search_level;
        search_col /= 2;
        search_row /= 2;
      }
    }

    // If no tile was found, then we can safely place the raw data
    // into the mosaic directly.
    if (!found) {

      // Write the un-composited tile to disk
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
  PlateCompositor<PixelGrayA<uint8> >::composite_mosaic_tile(boost::shared_ptr<PlateFile> platefile, 
                        ImageView<PixelGrayA<uint8> > tile,
                        int col, int row, int level,
                        int max_level, int transaction_id,
                        const ProgressCallback &progress_callback);

  template ImageView<PixelRGBA<uint8> > 
  PlateCompositor<PixelRGBA<uint8> >::composite_mosaic_tile(boost::shared_ptr<PlateFile> platefile, 
                        ImageView<PixelRGBA<uint8> > tile,
                        int col, int row, int level,
                        int max_level, int transaction_id,
                        const ProgressCallback &progress_callback);

}}
