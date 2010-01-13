#include <vw/Plate/PlateManager.h>
#include <vw/Image/Transform.h>

// Given a bbox, returns a list of smaller bboxes that perfectly
// tile the space of the larger bbox.
std::list<vw::BBox2i> vw::platefile::bbox_tiles(vw::BBox2i const& bbox, int width, int height) {
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
void vw::platefile::PlateManager::mipmap(int starting_level, vw::BBox2i const& bbox, 
                                         int transaction_id, const ProgressCallback &progress_callback) const {


  // Adjust the size of the bbox for the first mipmapping level, which
  // is one level down from the starting_level.  We also pad it
  // slightly here, since we want to make sure we don't trucate at
  // tiles when we divide by 2.
  BBox2i level_bbox = bbox;
  level_bbox.min().x() = floor( float(level_bbox.min().x()) / 2.0 );
  level_bbox.min().y() = floor( float(level_bbox.min().y()) / 2.0 );
  level_bbox.max().x() = ceil( float(level_bbox.max().x()+1) / 2.0 );
  level_bbox.max().y() = ceil( float(level_bbox.max().y()+1) / 2.0 );
    
  // Compute the range of progress for this SubProgressCallback.  The
  // geometric sum below probably could be computed more easily, but
  // I'm tired and a little lazy at the moment.
  float total_num_tiles = 0.0;
  float sum_denom = 4.0;
  for ( int level = starting_level-1; level >= 0; --level) {
    total_num_tiles += (bbox.width() * bbox.height()) / sum_denom;
    sum_denom *= 4.0;
  }
  
  float current_num_tiles = 0;
  sum_denom = 4.0;
  float prev_num_tiles = 0;
  for ( int level = starting_level-1; level >= 0; --level) {

    // Do a little progress callback math.
    current_num_tiles += (bbox.width() * bbox.height()) / sum_denom;
    sum_denom *= 4.0;
    SubProgressCallback sub_progress(progress_callback, 
                                     prev_num_tiles / total_num_tiles,
                                     current_num_tiles / total_num_tiles);
    prev_num_tiles = current_num_tiles;

    // Subdivide the bbox into smaller workunits if necessary.
    // This helps to keep operations efficient.
    std::list<BBox2i> tile_workunits = bbox_tiles(level_bbox, 256, 256);
    int prog_counter = 0;
    for ( std::list<BBox2i>::iterator iter = tile_workunits.begin(); iter != tile_workunits.end(); ++iter) {
      SubProgressCallback sub_sub_progress(sub_progress, 
                                           float(prog_counter)/tile_workunits.size(),
                                           float(prog_counter+1)/tile_workunits.size());
      prog_counter++;

      // The original bbox passed into the mipmap function only serves
      // as a hint to where we might find tiles.  Some tiles will be
      // missing due to transparency, and it's also possible that a
      // very large number of tiles might be missing if the bounding
      // box happens to cross one edge of the mosaic.  In both cases,
      // the call to valid_tiles() here fetches the tiles that
      // actually contain good data.
      BBox2i parent_region = *iter;
      parent_region.min() *= 2;
      parent_region.max() *= 2;
      std::list<TileHeader> valid_tile_records = m_platefile->valid_tiles(level+1, parent_region,
                                                                          transaction_id,
                                                                          transaction_id, 
                                                                          1);

      if (valid_tile_records.size() != 0) {

        // Once we compute the valid tiles at the parent level, we
        // translate that back down into valid tiles at this level.
        BBox2i trimmed_region;
        for (std::list<TileHeader>::iterator trim_iter = valid_tile_records.begin(); 
             trim_iter != valid_tile_records.end(); ++trim_iter) {
          trimmed_region.grow(Vector2i(trim_iter->col()/2, trim_iter->row()/2));
        }

        // Pad the bottom right edge of the bbox to make sure we get
        // tiles from the parent that only hald overlapped with the
        // child.
        trimmed_region.max().x() += 1;
        trimmed_region.max().y() += 1;
        
        for (int j = trimmed_region.min().y(); j < trimmed_region.max().y(); ++j) {
          for (int i = trimmed_region.min().x(); i < trimmed_region.max().x(); ++i) {
            this->generate_mipmap_tile(i,j,level,transaction_id);
            sub_sub_progress.report_incremental_progress(1.0/(trimmed_region.width() * trimmed_region.height()));
          }
        }
      }
      sub_sub_progress.report_finished();
    }
    sub_progress.report_finished();
    
    // Adjust the size of the bbox for this level
    level_bbox.min().x() = floor( float(level_bbox.min().x()) / 2 );
    level_bbox.min().y() = floor( float(level_bbox.min().y()) / 2 );
    level_bbox.max().x() = ceil( float(level_bbox.max().x()) / 2 );
    level_bbox.max().y() = ceil( float(level_bbox.max().y()) / 2 );        
  }
  progress_callback.report_finished();
}
