// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/gui/TileGenerator.h>

using namespace vw;
using namespace vw::gui;


BBox2i vw::gui::tile_to_bbox(int32 tile_size, int col, int row, int level, int max_level) {
  if (col < 0 || row < 0 || col >= pow(2, max_level) || row >= pow(2, max_level) ) {
    return BBox2i();
  } else {
    BBox2i result(tile_size*col, tile_size*row, tile_size, tile_size);
    return result * pow(2,max_level - level);
  }
}
  
std::list<TileLocator> vw::gui::bbox_to_tiles(int32 tile_size, BBox2i bbox, int level, int max_level) {

  std::list<TileLocator> results;

  // Compute the bounding box at the current level.
  BBox2i level_bbox = bbox / pow(2,max_level - level);

  // Grow that bounding box to align with tile boundaries
  BBox2i aligned_level_bbox = level_bbox;
  aligned_level_bbox.min() = ( (level_bbox.min() / tile_size) * tile_size );
  aligned_level_bbox.max().x() = ( int(ceilf( float(level_bbox.max().x()) / tile_size ))
                                   * tile_size );
  aligned_level_bbox.max().y() = ( int(ceilf( float(level_bbox.max().y()) / tile_size ))
                                   * tile_size );

  int tile_y = aligned_level_bbox.min().y() / tile_size;
  int dest_row = 0;
  while ( tile_y < aligned_level_bbox.max().y() / tile_size ) {
      
    int tile_x = aligned_level_bbox.min().x() / tile_size;
    int dest_col = 0;
    while ( tile_x < aligned_level_bbox.max().x() / tile_size ) {
      BBox2i tile_bbox(dest_col, dest_row, tile_size, tile_size);
      TileLocator loc;
      loc.col = tile_x;
      loc.row = tile_y;
      loc.level = level;
      results.push_back(loc);
        
      ++tile_x;
      dest_col += tile_size;
    }
    ++tile_y;
    dest_row += tile_size;
  }
  return results;
}


boost::shared_ptr<ViewImageResource> DummyTileGenerator::generate_tile(TileLocator const& tile_info) {
  ImageView<PixelRGBA<float> > tile(m_tile_size, m_tile_size);
  for (unsigned j = 0; j < m_tile_size; ++j){
    for (unsigned i = 0; i < m_tile_size; ++i){
      if (abs(i - j) < 10 || abs(i - (m_tile_size - j)) < 10)
        tile(i,j) = PixelRGBA<float>(1.0,0.0,0.0,1.0);
      else 
        tile(i,j) = PixelRGBA<float>(0.0,0.0,0.0,1.0);
    }
  }
  boost::shared_ptr<ViewImageResource> result( new ViewImageResource(tile) );
  return result;
}
