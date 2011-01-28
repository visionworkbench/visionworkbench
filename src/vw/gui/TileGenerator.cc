// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/gui/TileGenerator.h>
#include <vw/gui/ImageTileGenerator.h>
#include <vw/gui/PlatefileTileGenerator.h>
#include <vw/gui/TestPatternTileGenerator.h>
#include <vw/gui/WebTileGenerator.h>

#include <vw/Plate/HTTPUtils.h>
#include <boost/filesystem/convenience.hpp>

namespace fs = boost::filesystem;
using namespace vw::platefile;

namespace vw { namespace gui {


BBox2i tile_to_bbox(Vector2i tile_size, int col, int row, int level, int max_level) {
  if (col < 0 || row < 0 || col >= (1 << max_level) || row >= (1 << max_level) ) {
    return BBox2i();
  } else {
    BBox2i result(tile_size[0]*col, tile_size[1]*row, tile_size[0], tile_size[1]);
    return result * (1 << (max_level - level));
  }
}

std::list<TileLocator> bbox_to_tiles(Vector2i tile_size, BBox2i bbox, int level, int max_level, int transaction_id, bool exact_transaction_id_match) {
  std::list<TileLocator> results;

  // Compute the bounding box at the current level.
  BBox2i level_bbox = bbox / (1 << (max_level - level));

  // Grow that bounding box to align with tile boundaries
  BBox2i aligned_level_bbox = level_bbox;
  aligned_level_bbox.min().x() = ( (level_bbox.min().x() / tile_size[0]) * tile_size[0] );
  aligned_level_bbox.min().y() = ( (level_bbox.min().y() / tile_size[1]) * tile_size[1] );
  aligned_level_bbox.max().x() = ( int(ceilf( float(level_bbox.max().x()) / float(tile_size[0]) ))
                                   * tile_size[0] );
  aligned_level_bbox.max().y() = ( int(ceilf( float(level_bbox.max().y()) / float(tile_size[1]) ))
                                   * tile_size[1] );

  int tile_y = aligned_level_bbox.min().y() / tile_size[1];
  int dest_row = 0;
  while ( tile_y < aligned_level_bbox.max().y() / tile_size[1] ) {

    int tile_x = aligned_level_bbox.min().x() / tile_size[0];
    int dest_col = 0;
    while ( tile_x < aligned_level_bbox.max().x() / tile_size[0] ) {
      BBox2i tile_bbox(dest_col, dest_row, tile_size[0], tile_size[1]);
      TileLocator loc;
      loc.col = tile_x;
      loc.row = tile_y;
      loc.level = level;
      loc.transaction_id = transaction_id;
      loc.exact_transaction_id_match = exact_transaction_id_match;
      results.push_back(loc);

      ++tile_x;
      dest_col += tile_size[0];
    }
    ++tile_y;
    dest_row += tile_size[1];
  }
  return results;
}

boost::shared_ptr<TileGenerator> TileGenerator::create(std::string filename_) {

  // Remove trailing /
  boost::trim_right_if(filename_, boost::is_any_of("/"));
  Url u(filename_);

  try {
    if (u.scheme() == "http") {
      return boost::shared_ptr<TileGenerator>( new WebTileGenerator(u.string(),17));
    } else if (u.scheme() == "file") {
      if (fs::extension(u.path()) == ".plate")
        return boost::shared_ptr<TileGenerator>( new PlatefileTileGenerator(u.path()) );
      else if (u.path() == "testpattern")
        return boost::shared_ptr<TileGenerator>( new TestPatternTileGenerator(256) );
      else
        return boost::shared_ptr<TileGenerator>( new ImageTileGenerator(u.path()) );
    } else {
      std::cerr << "Could not open " << u << ":\n\t" << "No handler for url scheme " << u.scheme() << std::endl;
    }
  } catch (const vw::Exception& e) {
    std::cerr << "Could not open " << u << ":\n\t" << e.what() << std::endl;
  }
  exit(EXIT_FAILURE);
}

}} // namespace vw::gui
