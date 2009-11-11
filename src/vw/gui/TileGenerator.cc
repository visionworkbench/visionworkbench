// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/gui/TileGenerator.h>
using namespace vw;
using namespace vw::gui;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

// --------------------------------------------------------------------------------
//                                Utility Functions
// --------------------------------------------------------------------------------

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



// --------------------------------------------------------------------------
//                              TILE GENERATOR
// --------------------------------------------------------------------------

boost::shared_ptr<TileGenerator> TileGenerator::create(std::string filename) {
  // If ends in .plate, then assume platefile.
  if ( fs::extension(filename) == ".plate") {

    return boost::shared_ptr<TileGenerator>( new PlatefileTileGenerator(filename) );

  // If testpattern, then we use the testpattern tile generator
  } else if (filename == "testpattern") {

    vw_out(0) << "\t--> Starting vwv in testpattern mode.\n";
    return boost::shared_ptr<TileGenerator>( new TestPatternTileGenerator(256) );
    
  // Otherwise, assume an image.
  } else {
    vw_out(0) << "\t--> Loading image: " << filename << ".\n";
    return boost::shared_ptr<TileGenerator>( new ImageTileGenerator(filename) );
  }
}


// --------------------------------------------------------------------------
//                         TEST PATTERN TILE GENERATOR
// --------------------------------------------------------------------------

boost::shared_ptr<ViewImageResource> TestPatternTileGenerator::generate_tile(TileLocator const& tile_info) {
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

int TestPatternTileGenerator::cols() const { return 2048; }
int TestPatternTileGenerator::rows() const { return 2048; }
PixelFormatEnum TestPatternTileGenerator::pixel_format() const { return VW_PIXEL_RGBA; }
ChannelTypeEnum TestPatternTileGenerator::channel_type() const { return VW_CHANNEL_FLOAT32; }
int TestPatternTileGenerator::tile_size() const { return m_tile_size; }

// --------------------------------------------------------------------------
//                         PLATE FILE TILE GENERATOR
// --------------------------------------------------------------------------

PlatefileTileGenerator::PlatefileTileGenerator(std::string platefile_name) :
  m_platefile(new vw::platefile::PlateFile(platefile_name)) {
  std::cout << "\t--> Loading platefile \"" << platefile_name << "\" with " 
            << m_platefile->depth() << " levels.\n";
}

boost::shared_ptr<ViewImageResource> PlatefileTileGenerator::generate_tile(TileLocator const& tile_info) {
  
  switch (this->pixel_format()) {
  case VW_PIXEL_GRAY:
  case VW_PIXEL_GRAYA:
    if (this->channel_type() == VW_CHANNEL_UINT8) {
      ImageView<PixelGrayA<uint8> > tile;
      m_platefile->read(tile, tile_info.col, tile_info.row, tile_info.level);
      return boost::shared_ptr<ViewImageResource>( new ViewImageResource(tile) );
    } if (this->channel_type() == VW_CHANNEL_INT16) {
      ImageView<PixelGrayA<int16> > tile;
      m_platefile->read(tile, tile_info.col, tile_info.row, tile_info.level);
      return boost::shared_ptr<ViewImageResource>( new ViewImageResource(tile) );
    } if (this->channel_type() == VW_CHANNEL_UINT16) {
      ImageView<PixelGrayA<uint16> > tile;
      m_platefile->read(tile, tile_info.col, tile_info.row, tile_info.level);
      return boost::shared_ptr<ViewImageResource>( new ViewImageResource(tile) );
    } if (this->channel_type() == VW_CHANNEL_FLOAT32) {
      ImageView<PixelGrayA<float> > tile;
      m_platefile->read(tile, tile_info.col, tile_info.row, tile_info.level);
      return boost::shared_ptr<ViewImageResource>( new ViewImageResource(tile) );
    } else {
      std::cout << "This platefile has a channel type that is not yet support by vwv.\n";
      std::cout << "Exiting...\n\n";
      exit(0);
    }
      
    break;

  case VW_PIXEL_RGB:
  case VW_PIXEL_RGBA:
    if (this->channel_type() == VW_CHANNEL_UINT8) {
      ImageView<PixelRGBA<uint8> > tile;
      m_platefile->read(tile, tile_info.col, tile_info.row, tile_info.level);
      return boost::shared_ptr<ViewImageResource>( new ViewImageResource(tile) );
    } if (this->channel_type() == VW_CHANNEL_FLOAT32) {
      ImageView<PixelRGBA<float> > tile;
      m_platefile->read(tile, tile_info.col, tile_info.row, tile_info.level);
      return boost::shared_ptr<ViewImageResource>( new ViewImageResource(tile) );
    } else {
      std::cout << "This platefile has a channel type that is not yet support by vwv.\n";
      std::cout << "Exiting...\n\n";
      exit(0);
    }
      
    break;

  default:
    std::cout << "This platefile has a pixel format that is not yet support by vwv.\n";
    std::cout << "Exiting...\n\n";
    exit(0);
  }
  
  vw_throw(NoImplErr() << "Unsupported pixel format or channel type in TileGenerator.\n");
}

int PlatefileTileGenerator::cols() const {
  return this->tile_size() * pow(2, m_platefile->depth());
}

int PlatefileTileGenerator::rows() const {
  return this->tile_size() * pow(2, m_platefile->depth());
}

PixelFormatEnum PlatefileTileGenerator::pixel_format() const {
  return m_platefile->pixel_format();
}

ChannelTypeEnum PlatefileTileGenerator::channel_type() const {
  return m_platefile->channel_type();
}

int PlatefileTileGenerator::tile_size() const {
  return m_platefile->default_tile_size();
}


// --------------------------------------------------------------------------
//                             IMAGE TILE GENERATOR
// --------------------------------------------------------------------------

boost::shared_ptr<ViewImageResource> ImageTileGenerator::generate_tile(TileLocator const& tile_info) {

}

int ImageTileGenerator::cols() const {
}

int ImageTileGenerator::rows() const {

}

PixelFormatEnum ImageTileGenerator::pixel_format() const {

}

ChannelTypeEnum ImageTileGenerator::channel_type() const {

}

int ImageTileGenerator::tile_size() const {

}
