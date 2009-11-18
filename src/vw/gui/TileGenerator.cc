// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/gui/TileGenerator.h>
#include <vw/Image/ViewImageResource.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/MaskViews.h>
using namespace vw;
using namespace vw::gui;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

// --------------------------------------------------------------------------------
//                                Utility Functions
// --------------------------------------------------------------------------------

BBox2i vw::gui::tile_to_bbox(Vector2i tile_size, int col, int row, int level, int max_level) {
  if (col < 0 || row < 0 || col >= pow(2, max_level) || row >= pow(2, max_level) ) {
    return BBox2i();
  } else {
    BBox2i result(tile_size[0]*col, tile_size[1]*row, tile_size[0], tile_size[1]);
    return result * pow(2,max_level - level);
  }
}
  
std::list<TileLocator> vw::gui::bbox_to_tiles(Vector2i tile_size, BBox2i bbox, int level, int max_level) {
  std::list<TileLocator> results;

  // Compute the bounding box at the current level.
  BBox2i level_bbox = bbox / pow(2,max_level - level);

  // Grow that bounding box to align with tile boundaries
  BBox2i aligned_level_bbox = level_bbox;
  aligned_level_bbox.min().x() = ( (level_bbox.min().x() / tile_size[0]) * tile_size[0] );
  aligned_level_bbox.min().y() = ( (level_bbox.min().y() / tile_size[1]) * tile_size[1] );
  aligned_level_bbox.max().x() = ( int(ceilf( float(level_bbox.max().x()) / tile_size[0] ))
                                   * tile_size[0] );
  aligned_level_bbox.max().y() = ( int(ceilf( float(level_bbox.max().y()) / tile_size[1] ))
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
      results.push_back(loc);
        
      ++tile_x;
      dest_col += tile_size[0];
    }
    ++tile_y;
    dest_row += tile_size[1];
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

Vector2 TestPatternTileGenerator::minmax() { return Vector2(0.0, 1.0); }

PixelRGBA<float32> TestPatternTileGenerator::sample(int x, int y) {
  PixelRGBA<float32> result;
  return result;
}

int TestPatternTileGenerator::cols() const { return 2048; }
int TestPatternTileGenerator::rows() const { return 2048; }
PixelFormatEnum TestPatternTileGenerator::pixel_format() const { return VW_PIXEL_RGBA; }
ChannelTypeEnum TestPatternTileGenerator::channel_type() const { return VW_CHANNEL_FLOAT32; }
Vector2i TestPatternTileGenerator::tile_size() const { 
  return Vector2i(m_tile_size, m_tile_size); 
}
int32 TestPatternTileGenerator::num_levels() const {
  return 4;
}

// --------------------------------------------------------------------------
//                         PLATE FILE TILE GENERATOR
// --------------------------------------------------------------------------

PlatefileTileGenerator::PlatefileTileGenerator(std::string platefile_name) :
  m_platefile(new vw::platefile::PlateFile(platefile_name)) {
  std::cout << "\t--> Loading platefile \"" << platefile_name << "\" with " 
            << m_platefile->depth() << " levels.\n";
}


#define VW_DELEGATE_BY_PIXEL_TYPE(func, arg1, arg2)                                  \
  switch (this->pixel_format()) {                                                    \
  case VW_PIXEL_GRAY:                                                                \
  case VW_PIXEL_GRAYA:                                                               \
      if (this->channel_type() == VW_CHANNEL_UINT8) {                                \
        return func<PixelGrayA<uint8> >(arg1, arg2);                                 \
      } else if (this->channel_type() == VW_CHANNEL_FLOAT32) {                       \
        return func<PixelGrayA<float> >(arg1, arg2);                                 \
      } else {                                                                       \
        std::cout << "This platefile has a channel type that is not yet support by vwv.\n"; \
        std::cout << "Exiting...\n\n";                                               \
        exit(0);                                                                     \
      }                                                                              \
      break;                                                                         \
    case VW_PIXEL_RGB:                                                               \
    case VW_PIXEL_RGBA:                                                              \
      if (this->channel_type() == VW_CHANNEL_UINT8) {                                \
        return func<PixelRGBA<uint8> >(arg1, arg2);                                  \
      } else {                                                                       \
        std::cout << "This platefile has a channel type that is not yet support by vwv.\n"; \
        std::cout << "Exiting...\n\n";                                               \
        exit(0);                                                                     \
      }                                                                              \
      break;                                                                         \
    default:                                                                         \
      std::cout << "This platefile has a pixel format that is not yet support by vwv.\n"; \
      std::cout << "Exiting...\n\n";                                                 \
      exit(0);                                                                       \
    }                                                                                \

template<class PixelT>
Vector2 minmax_impl(TileLocator const& tile_info,
                    boost::shared_ptr<vw::platefile::PlateFile> platefile) {
  ImageView<PixelT> tile;
  platefile->read(tile, tile_info.col, tile_info.row, tile_info.level);
  typename PixelChannelType<PixelT>::type min, max;
  min_max_channel_values(alpha_to_mask(tile), min, max);
  Vector2 result(min, max);
  std::cout << "Here is the original answer: " << result << "\n";
  result /= ChannelRange<typename PixelChannelType<PixelT>::type>::max();
  std::cout << "NEW MIN AND MAX: " << result << "\n";
  return result;
}

Vector2 PlatefileTileGenerator::minmax() {
  try {
    TileLocator loc;
    loc.col = 0;
    loc.row = 0;
    loc.level = 0;
    VW_DELEGATE_BY_PIXEL_TYPE(minmax_impl, loc, m_platefile)
  } catch (platefile::TileNotFoundErr &e) {
    return Vector2(0,1);
  }
}

// template<class PixelT>
// std::string sample_impl(TileLocator const& tile_info,
//                         Vector2 const& px_loc) {
//   ImageView<PixelT> tile;
//   platefile->read(tile, tile_info.col, tile_info.row, tile_info.level);
//   VW_ASSERT(px_loc[0] >= 0 && px_loc[0] < tile.cols() &&
//             px_loc[1] >= 0 && px_loc[1] < tile.rows(),
//             ArgumentErr() << "sample_impl() invalid pixel location");
//   return tile(px_loc[0], px_loc[1]);
// }

PixelRGBA<float32> PlatefileTileGenerator::sample(int x, int y) {
  // TileLocator tile_loc;
  // tile_loc.col = floor(x/this->tile_size[0]);
  // tile_loc.row = floor(y/this->tile_size[1]);
  // tile_loc.level = this->num_levels();
  // px_loc = Vector2(x % this->tile_size[0], 
  //                  y % this->tile_size[1]);

  try {
    return PixelRGBA<float32>();
    //    VW_DELEGATE_BY_PIXEL_TYPE(sample_tile_impl, tile_loc, px_loc)
  } catch (platefile::TileNotFoundErr &e) {
    ImageView<PixelGrayA<uint8> > blank_tile(1,1);
    return PixelRGBA<float32>();
  }
}

template <class PixelT>
boost::shared_ptr<ViewImageResource> generate_tile_impl(TileLocator const& tile_info,
                                      boost::shared_ptr<vw::platefile::PlateFile> platefile) {
  ImageView<PixelT> tile;
  platefile->read(tile, tile_info.col, tile_info.row, tile_info.level);
  return boost::shared_ptr<ViewImageResource>( new ViewImageResource(tile) );
}

boost::shared_ptr<ViewImageResource> PlatefileTileGenerator::generate_tile(TileLocator const& tile_info) {

  vw_out(DebugMessage, "gui") << "Request to generate platefile tile " 
                              << tile_info.col << " " << tile_info.row 
                              << " @ " << tile_info.level << "\n";
  
  try {
    VW_DELEGATE_BY_PIXEL_TYPE(generate_tile_impl, tile_info, m_platefile)
  } catch (platefile::TileNotFoundErr &e) {
    ImageView<PixelGrayA<uint8> > blank_tile(1,1);
    return boost::shared_ptr<ViewImageResource>( new ViewImageResource(blank_tile) );    
  }
  
  vw_throw(NoImplErr() << "Unsupported pixel format or channel type in TileGenerator.\n");
}

int PlatefileTileGenerator::cols() const {
  return this->tile_size()[0] * pow(2, m_platefile->depth());
}

int PlatefileTileGenerator::rows() const {
  return this->tile_size()[1] * pow(2, m_platefile->depth());
}

PixelFormatEnum PlatefileTileGenerator::pixel_format() const {
  return m_platefile->pixel_format();
}

ChannelTypeEnum PlatefileTileGenerator::channel_type() const {
  return m_platefile->channel_type();
}

Vector2i PlatefileTileGenerator::tile_size() const {
  return Vector2i(m_platefile->default_tile_size(),
                  m_platefile->default_tile_size());
}

int32 PlatefileTileGenerator::num_levels() const {
  return m_platefile->depth();
}


// --------------------------------------------------------------------------
//                             IMAGE TILE GENERATOR
// --------------------------------------------------------------------------

ImageTileGenerator::ImageTileGenerator(std::string filename) : 
  m_filename(filename), m_rsrc( DiskImageResource::open(filename) ) {
  vw_out(0) << "\t--> Loading image: " << filename << ".\n";
}


// This little template makes the code below much cleaner.
template <class PixelT>
boost::shared_ptr<ViewImageResource> do_image_tilegen(boost::shared_ptr<ImageResource> rsrc,
                                                      BBox2i tile_bbox, 
                                                      int level, int max_levels) {
  ImageView<PixelT> tile(tile_bbox.width(), tile_bbox.height());
  rsrc->read(tile.buffer(), tile_bbox);
  ImageView<PixelT> reduced_tile = subsample(tile, pow(2,max_levels - level));
  return boost::shared_ptr<ViewImageResource>( new ViewImageResource(reduced_tile) );
}

boost::shared_ptr<ViewImageResource> ImageTileGenerator::generate_tile(TileLocator const& tile_info) {
  
  // Compute the bounding box of the image and the tile that is being
  // requested.  The bounding box of the tile depends on the pyramid
  // level we are looking at.
  BBox2i image_bbox(0,0,m_rsrc->cols(),m_rsrc->rows());
  BBox2i tile_bbox = tile_to_bbox(this->tile_size(), tile_info.col, 
                                  tile_info.row, tile_info.level, this->num_levels());

  // Check to make sure the image intersects the bounding box.  Print
  // an error to screen and return an empty tile if it does not.
  if (!image_bbox.intersects(tile_bbox)) {
    vw_out() << "WARNING in ImageTileGenerator: a tile was requested that doesn't exist.";
    ImageView<PixelGray<uint8> > blank_tile(this->tile_size()[0], this->tile_size()[1]);
    return boost::shared_ptr<ViewImageResource>( new ViewImageResource(blank_tile) );
  }

  // Make sure we don't access any pixels outside the image boundary
  // by cropping the tile to the image dimensions.
  tile_bbox.crop(image_bbox);

  switch (this->pixel_format()) {
  case VW_PIXEL_GRAY:
    if (this->channel_type() == VW_CHANNEL_UINT8) {
      return do_image_tilegen<PixelGray<uint8> >(m_rsrc, tile_bbox, 
                                                 tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_INT16) {
      return do_image_tilegen<PixelGray<int16> >(m_rsrc, tile_bbox, 
                                                 tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_UINT16) {
      return do_image_tilegen<PixelGray<uint16> >(m_rsrc, tile_bbox, 
                                                  tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_FLOAT32) {
      return do_image_tilegen<PixelGray<float> >(m_rsrc, tile_bbox, 
                                                  tile_info.level, this->num_levels());
    } else {
      std::cout << "This platefile has a channel type that is not yet support by vwv.\n";
      std::cout << "Exiting...\n\n";
      exit(0);
  }    
  break;

  case VW_PIXEL_GRAYA:
    if (this->channel_type() == VW_CHANNEL_UINT8) {
      return do_image_tilegen<PixelGrayA<uint8> >(m_rsrc, tile_bbox, 
                                                  tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_INT16) {
      return do_image_tilegen<PixelGrayA<int16> >(m_rsrc, tile_bbox, 
                                                  tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_UINT16) {
      return do_image_tilegen<PixelGrayA<uint16> >(m_rsrc, tile_bbox, 
                                                   tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_FLOAT32) {
      return do_image_tilegen<PixelGrayA<float> >(m_rsrc, tile_bbox, 
                                                  tile_info.level, this->num_levels());
    } else {
      std::cout << "This platefile has a channel type that is not yet support by vwv.\n";
      std::cout << "Exiting...\n\n";
      exit(0);
      }
      
    break;

  case VW_PIXEL_RGB:
    if (this->channel_type() == VW_CHANNEL_UINT8) {
      return do_image_tilegen<PixelRGB<uint8> >(m_rsrc, tile_bbox, 
                                                tile_info.level, this->num_levels());
    } else {
      std::cout << "This platefile has a channel type that is not yet support by vwv.\n";
      std::cout << "Exiting...\n\n";
      exit(0);
    }
      
    break;

  case VW_PIXEL_RGBA:
    if (this->channel_type() == VW_CHANNEL_UINT8) {
      return do_image_tilegen<PixelRGBA<uint8> >(m_rsrc, tile_bbox, 
                                                 tile_info.level, this->num_levels());
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

Vector2 ImageTileGenerator::minmax() {
  return Vector2(0,1.0); // TODO: Implement this properly
}

PixelRGBA<float32> ImageTileGenerator::sample(int x, int y) {
  PixelRGBA<float32> result; // TODO: Implement this properly
  return result;
}


int ImageTileGenerator::cols() const {
  return m_rsrc->cols();
}

int ImageTileGenerator::rows() const {
  return m_rsrc->rows();
}

PixelFormatEnum ImageTileGenerator::pixel_format() const {
  return m_rsrc->pixel_format();
}

ChannelTypeEnum ImageTileGenerator::channel_type() const {
  return m_rsrc->channel_type();
}

Vector2i ImageTileGenerator::tile_size() const {
  return m_rsrc->block_size();
}

int32 ImageTileGenerator::num_levels() const {
  int32 max_dimension = std::max(this->cols(), this->rows());
  int32 max_tilesize = std::max(this->tile_size()[0], this->tile_size()[1]);
  return ceil(log(float(max_dimension) / max_tilesize) / log(2));
}
