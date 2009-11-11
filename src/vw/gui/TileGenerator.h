// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_GUI_TILEGENERATOR_H__
#define __VW_GUI_TILEGENERATOR_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/BBox.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>

namespace vw {
namespace gui {

  struct TileLocator {
    int col;
    int row;
    int level;

    bool is_valid() {
      return (col >= 0 && row >= 0 && col < pow(2, level) && row < pow(2, level));
    }
  };

  BBox2i tile_to_bbox(int32 tile_size, int col, int row, int level, int max_level);
  
  std::list<TileLocator> bbox_to_tiles(int32 tile_size, BBox2i bbox, int level, int max_level);

  // --------------------------------------------------------------------------
  //                              TILE GENERATOR
  // --------------------------------------------------------------------------
  class TileGenerator {
  public:
    virtual ~TileGenerator() {}
    virtual ImageView<PixelRGBA<float> > generate_tile(TileLocator const& tile_info) = 0;

    static boost::shared_ptr<TileGenerator> create(std::string filename) {
      // If ends in .plate, then assume platefile.
      
      // Otherwise, assume an image.
    }
  };

  class DummyTileGenerator {
    int m_tile_size;

  public:
    DummyTileGenerator(int tile_size) : m_tile_size(tile_size) {}
    virtual ~DummyTileGenerator() {}

    virtual ImageView<PixelRGBA<float> > generate_tile(TileLocator const& tile_info);
  };

  // class PlatefileTileGenerator {
  //   boost::shared_ptr<PlateFile> m_platefile;

  // public:
  //   PlatefileTileGenerator(std::string platefile_name) :
  //     m_platefile(new PlateFile(platefile_name) {}
  //   virtual ~PlatefileTileGenerator() {}

  //   virtual ImageView<PixelRGBA<float> > generate_tile(TileLocator const& index) { 

  //   }
  // };


}} // namespace vw::gui

#endif // __VW_GUI_TILEGENERATOR_H__

