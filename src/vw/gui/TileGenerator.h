// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_GUI_TILEGENERATOR_H__
#define __VW_GUI_TILEGENERATOR_H__

#include <vw/Image/ViewImageResource.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/ImageView.h>
#include <vw/Math/BBox.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/FundamentalTypes.h>

namespace vw {
namespace gui {

  struct TileLocator {
    int col;
    int row;
    int level;
    int transaction_id;
    bool exact_transaction_id_match;

    bool is_valid() const {
      return col >= 0 && row >= 0 && col < (1 << level) && row < (1 << level);
    }
  };

  // Given a tile index, return the bounding box of that tile coverage
  // in the bottom (i.e. highest resolution) level of the source image
  // pyramid.
  BBox2i tile_to_bbox(Vector2i tile_size, int col, int row, int level, int max_level);

  std::list<TileLocator> bbox_to_tiles(Vector2i tile_size, BBox2i bbox, int level,
                                       int max_level, int transaction_id,
                                       bool exact_transaction_id_match);

  // --------------------------------------------------------------------------
  //                              TILE GENERATOR
  // --------------------------------------------------------------------------

  class TileGenerator {
  public:
    virtual ~TileGenerator() {}
    virtual boost::shared_ptr<ViewImageResource> generate_tile(TileLocator const& tile_info) = 0;
    virtual Vector2 minmax() = 0;
    virtual PixelRGBA<float> sample(int x, int y, int level, int transaction_id) = 0;

    virtual int cols() const = 0;
    virtual int rows() const = 0;
    virtual PixelFormatEnum pixel_format() const = 0;
    virtual ChannelTypeEnum channel_type() const = 0;
    virtual Vector2i tile_size() const = 0;
    virtual int32 num_levels() const = 0;

    // Use this method to generate the correct type of TileGenerator
    // for a given filename.
    static boost::shared_ptr<TileGenerator> create(std::string filename);
  };

}} // namespace vw::gui

#endif

