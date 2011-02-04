// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_GUI_TESTPATTERNTILEGENERATOR_H__
#define __VW_GUI_TESTPATTERNTILEGENERATOR_H__

#include <vw/gui/TileGenerator.h>

namespace vw {
namespace gui {

  class TestPatternTileGenerator : public TileGenerator {
    int m_tile_size;

  public:
    TestPatternTileGenerator(int tile_size) : m_tile_size(tile_size) {}
    virtual ~TestPatternTileGenerator() {}

    virtual boost::shared_ptr<SrcImageResource> generate_tile(TileLocator const& tile_info);
    virtual Vector2 minmax();
    virtual PixelRGBA<float> sample(int x, int y, int level, int transaction_id);

    virtual int cols() const;
    virtual int rows() const;
    virtual PixelFormatEnum pixel_format() const;
    virtual ChannelTypeEnum channel_type() const;
    virtual Vector2i tile_size() const;
    virtual int32 num_levels() const;
  };

}} // namespace vw::gui

#endif

