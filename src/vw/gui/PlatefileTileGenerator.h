// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_GUI_PLATEFILETILEGENERATOR_H__
#define __VW_GUI_PLATEFILETILEGENERATOR_H__

#include <vw/gui/TileGenerator.h>

namespace vw {

  namespace platefile { class PlateFile; }

namespace gui {

  class PlatefileTileGenerator : public TileGenerator {
    boost::shared_ptr<vw::platefile::PlateFile> m_platefile;
    int m_num_levels;

  public:
    PlatefileTileGenerator(const std::string& platefile_name);
    virtual ~PlatefileTileGenerator() {}

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

