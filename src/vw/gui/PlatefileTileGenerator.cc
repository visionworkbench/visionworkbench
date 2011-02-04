// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/gui/PlatefileTileGenerator.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/ViewImageResource.h>

namespace vw { namespace gui {

PlatefileTileGenerator::PlatefileTileGenerator(const std::string& platefile_name) :
  m_platefile(new vw::platefile::PlateFile(platefile_name)) {
  m_num_levels = m_platefile->num_levels();
  std::cout << "\t--> Loading platefile \"" << platefile_name << "\" with " << m_num_levels << " levels.\n";
}

#define VW_DELEGATE_BY_PIXEL_TYPE(func, arg1, arg2)                                  \
  switch (this->pixel_format()) {                                                    \
  case VW_PIXEL_GRAY:                                                                \
  case VW_PIXEL_GRAYA:                                                               \
      if (this->channel_type() == VW_CHANNEL_UINT8) {                                \
        return func<PixelGrayA<uint8> >(arg1, arg2);                                 \
      } else if (this->channel_type() == VW_CHANNEL_INT16) {                         \
        return func<PixelGrayA<int16> >(arg1, arg2);                                 \
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
      } else if (this->channel_type() == VW_CHANNEL_UINT16) {                        \
        return func<PixelRGBA<uint16> >(arg1, arg2);                                 \
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
  platefile->read(tile, tile_info.col, tile_info.row, tile_info.level, tile_info.transaction_id);
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
//   platefile->read(tile, tile_info.col, tile_info.row, tile_info.level, tile_info.transaction_id);
//   VW_ASSERT(px_loc[0] >= 0 && px_loc[0] < tile.cols() &&
//             px_loc[1] >= 0 && px_loc[1] < tile.rows(),
//             ArgumentErr() << "sample_impl() invalid pixel location");
//   return tile(px_loc[0], px_loc[1]);
// }

PixelRGBA<float32> PlatefileTileGenerator::sample(int /*x*/, int /*y*/, int /*level*/, int /*transaction_id*/) {
  // TileLocator tile_loc;
  // tile_loc.col = floor(x/this->tile_size[0]);
  // tile_loc.row = floor(y/this->tile_size[1]);
  // tile_loc.level = this->num_levels();
  // px_loc = Vector2(x % this->tile_size[0],
  //                  y % this->tile_size[1]);

  try {
    return PixelRGBA<float32>(1.0, 0.0, 0.0, 1.0);
    //    VW_DELEGATE_BY_PIXEL_TYPE(sample_tile_impl, tile_loc, px_loc)
  } catch (platefile::TileNotFoundErr &e) {
    ImageView<PixelGrayA<uint8> > blank_tile(1,1);
    return PixelRGBA<float32>();
  }
}

template <class PixelT>
boost::shared_ptr<SrcImageResource> generate_tile_impl(TileLocator const& tile_info,
                                      boost::shared_ptr<vw::platefile::PlateFile> platefile) {
  ImageView<PixelT> tile(1,1);
  try {

    platefile->read(tile, tile_info.col, tile_info.row,
                    tile_info.level, tile_info.transaction_id,
                    tile_info.exact_transaction_id_match);

  } catch (platefile::TileNotFoundErr &e) {
    return boost::shared_ptr<SrcImageResource>(make_point_src(PixelRGBA<uint8>(0, 20, 0, 255)));
  } catch (vw::IOErr &e) {
    std::cout << "WARNING: AMQP ERROR -- " << e.what() << "\n";
  }
  return boost::shared_ptr<SrcImageResource>( new ViewImageResource(tile) );
}

boost::shared_ptr<SrcImageResource> PlatefileTileGenerator::generate_tile(TileLocator const& tile_info) {

  vw_out(DebugMessage, "gui") << "Request to generate platefile tile "
                              << tile_info.col << " " << tile_info.row
                              << " @ " << tile_info.level << "\n";

  VW_DELEGATE_BY_PIXEL_TYPE(generate_tile_impl, tile_info, m_platefile)

  // If we get to here, then there was no support for the pixel format.
  vw_throw(NoImplErr() << "Unsupported pixel format or channel type in TileGenerator.\n");
}

int PlatefileTileGenerator::cols() const {
  return this->tile_size()[0] * (1 << (m_num_levels-1));
}

int PlatefileTileGenerator::rows() const {
  return this->tile_size()[1] * (1 << (m_num_levels-1));
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
  return m_num_levels;
}

}} // namespace vw::gui
