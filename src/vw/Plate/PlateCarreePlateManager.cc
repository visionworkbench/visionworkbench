// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/PlateCarreePlateManager.h>
using namespace vw::platefile;
using namespace vw;

template <class PixelT> std::vector<TileInfo>
PlateCarreePlateManager<PixelT>::kml_image_tiles( BBox2i const& image_bbox,
                                  int32 const /*resolution*/,
                                  int32 const tile_size) {
  std::vector<TileInfo> result;

  // There's no point in starting the search before there is good
  // image data, so we adjust the start point here.
  int32 minx = int(floor(image_bbox.min().x() / tile_size) * tile_size);
  int32 miny = int(floor(image_bbox.min().y() / tile_size) * tile_size);
  int x = minx / tile_size;
  int y = miny / tile_size;

  // Iterate over the bounding boxes in the entire KML space...
  int curx = minx;
  int cury = miny;
  while (cury < image_bbox.max().y()) {
    while (curx <= image_bbox.max().x()) {

      TileInfo be(x, y, BBox2i(curx, cury, tile_size, tile_size));

      // ...but only add bounding boxes that overlap with the image.
      if (image_bbox.intersects(be.bbox))
        result.push_back(be);

      curx += tile_size;
      ++x;
    }
    curx = minx;
    x = minx / tile_size;
    cury += tile_size;
    ++y;
  }
  return result;
}

/// This function generates a specific mipmap tile at the given
/// col, row, and level, and transaction_id.
template <class PixelT>
void PlateCarreePlateManager<PixelT>::generate_mipmap_tile(int col,
                                                   int row,
                                                   int level,
                                                   int transaction_id) const {

  // Create an image large enough to store all of the child nodes
  int tile_size = m_platefile->default_tile_size();
  ImageView<PixelT> super(2*tile_size, 2*tile_size);

  // Iterate over the children, gathering them and (recursively)
  // regenerating them if necessary.
  for( int j=0; j<2; ++j ) {
    for( int i=0; i<2; ++i ) {
      try {
        int child_col = 2*col+i;
        int child_row = 2*row+j;
        vw_out(VerboseDebugMessage, "platefile") << "Reading tile "
                                                 << child_col << " " << child_row
                                                 << " @  " << (level+1) << "\n";
        ImageView<PixelT> child;
        m_platefile->read(child, child_col, child_row, level+1,
                          transaction_id, true); // exact_transaction_match == true
        crop(super,tile_size*i,tile_size*j,tile_size,tile_size) = child;
      } catch (TileNotFoundErr &e) {
        // If that fails, then there is no tile.  Do nothing to this quadrant of super.
      }
    }
  }

  // We subsample after blurring with a standard 2x2 box filter.
  std::vector<float> kernel(2);
  kernel[0] = kernel[1] = 0.5;

  ImageView<PixelT> new_tile = subsample( separable_convolution_filter( super,
                                                                        kernel,
                                                                        kernel,
                                                                        1, 1,
                                                                        ConstantEdgeExtension() ), 2);

  if (!is_transparent(new_tile)) {
    vw_out(VerboseDebugMessage, "platefile") << "Writing " << col << " " << row
                                             << " @ " << level << "\n";
    m_platefile->write_update(new_tile, col, row, level, transaction_id);
  }
}

// Explicit template instatiation
namespace vw {
namespace platefile {

#define VW_INSTANTIATE_TOAST_PLATEMANAGER_TYPES(PIXELT)                 \
  template std::vector<TileInfo>                                        \
  PlateCarreePlateManager<PIXELT >::kml_image_tiles( BBox2i const& image_bbox,  \
                                             int32 const resolution,    \
                                             int32 const tile_size);    \
  template void                                                         \
  PlateCarreePlateManager<PIXELT >::generate_mipmap_tile(int col,               \
                                                 int row,               \
                                                 int level,             \
                                                 int transaction_id) const; \

  VW_INSTANTIATE_TOAST_PLATEMANAGER_TYPES(PixelGrayA<uint8>)
  VW_INSTANTIATE_TOAST_PLATEMANAGER_TYPES(PixelGrayA<int16>)
  VW_INSTANTIATE_TOAST_PLATEMANAGER_TYPES(PixelGrayA<float32>)
  VW_INSTANTIATE_TOAST_PLATEMANAGER_TYPES(PixelRGBA<uint8>)
}}


