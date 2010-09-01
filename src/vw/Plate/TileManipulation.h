// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_TILEMANIPULATION_H__
#define __VW_PLATE_TILEMANIPULATION_H__

#include <vw/Plate/PlateFile.h>
#include <vw/Image/Transform.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Interpolation.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Algorithms.h>
#include <vw/Math/BBox.h>
#include <list>

namespace vw {
namespace platefile {

  // Given a bbox, returns a list of smaller bboxes that perfectly
  // tile the space of the larger bbox.
  std::list<vw::BBox2i> bbox_tiles(vw::BBox2i const& bbox, int width, int height);

  // Resample image by reaching up a few levels and using the data there.
  template <typename PixelT>
  ImageView<PixelT> resample_img_from_level(const ImageView<PixelT> &src_tile, int src_col, int src_row,
                                            int src_level, int dst_col, int dst_row, int dst_level)
  {
    int level_diff = dst_level - src_level;
    VW_ASSERT(level_diff > 0, LogicErr() << "resample_img_from_level: dst_level must be > src_level");

    int scaling_factor = 1 << level_diff;

    int subtile_u = dst_col - src_col * scaling_factor;
    int subtile_v = dst_row - src_row * scaling_factor;
    BBox2i subtile_bbox(src_tile.cols() * subtile_u,
        src_tile.rows() * subtile_v,
        src_tile.cols(), src_tile.rows());

    // Scale up and interpolate the src_tile, then crop out the
    // subtile that we need.
    // XXX: ConstantEdgeExtension? It shouldn't ever go outside the box. Maybe NoEdgeExtension?
      return crop(
          transform(src_tile, ResampleTransform(scaling_factor, scaling_factor),
                    ConstantEdgeExtension(), BilinearInterpolation()),
          subtile_bbox);
  }

  // Create a tile filled with "value" at the given location
  // DOES NOT ISSUE A write_request FIRST! That's the caller's job.
  template <typename PixelT>
  void create_uniform_tile(PixelT value, PlateFile& output, int col, int row, int level, int transaction_id) {
    ImageView<PixelT> img(output.default_tile_size(), output.default_tile_size());
    vw::fill( img, value );
    output.write_update(img, col, row, level, transaction_id);
  }

}} // vw::plate

#endif
