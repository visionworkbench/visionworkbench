// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifndef __VW_PLATE_TILEMANIPULATION_H__
#define __VW_PLATE_TILEMANIPULATION_H__

#include <vw/Image/Transform.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Interpolation.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Filter.h>
#include <vw/Math/BBox.h>
#include <list>

namespace vw {
namespace platefile {

  class PlateFile;

  // Given a bbox, returns a list of smaller bboxes that perfectly
  // tile the space of the larger bbox.
  std::list<vw::BBox2i> bbox_tiles(vw::BBox2i const& bbox, int width, int height);

  template <class PixelT>
  void mipmap_one_tile(ImageView<PixelT>& dest, uint32 tile_size, const ImageView<PixelT>& UL, const ImageView<PixelT>& UR, const ImageView<PixelT>& LL, const ImageView<PixelT>& LR, bool blur = true)
  {
    VW_ASSERT(!UL || (UL.cols() == int32(tile_size) && UL.rows() == int32(tile_size)), LogicErr() << "Tiles must be the same size as tile_size");
    VW_ASSERT(!UR || (UR.cols() == int32(tile_size) && UR.rows() == int32(tile_size)), LogicErr() << "Tiles must be the same size as tile_size");
    VW_ASSERT(!LL || (LL.cols() == int32(tile_size) && LL.rows() == int32(tile_size)), LogicErr() << "Tiles must be the same size as tile_size");
    VW_ASSERT(!LR || (LR.cols() == int32(tile_size) && LR.rows() == int32(tile_size)), LogicErr() << "Tiles must be the same size as tile_size");
    VW_ASSERT(UL || UR || LL || LR, LogicErr() << "Must compose at least one tile");

    ImageView<PixelT> super(2*tile_size, 2*tile_size);

    if (UL) crop(super, 0,         0,         tile_size, tile_size) = UL;
    if (UR) crop(super, tile_size, 0,         tile_size, tile_size) = UR;
    if (LL) crop(super, 0,         tile_size, tile_size, tile_size) = LL;
    if (LR) crop(super, tile_size, tile_size, tile_size, tile_size) = LR;

    // We subsample after blurring with a standard 2x2 box filter.
    std::vector<float> kernel(2);
    kernel[0] = kernel[1] = 0.5;

    if (blur)
      dest = subsample( separable_convolution_filter( super, kernel, kernel, 1, 1, ConstantEdgeExtension() ), 2);
    else
      dest = subsample( super, 2 );
  }

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

}} // vw::plate

#endif
