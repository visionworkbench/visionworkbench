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

#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Core/Exception.h>

namespace vw {

  // *******************************************************************
  // subdivide_bbox()
  // *******************************************************************
  // See the .h file for documentation.
  std::vector<BBox2i>
  subdivide_bbox(BBox2i const& object, int32 block_width, int32 block_height,
                 bool include_partials, bool full_size) {

    // This check makes the logic below unambiguous.
    if (!include_partials && full_size)
      vw::vw_throw(vw::LogicErr() << "subdivide_bbox: Cannot have include_partials be false and full_size be true.");

    std::vector<BBox2i> bboxes;

    for (int j_offset = 0; j_offset < object.height(); j_offset += block_height) {
      int32 j_dim = block_height;
      if ((object.height() - j_offset) < block_height) {
        if (!include_partials)
          continue; // Skip row of non-full size boxes.
        else
          j_dim = object.height() - j_offset;

        if (full_size) {
            // Grow the box inwards to make it full size. This will make it
            // overlap with the previous boxes.
            j_dim = block_height;
            j_offset = int32(object.height()) - block_height;
        }
      }
      
      for (int i_offset = 0; i_offset < object.width(); i_offset += block_width) {
        int32 i_dim = block_width;
        if ((object.width() - i_offset) < block_width) {
          if (!include_partials)
            continue; // Skip non-full size boxes.
          else
            i_dim = object.width() - i_offset;

          if (full_size) {
            // Grow the box inwards to make it full size. This will make it
            // overlap with the previous boxes.
            i_dim = block_width;
            i_offset = int32(object.width()) - block_width;
          }
        }
        
        bboxes.push_back(BBox2i(i_offset + object.min().x(),
                                j_offset + object.min().y(),
                                i_dim,j_dim));
      }
    }
    return bboxes;
  }

/// A weight at a given pixel, based on an image row. Return
/// zero where image values are not valid, and positive where valid.
/// - hCenterLine contains the center column at each row/col
/// - hMaxDistArray contains the width of the column at each row/col
double compute_line_weights(Vector2 const& pix, bool horizontal,
                            std::vector<double> const& centers,
                            std::vector<double> const& widths) {
  
  int primary_axis = 0, secondary_axis = 1; // Vertical
  if (horizontal) {
    primary_axis   = 1;
    secondary_axis = 0;
  }
  
  // We round below, to avoid issues when we are within numerical value
  // to an integer value for row/col.
  // To do: Need to do interpolation here.

  int pos = (int)round(pix[primary_axis]); // The row or column
  if (pos < 0                          ||
      pos >= (int)widths.size() ||
      pos >= (int)centers.size() )
    return 0;
  
  double max_dist = widths[pos]/2.0; // Half column width
  double center   = centers[pos];
  double dist     = fabs(pix[secondary_axis]-center); // Pixel distance from center column

  if (max_dist <= 0 || dist < 0)
    return 0;

  // We want to make sure the weight is positive (even if small) at
  // the first/last valid pixel.
  double tol = 1e-8*max_dist;
  
  // The weight is just a fraction of the distance from the centerline.
  double weight = std::max(double(0.0), (max_dist - dist + tol)/max_dist);

  return weight;
}

}

