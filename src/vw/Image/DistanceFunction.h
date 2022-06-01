// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file DistanceFunction.h
///

#ifndef __VW_IMAGE_BOUNDEDSIGNEDDIST_H__
#define __VW_IMAGE_BOUNDEDSIGNEDDIST_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/MaskViews.h>

namespace vw {

  // Find the Euclidean distance function to boundary of valid pixels
  // of an image of Vector3 pixels. An invalid pixel equals Vector3().
  // - Pixels at the edges of the image are considered invalid.
  // - Invalid pixels have a distance of 0. The distance is
  //   positive for all valid pixels, increasing when moving away from
  //   invalid pixels, until the magnitude reaches max_dist where it
  //   stops growing.
  // - The complexity of this function is length of boundary times
  //   max_dist * max_dist, so it can be slow.
  // - It is assumed that the output image fits fully in memory.
  //   TODO(oalexan1): Make this more generic.
  void bounded_boundary_dist(ImageViewRef<Vector3> image, double max_dist,
                             ImageView<double> & dist);
  
  // Find the signed Euclidean distance function to boundary of invalid
  // pixels. Invalid pixels neighboring valid pixels have a distance of
  // 0. The distance is negative for other invalid pixels, and is
  // positive for all valid pixels, with both increasing in magnitude
  // away from the boundary, till the magnitude reaches current value
  // where it stops growing.
  // This is a single-threaded function which does all the calculation
  // in memory. The complexity is proportional to
  // boundary_len * max_dist * max_dist.
  template<class PixelT>
  void bounded_signed_dist(ImageViewRef<PixelMask<PixelT>> image,
                           double max_dist, ImageView<double> & dist) {
    
    if (max_dist < 0.0) 
      vw_throw(ArgumentErr() << "Expecting positive max_dist.");

    int max_dist_int = ceil(max_dist); // an int overestimate
    double max_dist_sq = max_dist * max_dist;
    int cols = image.cols(), rows = image.rows();
    dist.set_size(cols, rows);

    for (int col = 0; col < cols; col++) {
      for (int row = 0; row < rows; row++) {
        // Initialize to the max, with the right sign
        if (is_valid(image(col, row)))
          dist(col, row) = max_dist;
        else
          dist(col, row) = -max_dist;
      }
    }

    // For each boundary pixel (col, row), adjust the values of dist
    // which are affected by it.
    for (int col = 0; col < cols; col++) {
      for (int row = 0; row < rows; row++) {

        // Look for a boundary pixel, which is an invalid pixel with valid neighbors
        if (is_valid(image(col, row)))
          continue; // valid
        bool is_bd_pix = false;
        for (int c = col - 1; c <= col + 1; c++) {
          for (int r = row - 1; r <= row + 1; r++) {
            if (c < 0 || c >= cols || r < 0 || r >= rows)
              continue;
            if (is_valid(image(c, r))) {
              is_bd_pix = true;
              break;
            }
          }
          if (is_bd_pix) 
            break;
        }
        if (!is_bd_pix) 
          continue;

        // Found a boundary pixel. Pixels closer to it than max_dist
        // may have their distance function modified.
        for (int c = col - max_dist_int; c <= col + max_dist_int; c++) {
          for (int r = row - max_dist_int; r <= row + max_dist_int; r++) {
            if (c < 0 || c >= cols || r < 0 || r >= rows)
              continue;

            // Cast to double before multiplying to avoid integer overflow
            double dsq = double(c - col) * double(c - col) + 
              double(r - row) * double(r - row);

            // Too far 
            if (dsq >= max_dist_sq) 
              continue;

            double d = sqrt(dsq);
            if (is_valid(image(c, r))) 
              dist(c, r) = std::min(dist(c, r), d);
            else
              dist(c, r) = std::max(dist(c, r), -d);
          }
        }
        
      }
    }

    return;
  }

} //end namespace vw

#endif//__VW_IMAGE_BOUNDEDSIGNEDDIST_H__
