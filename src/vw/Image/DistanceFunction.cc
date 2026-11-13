// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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


/// \file DistanceFunction.cc

#include <vw/Image/DistanceFunction.h>

namespace vw {

  // Find the Euclidean distance function to boundary of valid pixels
  // of an image. An invalid pixel equals 0.
  // - Pixels at the edges of the image are considered invalid.
  // - Invalid pixels have a distance of 0. The distance is
  //   positive for all valid pixels, increasing when moving away from
  //   invalid pixels, until the magnitude reaches max_dist where it
  //   stops growing.
  // - The complexity of this function is length of boundary times
  //   max_dist * max_dist, so it can be slow.
  // - It is assumed that the output image fits fully in memory.

void bounded_dist(ImageViewRef<int> image, double max_dist, ImageView<double> & dist) {
  
  if (max_dist < 0.0) 
    vw_throw(ArgumentErr() << "Expecting positive max_dist.");
  
  int max_dist_int = ceil(max_dist); // an int overestimate
  double max_dist_sq = max_dist * max_dist;
  int cols = image.cols(), rows = image.rows();
  dist.set_size(cols, rows);

  for (int col = 0; col < cols; col++) {
    for (int row = 0; row < rows; row++) {
      // Initialize to the max for valid pixels and not at image edges.
      if (image(col, row) != 0 && col > 0 && col < cols - 1 && row > 0 && row < rows - 1)
        dist(col, row) = max_dist;
      else
        dist(col, row) = 0;
    }
  }

  // For each boundary pixel (col, row), which is an invalid pixel
  // with a valid neighbors, adjust the values of dist which are
  // affected by it.
  for (int col = 0; col < cols; col++) {
    for (int row = 0; row < rows; row++) {

      // Look for a boundary pixel
      if (dist(col, row) > 0)
        continue; // valid, but need invalid
      bool is_bd_pix = false;
      for (int c = col - 1; c <= col + 1; c++) {
        for (int r = row - 1; r <= row + 1; r++) {
          if (c < 0 || c >= cols || r < 0 || r >= rows)
            continue;
          if (dist(c, r) > 0) {
            is_bd_pix = true; // has valid neighbor
            break;
          }
        }
        if (is_bd_pix) 
          break;
      }
      if (!is_bd_pix) 
        continue; // does not have valid neighbors

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

          if (dist(c, r) > 0) {
            double d = sqrt(dsq);
            // If this was positive before, it stays positive, as d > 0, since
            // it is the distance between a valid and an invalid pixel
            dist(c, r) = std::min(dist(c, r), d);
          }
        }
      }
        
    }
  }

  return;
}
   
// Clamped signed distance to the boundary between a positive (lit) region and
// its complement (shadow). Here it is assumed that there is a region in which
// lit_grass_dist is positive, while in its complement the shadow_grass_dist is
// positive. The boundary is where both are no more than 1 in value.
ImageView<double> clamped_blend_dist(ImageView<float> const& lit_grass_dist,
                                     ImageView<float> const& shadow_grass_dist,
                                     double lit_blend_length,
                                     double shadow_blend_length) {

  if (lit_grass_dist.cols() != shadow_grass_dist.cols() ||
      lit_grass_dist.rows() != shadow_grass_dist.rows())
    vw_throw(ArgumentErr()
             << "Input images to clamped_blend_dist must have the same dimensions.");

  ImageView<double> dist_to_bd;
  dist_to_bd.set_size(lit_grass_dist.cols(), lit_grass_dist.rows());
  for (int col = 0; col < lit_grass_dist.cols(); col++) {
    for (int row = 0; row < lit_grass_dist.rows(); row++) {

      if (lit_grass_dist(col, row) > 1.5 * lit_blend_length) {
        // Too far in the lit region
        dist_to_bd(col, row) = lit_blend_length; // clamp at the blending length
        continue;
      }

      if (shadow_grass_dist(col, row) > 1.5 * shadow_blend_length) {
        // Too far in the shadow region
        dist_to_bd(col, row) = -shadow_blend_length;
        continue;
      }

      // Find the shortest Euclidean distance to the no-data region.
      double max_dist = std::max(lit_blend_length, shadow_blend_length);
      double signed_dist = 0.0;
      if (lit_grass_dist(col, row) > 0) {
        signed_dist = lit_blend_length;
      } else if (shadow_grass_dist(col, row) > 0) {
        signed_dist = -shadow_blend_length;
      }

      for (int col2 = std::max(0.0, col - max_dist);
           col2 <= std::min(lit_grass_dist.cols() - 1.0, col + max_dist);
           col2++) {

        // Estimate the range of rows for the circle with given radius
        // at given col value.
        double ht_val = ceil(sqrt(double(max_dist * max_dist) -
                                  double((col - col2) * (col - col2))));

        for (int row2 = std::max(0.0, row - ht_val);
             row2 <= std::min(lit_grass_dist.rows() - 1.0, row + ht_val);
             row2++) {

          if (lit_grass_dist(col2, row2) > 1 || shadow_grass_dist(col2, row2) > 1)
            continue; // not at the boundary

          // See if the current point is closer than anything so far
          double curr_dist = sqrt(double(col - col2) * (col - col2) +
                                  double(row - row2) * (row - row2));
          if (lit_grass_dist(col, row) > 0) {
            if (curr_dist < signed_dist)
              signed_dist = curr_dist;
          } else if (shadow_grass_dist(col, row) > 0) {
            if (curr_dist < -signed_dist)
              signed_dist = -curr_dist;
          }

        }
      }

      // The closest we've got
      dist_to_bd(col, row) = signed_dist;
    }
  }
  return dist_to_bd;
}

} //end namespace vw
