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


/// \file InpaintView.cc

#include <vw/Image/InpaintView.h>

namespace vw {

// Fill nodata values in the input image using a search radius and a fill percentage.
// This assumes that the image fits fully in memory.
ImageView<PixelMask<double>> 
fillNodataWithSearchRadius(ImageView<PixelMask<double>> const& image,
                         double search_radius, double fill_power,
                         double fill_percent) { // percent is between 0 and 100.

  int search_radius_int = ceil(search_radius); // an int overestimate
  double search_radius_sq = search_radius * search_radius;
  int cols = image.cols(), rows = image.rows();
  
  // Initialize the output image as invalid
  ImageView<PixelMask<double>> out;
  out.set_size(cols, rows);
  for (int col = 0; col < cols; col++) {
    for (int row = 0; row < rows; row++) {
      out(col, row) = PixelMask<double>();
      out(col, row).invalidate();
    }
  }

  // For each pixel in the input image, if it is valid, copy it to the output.
  // Otherwise, fill from neighbors.
  for (int col = 0; col < cols; col++) {
    for (int row = 0; row < rows; row++) {

      if (is_valid(image(col, row))) {
        out(col, row) = image(col, row);
        continue;
      }

      double val = 0.0, weight_sum = 0.0, valid_count = 0.0, total_count = 0.0; 
      for (int c = col - search_radius_int; c <= col + search_radius_int; c++) {
        for (int r = row - search_radius_int; r <= row + search_radius_int; r++) {
          if (c < 0 || c >= cols || r < 0 || r >= rows)
            continue;

          // Cast to double before multiplying to avoid integer overflow
          double dsq = double(c - col) * double(c - col) + 
            double(r - row) * double(r - row);
   
          if (dsq > search_radius_sq)
            continue; // too far

          // We are within the search radius, so update the count
          total_count += 1.0;

          if (!is_valid(image(c, r)))
            continue;

          double weight = 1.0 / (1.0 + pow(dsq, fill_power/2.0));
          val += weight * image(c, r).child();
          weight_sum += weight;
          valid_count += 1.0;
        }
      }

      if (total_count == 0.0 || weight_sum == 0.0) 
        continue; // no valid pixels found

      // Note that both valid_count and total_count are double, so the division
      // is not an integer division.
      if (valid_count / total_count < fill_percent / 100.0)
        continue; // not enough valid pixels found

      out(col, row) = PixelMask<double>(val / weight_sum);
      out(col, row).validate();
    }
  }

  return out;                                  
}

// Fill nodata values in the input image using a search radius and a fill percentage.
// This assumes that the image fits fully in memory. Use multiple passes.
ImageView<PixelMask<double>> 
fillNodataWithSearchRadius(ImageView<PixelMask<double>> const& image,
                         double search_radius, double fill_power,
                         double fill_percent, int num_passes) { 
  ImageView<PixelMask<double>> out = image;
  for (int pass = 0; pass < num_passes; pass++) {
    out = fillNodataWithSearchRadius(out, search_radius, fill_power, fill_percent);
  }

  return out;
}
   
} //end namespace vw
