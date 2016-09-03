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

#include <vw/Image/ImageView.h>

// TODO: Use SSE to accelerate this operation
//#if defined(VW_ENABLE_SSE) && (VW_ENABLE_SSE==1)
//  #include <emmintrin.h>
//  #include <smmintrin.h> // SSE4.1
//#endif

/**
  Tools for computing the Census Transform of an image and comparing transformed pixels
*/

namespace vw {

//============================================================================
// Function declarations


/// Functions to compute the Census transform at one location in an image.
/// It is up to the user to perform bounds checking before using these functions!
// TODO: Move them somewhere, optimize the speed.
uint8  get_census_value_3x3(ImageView<uint8> const& image, int col, int row);
uint32 get_census_value_5x5(ImageView<uint8> const& image, int col, int row);
uint64 get_census_value_7x7(ImageView<uint8> const& image, int col, int row);


// TODO: Consolidate with code in Matcher.h!

/// Simple, unoptimized code for computing the hamming distance of two integers.
size_t hamming_distance(uint8  a, uint8  b);
size_t hamming_distance(uint32 a, uint32 b);
size_t hamming_distance(uint64 a, uint64 b);

// If we need a function to compute the costs of an entire image, this is a good
//  location for it!


//============================================================================
//============================================================================
//============================================================================
// Function definitions


uint8 get_census_value_3x3(ImageView<uint8> const& image, int col, int row) {
  // This will be an 8 bit sequence.
  uint8 output = 0;
  uint8 center = image(col, row);
  if (image(col-1, row-1) > center) output += 128;
  if (image(col  , row-1) > center) output +=  64;
  if (image(col+1, row-1) > center) output +=  32;
  if (image(col-1, row  ) > center) output +=  16;
  if (image(col+1, row  ) > center) output +=   8;
  if (image(col-1, row+1) > center) output +=   4;
  if (image(col  , row+1) > center) output +=   2;
  if (image(col+1, row+1) > center) output +=   1;
  return output;
}
uint32 get_census_value_5x5(ImageView<uint8> const& image, int col, int row) {
  // This will be a 24 bit sequence.
  uint32 output = 0;
  uint32 addend = 1;
  uint32 center = image(col, row);
  for (int r=row+2; r>=row-2; --r) {
    for (int c=col+2; c>=col-2; --c) {
      if ((r == row) && (c==col)) // Skip the central pixel
        continue;
      if (image(c,r) > center)
        output += addend;
      addend *=2;
    }
  }
  return output;
}

uint64 get_census_value_7x7(ImageView<uint8> const& image, int col, int row) {
  // This will be a 48 bit sequence.
  uint64 output = 0;
  uint64 addend = 1;
  uint64 center = image(col, row);
  for (int r=row+3; r>=row-3; --r) {
    for (int c=col+3; c>=col-3; --c) {
      if ((r == row) && (c==col)) // Skip the central pixel
        continue;
      if (image(c,r) > center)
        output += addend;
      addend *=2;
    }
  }
  return output;
}

// TODO: Consolidate with code in Matcher.h!


size_t hamming_distance(uint8 a, uint8 b) {
    uint8 dist = 0;
    uint8 val = a ^ b; // XOR

    // Count the number of bits set
    while (val != 0) {
        // A bit is set, so increment the count and clear the bit
        ++dist;
        val &= val - 1;
    }
    return static_cast<size_t>(dist); // Return the number of differing bits
}

size_t hamming_distance(uint32 a, uint32 b) {
    uint32 dist = 0;
    uint32 val = a ^ b; // XOR

    // Count the number of bits set
    while (val != 0) {
        // A bit is set, so increment the count and clear the bit
        ++dist;
        val &= val - 1;
    }
    return static_cast<size_t>(dist); // Return the number of differing bits
}

size_t hamming_distance(uint64 a, uint64 b) {
    uint64 dist = 0;
    uint64 val = a ^ b; // XOR

    // Count the number of bits set
    while (val != 0) {
        // A bit is set, so increment the count and clear the bit
        ++dist;
        val &= val - 1;
    }
    return static_cast<size_t>(dist); // Return the number of differing bits
}


} // end namespace vw

