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

#include <vw/Math/Functions.h>
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


// TODO: Move to a cpp file

/// Functions to compute the Census transform at one location in an image.
/// It is up to the user to perform bounds checking before using these functions!
/// - The ternary transform is from the paper "TEXTURE-AWARE DENSE IMAGE MATCHING USING TERNARY CENSUS TRANSFORM"
///   by Han Hu,Chongtai Chen, Bo Wu, Xiaoxia Yang, Qing Zhu, Yulin Ding
// TODO: Move them somewhere, optimize the speed.
inline uint8  get_census_value_3x3        (ImageView<uint8> const& image, int col, int row);
inline uint32 get_census_value_5x5        (ImageView<uint8> const& image, int col, int row);
inline uint64 get_census_value_7x7        (ImageView<uint8> const& image, int col, int row);
inline uint32 get_census_value_9x9        (ImageView<uint8> const& image, int col, int row);
inline uint16 get_census_value_ternary_3x3(ImageView<uint8> const& image, int col, int row, int diff_threshold=2);
inline uint64 get_census_value_ternary_5x5(ImageView<uint8> const& image, int col, int row, int diff_threshold=2);
inline uint64 get_census_value_ternary_7x7(ImageView<uint8> const& image, int col, int row, int diff_threshold=2);
inline uint64 get_census_value_ternary_9x9(ImageView<uint8> const& image, int col, int row, int diff_threshold=2);


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
  uint8 center = image(col, row);
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
  uint8 center = image(col, row);
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

uint32 get_census_value_9x9(ImageView<uint8> const& image, int col, int row) {
  // This will be a 32 bit sequence.
  // - Could expand this out to a 64 bit sequence
  uint32 output = 0;
  uint32 addend = 1;
  
  /* Pixels are checked in this pattern:
  *...*...*
  .*.*.*.*.
  ..*.*.*..
  .*..*..*.
  *.**.**.*
  .*..*..*.
  ..*.*.*..
  .*.*.*.*.
  *...*...*
  Taken from the paper "TEXTURE-AWARE DENSE IMAGE MATCHING USING TERNARY CENSUS TRANSFORM"
  by Han Hu,Chongtai Chen, Bo Wu, Xiaoxia Yang, Qing Zhu, Yulin Ding
  */
  
  const int NUM_POSITIONS = 32; // * 1 bits per position = 32 bits
  int cols[NUM_POSITIONS];
  int rows[NUM_POSITIONS];
  cols[ 0] = 0;  cols[ 1] = 4;  cols[ 2] = 8;
  cols[ 3] = 1;  cols[ 4] = 3;  cols[ 5] = 5;  cols[ 6] = 7;
  cols[ 7] = 2;  cols[ 8] = 4;  cols[ 9] = 6;
  cols[10] = 1;  cols[11] = 4;  cols[12] = 7;
  cols[13] = 0;  cols[14] = 2;  cols[15] = 3; cols[16] = 5;  cols[17] = 6;  cols[18] = 8;
  cols[19] = 1;  cols[20] = 4;  cols[21] = 7;
  cols[22] = 2;  cols[23] = 4;  cols[24] = 6;
  cols[25] = 1;  cols[26] = 3;  cols[27] = 5;  cols[28] = 7;
  cols[29] = 0;  cols[30] = 4;  cols[31] = 8;

  rows[ 0] = rows[ 1] = rows[ 2] = 0;
  rows[ 3] = rows[ 4] = rows[ 5] = rows[ 6] = 1;
  rows[ 7] = rows[ 8] = rows[ 9] = 2;
  rows[10] = rows[11] = rows[12] = 3;
  rows[13] = rows[14] = rows[15] = rows[16] = rows[17] = rows[18] = 4;
  rows[19] = rows[20] = rows[21] = 5;
  rows[22] = rows[23] = rows[24] = 6;
  rows[25] = rows[26] = rows[27] = rows[28] = 7;
  rows[29] = rows[30] = rows[31] = 8;
  
  const int center = image(col, row);
  
  for (int i=0; i<NUM_POSITIONS; ++i) {
    int val = image(col+cols[i]-4,row+rows[i]-4);
    if (val > center)
      output += addend;
    addend *=2;
  }
  return output;
}


uint16 get_census_value_ternary_3x3(ImageView<uint8> const& image, int col, int row,
                                    int diff_threshold) {
  // This will be a 16 bit sequence.
  uint16 output = 0;
  uint16 addend = 1;
  const int center = image(col, row);
  const int low_thresh  = center - diff_threshold;
  const int high_thresh = center + diff_threshold;
  for (int r=row+1; r>=row-1; --r) {
    for (int c=col+1; c>=col-1; --c) {
      if ((r == row) && (c==col)) // Skip the central pixel
        continue;
      int val = image(c,r);
      if (val >= low_thresh) {
        output += addend;
        if (val > high_thresh) // Greater, += 11
          output += addend*2;
        // else Middle range, += 01
      }
      // For low range, += 00

      addend *=4;
    }
  }
  return output;
}

uint64 get_census_value_ternary_5x5(ImageView<uint8> const& image, int col, int row,
                                    int diff_threshold) {
  // This will be a 48 bit sequence.
  uint64 output = 0;
  uint64 addend = 1;
  const int center = image(col, row);
  const int low_thresh  = center - diff_threshold;
  const int high_thresh = center + diff_threshold;
  for (int r=row+2; r>=row-2; --r) {
    for (int c=col+2; c>=col-2; --c) {
      if ((r == row) && (c==col)) // Skip the central pixel
        continue;
      int val = image(c,r);
      if (val >= low_thresh) {
        output += addend;
        if (val > high_thresh) // Greater, += 11
          output += addend*2;
        // else Middle range, += 01
      }
      // For low range, += 00

      addend *=4;
    }
  }
  return output;
}


uint64 get_census_value_ternary_7x7(ImageView<uint8> const& image, int col, int row,
                                    int diff_threshold) {
  // This will be a 64 bit sequence.
  uint64 output = 0;
  uint64 addend = 1;
  
  /* Pixels are checked in this pattern:
  *.***.*  
  .*.*.*.
  *.***.*
  ***.***
  *.***.*
  .*.*.*.
  *.***.*
  Custom pattern.
  */
  
  const int NUM_POSITIONS = 32; // * 2 bits per position = 64 bits
  int cols[NUM_POSITIONS];
  int rows[NUM_POSITIONS];
  cols[ 0] = 0;  cols[ 1] = 2;  cols[ 2] = 3;  cols[ 3] = 4;  cols[ 4] = 6;
  cols[ 5] = 1;  cols[ 6] = 3;  cols[ 7] = 5;
  cols[ 8] = 0;  cols[ 9] = 2;  cols[10] = 3;  cols[11] = 4;  cols[12] = 6;
  cols[13] = 0;  cols[14] = 1;  cols[15] = 2;  cols[16] = 4;  cols[17] = 5;  cols[18] = 6;
  cols[19] = 0;  cols[20] = 2;  cols[21] = 3;  cols[22] = 4;  cols[23] = 6;
  cols[24] = 1;  cols[25] = 3;  cols[26] = 5;
  cols[27] = 0;  cols[28] = 2;  cols[29] = 3;  cols[30] = 4;  cols[31] = 6;

  rows[ 0] = rows[ 1] = rows[ 2] = rows[ 3] = rows[ 4] = 0;
  rows[ 5] = rows[ 6] = rows[ 7] = 1;
  rows[ 8] = rows[ 9] = rows[10] = rows[11] = rows[12] = 2;
  rows[13] = rows[14] = rows[15] = rows[16] = rows[17] = rows[18] = 3;
  rows[19] = rows[20] = rows[21] = rows[22] = rows[23] = 4;
  rows[24] = rows[25] = rows[26] = 5;
  rows[27] = rows[28] = rows[29] = rows[30] = rows[31] = 6;
  
  const int center = image(col, row);
  const int low_thresh  = center - diff_threshold;
  const int high_thresh = center + diff_threshold;
  
  for (int i=0; i<NUM_POSITIONS; ++i) {
    int val = image(col+cols[i]-3,row+rows[i]-3);
    if (val >= low_thresh) {
      output += addend;
      if (val > high_thresh) // Greater, += 11
        output += addend*2;
      // else Middle range, += 01
    }
    // For low range, += 00

    addend *=4;
  }
  return output;
}

uint64 get_census_value_ternary_9x9(ImageView<uint8> const& image, int col, int row,
                                    int diff_threshold) {
  // This will be a 64 bit sequence.
  uint64 output = 0;
  uint64 addend = 1;
  
  /* Pixels are checked in this pattern:
  *...*...*
  .*.*.*.*.
  ..*.*.*..
  .*..*..*.
  *.**.**.*
  .*..*..*.
  ..*.*.*..
  .*.*.*.*.
  *...*...*
  Taken from the paper "TEXTURE-AWARE DENSE IMAGE MATCHING USING TERNARY CENSUS TRANSFORM"
  by Han Hu,Chongtai Chen, Bo Wu, Xiaoxia Yang, Qing Zhu, Yulin Ding
  */
  
  const int NUM_POSITIONS = 32; // * 2 bits per position = 64 bits
  int cols[NUM_POSITIONS];
  int rows[NUM_POSITIONS];
  cols[ 0] = 0;  cols[ 1] = 4;  cols[ 2] = 8;
  cols[ 3] = 1;  cols[ 4] = 3;  cols[ 5] = 5;  cols[ 6] = 7;
  cols[ 7] = 2;  cols[ 8] = 4;  cols[ 9] = 6;
  cols[10] = 1;  cols[11] = 4;  cols[12] = 7;
  cols[13] = 0;  cols[14] = 2;  cols[15] = 3; cols[16] = 5;  cols[17] = 6;  cols[18] = 8;
  cols[19] = 1;  cols[20] = 4;  cols[21] = 7;
  cols[22] = 2;  cols[23] = 4;  cols[24] = 6;
  cols[25] = 1;  cols[26] = 3;  cols[27] = 5;  cols[28] = 7;
  cols[29] = 0;  cols[30] = 4;  cols[31] = 8;

  rows[ 0] = rows[ 1] = rows[ 2] = 0;
  rows[ 3] = rows[ 4] = rows[ 5] = rows[ 6] = 1;
  rows[ 7] = rows[ 8] = rows[ 9] = 2;
  rows[10] = rows[11] = rows[12] = 3;
  rows[13] = rows[14] = rows[15] = rows[16] = rows[17] = rows[18] = 4;
  rows[19] = rows[20] = rows[21] = 5;
  rows[22] = rows[23] = rows[24] = 6;
  rows[25] = rows[26] = rows[27] = rows[28] = 7;
  rows[29] = rows[30] = rows[31] = 8;
  
  const int center = image(col, row);
  const int low_thresh  = center - diff_threshold;
  const int high_thresh = center + diff_threshold;
  
  for (int i=0; i<NUM_POSITIONS; ++i) {
    int val = image(col+cols[i]-4,row+rows[i]-4);
    if (val >= low_thresh) {
      output += addend;
      if (val > high_thresh) // Greater, += 11
        output += addend*2;
      // else Middle range, += 01
    }
    // For low range, += 00

    addend *=4;
  }
  return output;
}

} // end namespace vw

