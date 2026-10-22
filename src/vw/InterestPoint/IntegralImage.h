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


/// IntegralImage.h
///
/// Provides support for algorithms that evaluate with Integral Images
/// such as SURF.  Contains a function to create an integral image, and
/// multiple functions to perform operations on them.
#ifndef __VW_INTERESTPOINT_INTEGRALIMAGE_H__
#define __VW_INTERESTPOINT_INTEGRALIMAGE_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Image/ImageView.h>

// TODO: Change the function names to meet the standard convention!

namespace vw {
namespace ip {

  /// Function to create an integral image of an input image.
  /// - Despite the caps, this is a function and IntegralImage is not a type!
  /// - An integral image can be used to quickly find regional sums using
  ///   the function below.
  /// - Returns ImageView<float>. All callers use float.
  template <class ViewT>
  inline ImageView<float>
  IntegralImage(ImageViewBase<ViewT> const& source) {

    // Allocating space
    ImageView<float> integral(source.impl().cols() + 1,
                              source.impl().rows() + 1);

    typedef typename ImageView<float>::pixel_accessor dst_accessor;
    typedef typename ViewT::pixel_accessor src_accessor;

    // Zero-ing the first row and col
    dst_accessor dest_row = integral.origin();
    *dest_row = 0;
    for (int i = 1; i < integral.cols(); i++) {
      dest_row.next_col();
      *dest_row = 0;
    }
    dest_row = integral.origin();
    for (int i = 1; i < integral.rows(); i++) {
      dest_row.next_row();
      *dest_row = 0;
    }

    // Pointer and pointer offset
    int offset = integral.cols();
    float* p = nullptr;

    // Now filling in integral image
    dest_row = integral.origin();
    dest_row.advance(1, 1);
    src_accessor src_row = source.impl().origin();
    for (int iy = 0; iy < source.impl().rows(); iy++) {
      dst_accessor dest_col = dest_row;
      src_accessor src_col = src_row;
      for (int ix = 0; ix < source.impl().cols(); ix++) {
        p = &(*dest_col);

        // Summing
        *dest_col = pixel_cast<PixelGray<float>>(*src_col).v()
                    + *(p - 1) + *(p - offset) - *(p - offset - 1);

        dest_col.next_col();
        src_col.next_col();
      }
      dest_row.next_row();
      src_row.next_row();
    }

    return integral;
  } // End IntegralImage function

  /// Using an integral image, compute the summed value of a region
  /// in the original image.
  float IntegralBlock(ImageView<float> const& integral,
                      Vector2i         const& top_left,
                      Vector2i         const& bottom_right);

  // Horizontal Wavelet
  // - integral  = Integral used for calculations
  // - x         = x location to evaluate at
  // - y         = y location to evaluate at
  // - size      = side of the square used for evaluate
  // For an interpolated integral view x and y may be fractional.
  template <class ViewT, class NumberT>
  inline float HHaarWavelet( ImageViewBase<ViewT> const& integral,
                             NumberT const& x, NumberT const& y,
                             float const& size ) {

    float response;
    float half_size = size / 2.0;
    float top  = float(y) - half_size;
    float left = float(x) - half_size;

    response = -integral.impl()(left, top);
    response += 2*integral.impl()(left+half_size, top);
    response -= integral.impl()(left+size, top);
    response += integral.impl()(left, top+size);
    response -= 2*integral.impl()(left+half_size, top+size);
    response += integral.impl()(left+size, top+size);

    return response;
  }

  // Vertical Wavelet
  // - integral  = Integral used for calculations
  // - x         = x location to evaluate at
  // - y         = y location to evaluate at
  // - size      = side of the square used for evaluate
  // For an interpolated integral view x and y may be fractional.
  template <class ViewT, class NumberT>
  inline float VHaarWavelet( ImageViewBase<ViewT> const& integral,
                             NumberT const& x, NumberT const& y,
                             float const& size ) {

    float response;
    float half_size = size / 2.0;
    float top  = float(y) - half_size;
    float left = float(x) - half_size;

    response = -integral.impl()(left, top);
    response += integral.impl()(left+size, top);
    response += 2*integral.impl()(left, top+half_size);
    response -= 2*integral.impl()(left+size, top+half_size);
    response -= integral.impl()(left, top+size);
    response += integral.impl()(left+size, top+size);

    return response;
  }

}} // end namespace vw::ip

#endif // __VW_INTERESTPOINT_INTEGRALIMAGE_H__
