// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

/// \file MTBAlign.h
///
/// Routine for finding the translational offset between two images
/// using the mean threshold bitmap algorithm.
#ifndef __VW_HDR_MTBALIGN_H__
#define __VW_HDR_MTBALIGN_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Transform.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/ImageMath.h>
#include <boost/type_traits.hpp>
#include <algorithm>
#include <vector>
#include <iostream>

namespace vw {
namespace hdr {

class BooleanThresholdFunctor : public ReturnFixedType<int> {
private:
  double m_thresh;
public:
  BooleanThresholdFunctor(double thresh) : m_thresh(thresh) { }

  inline int operator()( PixelGray<float> const& pixel ) const {
    return (pixel[0] > m_thresh) ? 1 : 0;
  }
};

class ExclusionBitmapFunctor : public ReturnFixedType<int> {
private:
  double m_median, m_tolerance;
public:
  ExclusionBitmapFunctor(double median, double tolerance) : m_median(median), m_tolerance(tolerance) { }

  inline int operator()( PixelGray<float> const& pixel ) const {
    return (std::abs(pixel[0] - m_median) < m_tolerance) ? 0 : 1;
  }
};


template <class ChannelT>
void compute_bitmaps(ImageView<PixelGray<ChannelT> > img, ImageView<int> &tb, ImageView<int> &eb, ChannelT thresh = -1) {
  if (thresh < 0) {
    thresh = mean_channel_value(img);
    std::cout << "threshold: " << thresh << "\n";
  }
  tb = channel_cast<int>(threshold(select_channel(img,0), thresh));
  eb = threshold(abs(select_channel(img,0) - thresh),0.02, 1.0, 0.0); 
//  tb = per_pixel_filter(img, BooleanThresholdFunctor(threshold));
//  eb = per_pixel_filter(img, ExclusionBitmapFunctor(threshold, 0.02));
}

// xor of two bitmaps
struct PixelPixelXorFunctor : ReturnFixedType<int> {
  inline int operator()( int const& pixel1, int const& pixel2 ) const {
    return ((bool)pixel1 != (bool)pixel2) ? 1 : 0;
  }
};

BinaryPerPixelView<ImageView<int>, ImageView<int>, PixelPixelXorFunctor>
inline bitmap_xor( ImageView<int> const& image1, ImageView<int> const& image2 )
{
  typedef PixelPixelXorFunctor functor_type;
  typedef BinaryPerPixelView<ImageView<int>, ImageView<int>, functor_type> result_type;
  return result_type( image1.impl(), image2.impl(), functor_type() );
}

// and of two bitmaps
struct PixelPixelAndFunctor : ReturnFixedType<int> {
  inline int operator()( int const& pixel1, int const& pixel2 ) const {
    return (pixel1 && pixel2) ? 1 : 0;
  }
};

BinaryPerPixelView<ImageView<int>, ImageView<int>, PixelPixelAndFunctor>
inline bitmap_and( ImageView<int> const& image1, ImageView<int> const& image2 )
{
  typedef PixelPixelAndFunctor functor_type;
  typedef BinaryPerPixelView<ImageView<int>, ImageView<int>, functor_type> result_type;
  return result_type( image1.impl(), image2.impl(), functor_type() );
}

template <class ChannelT>
void get_exp_shift(ImageView<PixelGray<ChannelT> > &img1, ImageView<PixelGray<ChannelT> > &img2, int shift_bits, int shift[2], ChannelT threshold1 = -1, ChannelT threshold2 = -1) {
  typedef ImageView<PixelGray<ChannelT> > Image;
  typedef ImageView<int> Bitmap; // mock-up for now

  int min_err;
  int cur_shift[2];
  Bitmap tb1, tb2;
  Bitmap eb1, eb2;

  std::cout << "Entering get_exp_shift level " << shift_bits << "\n";
  if (shift_bits > 0) {
    std::cout << "\tSubsampling... ";
    Image sml_img1 = subsample(gaussian_filter(img1, 1), 2);
    Image sml_img2 = subsample(gaussian_filter(img2, 1), 2);
    std::cout << "done\n";
    get_exp_shift(sml_img1, sml_img2, shift_bits - 1, cur_shift);
    cur_shift[0] *= 2;
    cur_shift[1] *= 2;
  } else {
    cur_shift[0] = cur_shift[1] = 0;
  }

  std::cout << "Computing bitmaps level " << shift_bits << "\n";
  compute_bitmaps(img1, tb1, eb1, threshold1);
  compute_bitmaps(img2, tb2, eb2, threshold2);
  min_err = img1.rows() * img1.cols();

//   write_image("tb1.jpg", tb1);
//   write_image("tb2.jpg", tb2);
//   write_image("eb1.jpg", eb1);
//   write_image("eb2.jpg", eb2);

  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      int xs = cur_shift[0] + i;
      int ys = cur_shift[1] + j;
      Bitmap shifted_tb2 = translate(tb2, xs, ys);
      Bitmap shifted_eb2 = translate(eb2, xs, ys);
      Bitmap diff_b = bitmap_xor(tb1, shifted_tb2);
      diff_b = bitmap_and(diff_b, eb1);
      diff_b = bitmap_and(diff_b, shifted_eb2);
      long int err = sum_of_channel_values(diff_b);
      std::cout << "Error: " << err << "  Min error: " << min_err << "\n";
      if (err < min_err) {
        shift[0] = xs;
        shift[1] = ys;
        min_err = err;
      }
    }
  }
  std::cout << "\tshift = " << shift[0] << " " << shift[1] << "\n";
}

template <class PixelT>
void mtb_align(std::vector<ImageView<PixelRGB<PixelT> > > &images) {

  for (int i = 0; i < images.size() - 1; i++) {
    ImageView<PixelGray<typename PixelChannelType<PixelT>::type > > base_image = images[i];
    ImageView<PixelGray<typename PixelChannelType<PixelT>::type > > image_to_align = images[i+1];

    int shift_bits = (int)roundf(log(images[0].cols()) / log(2)) - 4;
    if (shift_bits < 0) shift_bits = 0;
    std::cout << "Shift bits: " << shift_bits << "\n";

    int shift[2];
    get_exp_shift(base_image, image_to_align, shift_bits, shift);
    std::cout << "Shift from " << i << " to " << i+1 << ": " << shift[0] << " " << shift[1] << "\n";
    images[i+1] = translate(copy(images[i+1]), shift[0], shift[1]);
    
    std::ostringstream os;
    os << "image" << i+1 << ".jpg";
    write_image(os.str(), images[i+1]);

  }
}

}} // namespace vw::hdr

#endif  // __VW_HDR_MTBAlign_H__
