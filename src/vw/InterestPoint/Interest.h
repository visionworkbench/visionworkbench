// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
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

/// \file Interest.h
/// 
/// Basic classes and functions for calculating interest images.
/// 
#ifndef _INTEREST_POINT_INTEREST_H_
#define _INTEREST_POINT_INTEREST_H_

#include <vector>
#include <stdio.h>

#include <vw/Image/Manipulation.h>
#include <vw/Image/Algorithms.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/FileIO.h>

namespace vw { namespace ip {

  /// Returns "cornerness" image, where the local maxima correspond to corners.
  /// By default uses Noble measure of corner strength (requires no tuning).
  /// Also supports Harris measure if positive k is specified (typical values:
  /// 0.04 <= k <= 0.15).
  template <class T>
  void harris_interest(ImageInterestData<T> const& data, T k = -1.0, T scale = 1.0) {
    typedef ImageView<T> Image;

    // Calculate elements of Harris matrix
    std::vector<T> kernel;
    generate_gaussian_kernel(kernel, scale, 0);
    Image Ix2 = separable_convolution_filter(data.grad_x * data.grad_x,
					     kernel, kernel);
    Image Iy2 = separable_convolution_filter(data.grad_y * data.grad_y,
					     kernel, kernel);
    Image Ixy = separable_convolution_filter(data.grad_x * data.grad_y,
					     kernel, kernel);

    // Estimate "cornerness"
    Image trace = Ix2 + Iy2;
    Image det = Ix2 * Iy2 - Ixy * Ixy;
    if (k < 0) {
      // Noble measure (preferred)
      data.interest = det / (trace + 0.000001);
    } else {
      // Standard Harris corner measure
      data.interest = det - k * trace * trace;
    }
  }

  /// Returns Harris "cornerness" image for a source image without the
  /// gradients precalculated.
  template <class T>
  ImageView<T> harris_interest(ImageView<T> const& image, T k = -1.0, T scale = 1.0) {
    ImageInterestData<T> data(image);
    harris_interest(data, k, scale);
    return data.interest;
  }

  /// Scale-normalized Laplacian of Gaussian interest measure.
  template <class T>
  void log_interest(ImageInterestData<T> const& data, T scale = 1.0) {
    data.interest = scale * laplacian_filter(data.src);
  }

  /// Scale-normalized Laplacian of Gaussian interest measure without
  /// the gradients precalculated.
  template <class T>
  ImageView<T> log_interest(ImageView<T> const& image, T scale = 1.0) {
    ImageInterestData<T> data(image);
    log_interest(data, scale);
    return data.interest;
  }

  /// Abstract base class for interest measures. An interest class should
  /// implement compute_interest and set its peak type in its constructor.
  template <class T>
  class InterestBase {
  protected:
    PeakType type;

  public:
    /// Store the interest image for data.src in data.interest, using the
    /// specified scale (if applicable).
    virtual int compute_interest(ImageInterestData<T> const& data, T scale = 1.0) {
      return 0;
    }

    PeakType peak_type() { return type; }
  };

  /// Class for computing Harris interest images.
  template <class T>
  class HarrisInterest : public InterestBase<T> {
  protected:
    T k;
    T v2; // Relative integration scale parameter (squared)

  public:
    HarrisInterest(T k_in = -1.0, T v2_in = 2.0) : k(k_in), v2(v2_in) { this->type = IP_MAX; }

    virtual int compute_interest(ImageInterestData<T> const& data, T scale = -1.0) {
      if (scale < 0)
	harris_interest(data, k);
      else
	harris_interest(data, k, scale / v2);
      return 0;
    }
  };

  /// Class for computing Laplacian of Gaussian interest images.
  template <class T>
  class LoGInterest : public InterestBase<T> {
  public:
    LoGInterest() { this->type = IP_MINMAX; }

    virtual int compute_interest(ImageInterestData<T> const& data, T scale = 1.0) {
      log_interest(data, scale);
      return 0;
    }
  };

} } //namespace vw::ip

#endif
