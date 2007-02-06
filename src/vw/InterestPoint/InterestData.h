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

/// \file InterestData.h
/// 
/// Basic classes and structures for storing image interest points.
/// 
#ifndef __INTEREST_DATA_H__
#define __INTEREST_DATA_H__

#include <vw/Math/Vector.h>
#include <vw/Image/Filter.h>

namespace vw { 
namespace ip {

  /// A class for storing information about an interest point.
  struct InterestPoint {

    /// subpixel (col,row) location of point
    float x,y;

    /// scale of point.  This may come from the pyramid level, from
    /// interpolating the interest function between levels, or from some
    /// other scale detector like the Laplace scale used by Mikolajczyk
    /// & Schmid
    float scale;

    /// Integer location (unnormalized), mainly for internal use.
    int ix;
    int iy;

    /// Since the orientation is not necessarily unique, we may have more
    /// than one hypothesis for the orientation of an interest point.  I
    /// considered making a vector of orientations for a single point.
    /// However, it is probably better to make more than one interest
    /// point with the same (x,y,s) since the descriptor will be unique
    /// for a given orientation anyway.
    float orientation;

    float interest;

    /// And finally the descriptor for the interest point.  SIFT points
    /// have a vector of integers, PCA-SIFT features have a vector of
    /// floats or doubles...
    vw::Vector<float> descriptor;

    int size() const { return 2; }
    float operator[] (int index) const { 
      if (index == 0) return x;
      else if (index == 1) return y;
      else vw_throw( vw::ArgumentErr() << "Interest Point: Invalid index" );
    }

    bool operator< (const InterestPoint& other) const {
      return (other.interest < interest);
    }
  };

  template <class PixelT>
  struct ImageInterestData {
    ImageView<PixelT> src;
    ImageView<PixelT> grad_x;
    ImageView<PixelT> grad_y;
    ImageView<PixelT> ori;
    ImageView<PixelT> mag;
    ImageView<PixelT> interest;

    ImageInterestData() { }

    ImageInterestData(const ImageView<PixelT>& img) {
      set_source(img);
    }

    void set_source(const ImageView<PixelT>& img) {
      src = img;
      grad_x = derivative_filter(src, 1, 0);
      grad_y = derivative_filter(src, 0, 1);
      ori = atan2(grad_y, grad_x);
      mag = hypot(grad_x, grad_y);
      interest.set_size(img.cols(), img.rows());
    }
  };

  enum PeakType { IP_MAX, IP_MIN, IP_MINMAX };

}} // namespace vw::ip

#endif //__INTEREST_DATA_H__
