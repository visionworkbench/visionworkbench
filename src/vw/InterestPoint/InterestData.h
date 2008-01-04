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

#include <vw/InterestPoint/InterestTraits.h>
#include <vw/Math/Vector.h>
#include <vw/Image.h>

#include <vector>
#include <list>
#include <algorithm>

namespace vw { 
namespace ip {

  /// A class for storing information about an interest point.
  struct InterestPoint {

    InterestPoint() {}

    InterestPoint(int x, int y, float scale=1.0, float interest=0.0, float ori=0.0)
      : x((float)x), y((float)y), scale(scale), ix(x), iy(y), orientation(ori), interest(interest) {}

    /// Subpixel (col,row) location of point
    float x,y;

    /// Scale of point.  This may come from the pyramid level, from
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

    /// The interest measure (could be Harris, LoG, etc.).
    float interest;

    /// And finally the descriptor for the interest point.  For example, 
    /// PCA descriptors would have a vector of floats or doubles...
    vw::Vector<float> descriptor;
    
    vw::Vector<float>::iterator begin(){
      return descriptor.begin();
    }
    vw::Vector<float>::iterator end(){
      return descriptor.end();
    }
    
    int size() const { return 2; }

    float operator[] (int index) const { 
      if (index == 0) return x;
      else if (index == 1) return y;
      else vw_throw( vw::ArgumentErr() << "Interest Point: Invalid index" );

      // Control should never reach this point
      return 0;
    }

    /// std::sort can be used to sort InterestPoints in descending
    /// order of interest.
    bool operator< (const InterestPoint& other) const {
      return (other.interest < interest);
    }
  };

  // Need to use list instead of vector for efficient thresholding.
  typedef std::list<InterestPoint> InterestPointList;

  /// Select only the interest points that fall within the specified bounding box.
  template <class RealT>
  InterestPointList crop(InterestPointList const& interest_points, BBox<RealT,2> const& bbox) {
    InterestPointList return_val;
    for (InterestPointList::iterator i = interest_points.begin(); i != interest_points.end(); ++i) {
      if (bbox.contains(Vector<RealT,2>(RealT((*i).x), RealT((*i).y))))
        return_val.push_back(*i);
    }
    return return_val;
  }

  /// ImageInterestData
  ///
  /// This struct encapsulates some basic and widely useful processed
  /// views of a source image: the horizontal and vertical gradients,
  /// the orientation image, the gradient magnitude image, and the
  /// interest image. This is useful to ensure that these images are
  /// not redundantly calculated by different steps of the feature
  /// detection algorithm.
  ///
  /// The interest type is used to determine at compile-time which processed
  /// views should be fully rasterized. For speed in feature detection, the
  /// source type should be ImageView<T> or a simple manipulation of it. 
  /// For memory efficiency, the source type should be ImageViewRef<T>.
  ///
  /// If some other sort of shared data is needed or any of the temporaries
  /// should be calculated in a different fashion, ImageInterestData can be
  /// partially specialized on InterestT.
  template <class SrcT, class InterestT>
  class ImageInterestData {
  public:
    /// The image types defined by InterestTraits control whether each processed
    /// view is fully rasterized or not. Only those used in calculating each
    /// pixel's interest measure should be fully rasterized. Later operations
    /// (thresholding, orientation assignment, etc.) require at most support
    /// regions around the interest points.
    typedef SrcT source_type;
    typedef typename SrcT::pixel_type pixel_type;

    typedef typename InterestOperatorTraits<SrcT, InterestT>::rasterize_type rasterize_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::gradient_type gradient_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::mag_type mag_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::ori_type ori_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::interest_type interest_type;

    static const int peak_type = InterestPeakType<InterestT>::peak_type;

    /// Constructor which sets the source image and creates the processed views.
    template <class ViewT>
    ImageInterestData(ImageViewBase<ViewT> const& img) :
      m_src(img.impl()),
      m_grad_x(derivative_filter(m_src, 1, 0).impl()),
      m_grad_y(derivative_filter(m_src, 0, 1).impl()),
      m_mag(hypot(m_grad_x, m_grad_y).impl()),
      m_ori(atan2(m_grad_y, m_grad_x).impl()),
      m_interest(NULL) {}

    ~ImageInterestData() { if (m_interest) delete m_interest; }

    /// Accessors to immutable processed views.
    inline source_type const& source() const { return m_src; }
    inline gradient_type const& gradient_x() const { return m_grad_x; }
    inline gradient_type const& gradient_y() const { return m_grad_y; }
    inline mag_type const& magnitude() const { return m_mag; }
    inline ori_type const& orientation() const { return m_ori; }

    /// Accessors to mutable interest image.
    inline interest_type& interest() const { 
      if (!m_interest) vw_throw(LogicErr() << "ImageInterestData::interest() Interest image has not yet been computed.");
      return *m_interest; 
    }

    template <class ViewT>
    inline void set_interest(ImageViewBase<ViewT> const& interest) {
      if (m_interest) delete m_interest;
      m_interest = new interest_type(interest.impl());
    }

  protected:
    /// Cached processed data
    source_type m_src;
    gradient_type m_grad_x, m_grad_y;
    mag_type m_mag;
    ori_type m_ori;
    interest_type *m_interest;
  };


}} // namespace vw::ip

#endif //__INTEREST_DATA_H__
