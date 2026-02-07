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

/// \file ImageInterestData.h
///
#ifndef __IMAGE_INTEREST_DATA_H__
#define __IMAGE_INTEREST_DATA_H__

#include <vector>
#include <list>
#include <algorithm>
#include <sstream>

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Functors.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/Algorithms.h>
#include <vw/InterestPoint/InterestTraits.h>

#if defined(VW_HAVE_PKG_OPENCV) && (VW_HAVE_PKG_OPENCV == 1)
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#endif

namespace vw { namespace ip {

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
    typedef          SrcT             source_type;
    typedef typename SrcT::pixel_type pixel_type;

    typedef typename InterestOperatorTraits<SrcT, InterestT>::rasterize_type rasterize_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::gradient_type  gradient_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::mag_type       mag_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::ori_type       ori_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::interest_type  interest_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::integral_type  integral_type;

    static const int peak_type = InterestPeakType<InterestT>::peak_type;

    /// Constructor which sets the source image and creates the processed views.
    /// This is a generic constructor that assumes the requires grad x/y.
    template <class ViewT>
    ImageInterestData(ImageViewBase<ViewT> const& img) :
      m_src(img.impl()),
      m_grad_x(derivative_filter(m_src, 1, 0).impl()),
      m_grad_y(derivative_filter(m_src, 0, 1).impl()),
      m_mag(hypot(m_grad_x, m_grad_y).impl()),
      m_ori(atan2(m_grad_y, m_grad_x).impl()),
      m_interest(NULL),
      m_integral(NULL) {}

    template <class ViewT>
    ImageInterestData(ImageViewBase<ViewT> const& img,
                      ImageViewBase<integral_type> const& integral ) :
      m_src(img.impl()),
      m_interest(NULL),
      m_integral(&integral.impl()) {}

    ~ImageInterestData() {
      if (m_interest) delete m_interest;
    }

    /// Accessors to immutable processed views.
    inline source_type   const& source     () const { return m_src;    }
    inline gradient_type const& gradient_x () const { return m_grad_x; }
    inline gradient_type const& gradient_y () const { return m_grad_y; }
    inline mag_type      const& magnitude  () const { return m_mag;    }
    inline ori_type      const& orientation() const { return m_ori;    }

    /// Accessors to mutable interest image.
    inline interest_type& interest() const {
      if (!m_interest) 
        vw_throw(LogicErr() << "ImageInterestData::interest() Interest image has not yet been computed.");
      return *m_interest;
    }

    template <class ViewT>
    inline void set_interest(ImageViewBase<ViewT> const& interest) {
      if (m_interest) 
        delete m_interest;
      m_interest = new interest_type(interest.impl());
    }

    /// Accessors to mutable integral image.
    inline integral_type const& integral() const {
      if (!m_integral) 
        vw_throw(LogicErr() << "ImageInterestData::integral() Integral image has not yet been computed.");
      return *m_integral;
    }

    template <class ViewT>
    inline void set_integral(ImageViewBase<ViewT> const& integral) {
      // I don't really recommend using this -ZMM
      m_integral = &integral;
    }

  protected:
    /// Cached processed data
    source_type   m_src;
    gradient_type m_grad_x, m_grad_y;
    mag_type      m_mag;
    ori_type      m_ori;
          interest_type *m_interest;
    const integral_type *m_integral;
  }; // End class ImageInterestData

}} // namespace vw::ip

#endif //__IMAGE_INTEREST_DATA_H__
