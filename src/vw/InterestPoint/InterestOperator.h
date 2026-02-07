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


/// \file InterestOperator.h
///
/// Basic classes and functions for calculating interest images.
///
#ifndef __VW_INTERESTPOINT_INTERESTVIEW_H__
#define __VW_INTERESTPOINT_INTERESTVIEW_H__

// STL Headers
#include <vector>

// Vision Workbench Headers
#include <vw/Image/ImageViewRef.h>
#include <vw/InterestPoint/InterestTraits.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/ImageInterestData.h>

namespace vw {
namespace ip {

  // These various InterestOperator classes (and the one on IntegralInterestOperator.h)
  //  all conform to a class pattern but do not derive from anything.

  // --------------------------- ------------------------ ----------------------------
  // --------------------------- Harris Interest Operator ----------------------------
  // --------------------------- ------------------------ ----------------------------

  /// Harris Corner interest operator
  class HarrisInterestOperator {
    double m_k, m_threshold;

  public:
    template <class ViewT> struct ViewType {
      typedef ImageViewRef<typename ViewT::pixel_type> type;
    };

    HarrisInterestOperator(double threshold = 1e-5, double k = -1.0) : m_k(k), m_threshold(threshold) {}

    /// Returns "cornerness" image, where the local maxima correspond to corners.
    /// By default uses Noble measure of corner strength (requires no tuning).
    /// Also supports Harris measure if positive k is specified (typical values:
    /// 0.04 <= k <= 0.15).
    template <class DataT>
    inline void operator() (DataT& data, float scale = 1.0) const {
      typedef typename DataT::source_type::pixel_type pixel_type;

      // Calculate elements of Harris matrix
      std::vector<float> kernel;
      generate_gaussian_kernel(kernel, scale, 0);
      ImageView<pixel_type> Ix2 = separable_convolution_filter(data.gradient_x() * data.gradient_x(),
                                                               kernel, kernel);
      ImageView<pixel_type> Iy2 = separable_convolution_filter(data.gradient_y() * data.gradient_y(),
                                                               kernel, kernel);
      ImageView<pixel_type> Ixy = separable_convolution_filter(data.gradient_x() * data.gradient_y(),
                                                               kernel, kernel);

      // Estimate "cornerness"
      ImageView<pixel_type> trace = Ix2 + Iy2;
      ImageView<pixel_type> det = Ix2 * Iy2 - Ixy * Ixy;
      if (m_k < 0) {
        // Noble measure (preferred)
        data.set_interest(det / (trace + 0.000001));
      } else {
        // Standard Harris corner measure
        data.set_interest(det - m_k * trace * trace);
      }
    }

    template <class ViewT>
    inline ImageViewRef<typename ViewT::pixel_type>
    operator() (ImageViewBase<ViewT> const& source, float scale = 1.0) const {
      ImageInterestData<ViewT, HarrisInterestOperator> data(source.impl());
      this->operator()(data, scale);
      return data.interest;
    }

    template <class DataT>
    inline bool threshold (InterestPoint const& pt, DataT const& /*data*/) const {
      return (pt.interest > m_threshold);
    }
  };

  /// Type traits for Harris interest
  template <> struct InterestPeakType <HarrisInterestOperator> { static const int peak_type = IP_MAX; };


//   /// Harris interest measure uses the image gradients.
//   template <class SrcT>
//   struct InterestOperatorTraits<SrcT, HarrisInterestOperator> {
//     typedef typename DefaultRasterizeT<SrcT>::type          rasterize_type;
//     typedef ImageView<typename SrcT::pixel_type>            gradient_type;
//     typedef typename DefaultMagT<SrcT>::type                mag_type;
//     typedef typename DefaultOriT<SrcT>::type                ori_type;
//     typedef rasterize_type                                  interest_type;
//   };

  // --------------------------- ------------------------ ----------------------------
  // --------------------- Laplacian of Gaussian Interest Operator -------------------
  // --------------------------- ------------------------ ----------------------------

  /// Log interest functor
  class LogInterestOperator {
    double m_threshold;

  public:
    template <class ViewT> struct ViewType {
      typedef ImageViewRef<typename ViewT::pixel_type> type;
    };

    LogInterestOperator(float threshold = 0.03) : m_threshold(threshold) {}

    template <class ViewT>
    inline ImageViewRef<typename ViewT::pixel_type>
    operator() (ImageViewBase<ViewT> const& source, float scale = 1.0) const {
      return scale * laplacian_filter(source.impl());
    }

    // TODO: this should return something
    template <class DataT>
    inline void operator() (DataT& data, float scale = 1.0) const {
      data.set_interest(ImageViewRef<typename DataT::source_type::pixel_type>( scale * laplacian_filter(data.source())));
    }

    template <class DataT>
    inline bool threshold (InterestPoint const& pt, DataT const& /*data*/) const {
      return (fabs(pt.interest) > m_threshold);
    }
  };

  /// Type traits for Log interest
  template <> struct InterestPeakType <LogInterestOperator> { static const int peak_type = IP_MINMAX; };

}} //namespace vw::ip

#endif // __VW_INTERESTPOINT_INTERESTVIEW_H__
