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


/// \file IntegralInterestOperator.h
///
/// Classes that detect interest by using integral images
///
#ifndef __VW_INTEGRAL_INTEREST_OPERATOR_H__
#define __VW_INTEGRAL_INTEREST_OPERATOR_H__

// STL
#include <vector>

#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Filter.h>

// Interest Point Headers
#include <vw/InterestPoint/BoxFilter.h>
#include <vw/InterestPoint/InterestTraits.h>
#include <vw/InterestPoint/InterestPoint.h>

namespace vw {
namespace ip {


// TODO: Should this move to InterestOperator.h?

  /// Functions adapted from the source code of the OBALoG
  /// Detector. OBALoG is release under the 'new' BSD License below.
  /*
   * Optimized Box Approximation of Laplacian of Gaussian (OBALoG)
   * Authors : Zachary Moratto, Vinayak Jakkula, Chris Lewis, Dale Schinstock
   *
   * Copyright (c) (2009): <Kansas State University,Manhattan,Kansas>
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *     * Redistributions of source code must retain the above copyright
   *       notice, this list of conditions and the following disclaimer.
   *     * Redistributions in binary form must reproduce the above copyright
   *       notice, this list of conditions and the following disclaimer in the
   *       documentation and/or other materials provided with the distribution.
   *     * Neither the name of Kansas State University nor the
   *       names of its contributors may be used to endorse or promote products
   *       derived from this software without specific prior written permission.
   */


  // OBALoG Interest Operator
  // _____________________________________________________________
  class OBALoGInterestOperator {

    // Predefined functions
    static const double SCALE_LOG_SIGMA [10];     // For reference
    static const int    SCALE_BOX_WIDTH [10][6];
    static const int    SCALE_BOX_HEIGHT[10][6];
    static const double SCALE_BOX_WEIGHT[10][6];

    double m_threshold;
    std::vector<double> m_gaussian_kernel;

  public:
    template <class ViewT> struct ViewType {
      typedef ImageViewRef<typename ViewT::pixel_type> type;
    };

    OBALoGInterestOperator(double threshold = 0.05) : m_threshold(threshold) {
      generate_gaussian_kernel(m_gaussian_kernel,2,9);
    }

    // A floating point scale provided will invoke OBALoG trying to
    // tailor a box filter for that scale.
    template <class DataT>
    inline void operator() (DataT& /*data*/, float /*scale*/) const {
      vw_throw( NoImplErr() << "OBALoG Box Filter creation on the fly has not been implemented\n" );
    }

    // A integer scale will invoke the standard OBALoG with precomputed box filters.
    template <class DataT>
    inline void operator() (DataT& data, int scale = 0 ) const {
      // 1.) Assemble Filter
      BoxFilter bfilter;
      for ( uint8 b = 0; b < 6; b++ ) {
        SumBox instance;
        instance.size = Vector2i(SCALE_BOX_WIDTH[scale][b],
                                 SCALE_BOX_HEIGHT[scale][b]);
        instance.start = instance.size;
        instance.start[0] = instance.start[0] >> 1; // Divide by 2, round down
        instance.start[1] = instance.start[1] >> 1;
        instance.start[0] = -instance.start[0];
        instance.start[1] = -instance.start[1];
        instance.weight = SCALE_BOX_WEIGHT[scale][b];
        bfilter.push_back(instance);
      }

      // 2.) Apply Filter ( DataT performs actual render )
      data.set_interest(ImageViewRef<typename DataT::integral_type::pixel_type>( abs(box_filter(data.integral(), bfilter))));
    }

    // Threshold will reassign the interest with the harris corner detector
    template <class DataT>
    inline bool threshold (InterestPoint const& ip,
                           DataT const& data, int scale) const {
      return threshold( ip, data, float(SCALE_LOG_SIGMA[scale]) );
    }

    template <class DataT>
    inline bool threshold (InterestPoint const& ip,
                           DataT const& data, float scale) const {
      if ( ip.interest <  m_threshold )
        return false;

      // Threshold secretly also applies the harris operation;
      int step = int(scale);
      int offset = 4*step;
      float sum_Lx_2 = 0;
      float sum_Ly_2 = 0;
      float sum_LxLy = 0;
      for ( int y = ip.iy - offset, iy = 0;
            y <= ip.iy + offset; y += step, iy++ ) {
        for ( int x = ip.ix - offset, ix = 0;
              x <= ip.ix + offset; x += step, ix++ ) {
          float Lx = data.interest()(x-1,y) - data.interest()(x+1,y);
          float Ly = data.interest()(x,y-1) - data.interest()(x,y+1);

          float weight = m_gaussian_kernel[ix]*m_gaussian_kernel[iy];
          sum_Lx_2 += Lx*Lx*weight;
          sum_Ly_2 += Ly*Ly*weight;
          sum_LxLy += Lx*Ly*weight;
        }
      }

      float trace = sum_Lx_2 + sum_Ly_2;
      float det = sum_Lx_2*sum_Ly_2 - sum_LxLy*sum_LxLy;

      return (trace*trace)/det < 12;
    }

    inline float float_scale( int const& scale ) const {
      return SCALE_LOG_SIGMA[scale];
    }

  };

  // Type traits for OBALoG Interest
  template <> struct InterestPeakType <OBALoGInterestOperator> { static const int peak_type = IP_MAX; };



}} // end vw::ip

#endif//__VW_INTEGRAL_INTEREST_OPERATOR_H__
