// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

/// \file IntegralDetector.h
///
/// Detectors that operate off of the use of integral images
///

#ifndef __VW_INTERESTPOINT_INTEGRAL_DETECTOR_H__
#define __VW_INTERESTPOINT_INTEGRAL_DETECTOR_H__

#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/IntegralImage.h>
#include <vw/InterestPoint/IntegralInterestOperator.h>
#include <deque>

namespace vw {
namespace ip {

/// InterestDetector implementation for all detectors with operate off of an integral image.
class IntegralInterestPointDetector:
  public InterestDetectorBase<IntegralInterestPointDetector>,
  private boost::noncopyable {

public:
  static const int IP_DEFAULT_SCALES = 8;
  typedef vw::ImageView<float> ImageT;
  typedef ImageInterestData<ImageT,OBALoGInterestOperator> DataT;

  /// Setting max_points = 0 will disable interest point culling.
  IntegralInterestPointDetector(int max_points = 1000, int num_scales = IP_DEFAULT_SCALES):
    m_interest(OBALoGInterestOperator()), m_scales(num_scales), m_max_points(max_points) {}

  IntegralInterestPointDetector(OBALoGInterestOperator const& interest, int max_points = 1000):
    m_interest(interest), m_scales(IP_DEFAULT_SCALES), m_max_points(max_points) {}

  IntegralInterestPointDetector(OBALoGInterestOperator const& interest, int scales, 
                                int max_points):
    m_interest(interest), m_scales(scales), m_max_points(max_points) {}

  /// Detect interest points in the source image
  InterestPointList process_image(vw::ImageViewRef<float> const& image,
                                  int desired_num_ip = 0) const;

protected:

  // Helper function
  void threshold(InterestPointList      & points,
                 DataT             const& img_data,
                 int               const& scale) const;

  OBALoGInterestOperator m_interest;
  int m_scales, m_max_points;

};

/// Implementation of IntegralInterestPointDetector based on OBALoGInterestOperator 
class IntegralAutoGainDetector: public InterestDetectorBase<IntegralAutoGainDetector>,
                                private boost::noncopyable {
public:

  static const int IP_DEFAULT_SCALES = 8;
  typedef vw::ImageView<float> ImageT;
  typedef vw::ip::ImageInterestData<ImageT, vw::ip::OBALoGInterestOperator> DataT;

  IntegralAutoGainDetector(size_t max_points = 200, size_t scales = IP_DEFAULT_SCALES):
    IntegralAutoGainDetector(OBALoGInterestOperator(0), scales, max_points) {}

  IntegralAutoGainDetector(OBALoGInterestOperator const& interest, int scales,
                            int max_points):
  m_interest(interest), m_scales(scales), m_max_points(max_points) {}

  // Detect interest points
  InterestPointList process_image(vw::ImageViewRef<float> const& image,
                                  int desired_num_ip = 0) const;

protected:

  OBALoGInterestOperator m_interest;
  int m_scales, m_max_points;

  // Helper function
  void threshold(vw::ip::InterestPointList& points,
                 DataT const& img_data,
                 int const& scale, float threshold_lvl) const;
  
}; // End class IntegralAutoGainDetector

}} // end vw::ip

#endif//__VW_INTERESTPOINT_INTEGRAL_DETECTOR_H__
