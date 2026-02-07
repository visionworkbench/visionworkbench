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

/// \file OpenCvDetector.h
///
/// Built-in classes and functions for interest point detection with OpenCV.
///

#ifndef __VW_INTEREST_POINT_OPENCV_DETECTOR_H__
#define __VW_INTEREST_POINT_OPENCV_DETECTOR_H__

#include <vw/vw_config.h> // defines VW_HAVE_PKG_OPENCV

#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1

#include <vw/Image/ImageViewRef.h>
#include <vw/InterestPoint/Detector.h>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

namespace vw {
namespace ip {

enum OpenCvIpDetectorType {OPENCV_IP_DETECTOR_TYPE_BRISK = 0,
                            OPENCV_IP_DETECTOR_TYPE_ORB   = 1,
                            OPENCV_IP_DETECTOR_TYPE_SIFT  = 2,
                            OPENCV_IP_DETECTOR_TYPE_SURF  = 3};

/// Struct to convert a basic type to a single channel OpenCV type
template <typename T>
struct GetOpenCvPixelType { static const int type=CV_8UC1; };
template <>
struct GetOpenCvPixelType<short> { static const int type=CV_16SC1; };
template <>
struct GetOpenCvPixelType<unsigned short> { static const int type=CV_16UC1; };
template <>
struct GetOpenCvPixelType<int> { static const int type=CV_32SC1; };
template <>
struct GetOpenCvPixelType<float> { static const int type=CV_32FC1; };
template <>
struct GetOpenCvPixelType<double> { static const int type=CV_64FC1; };

/// Get an OpenCV wrapper, rasterizing the VW image to a provided buffer.
void get_opencv_wrapper(vw::ImageViewRef<float> const& input_image,
                        cv::Mat & cv_image,
                        ImageView<vw::uint8> &image_buffer,
                        cv::Mat & cv_mask,
                        bool normalize = true);

/// Interest point detector build using OpenCV functions
class OpenCvInterestPointDetector:
  public InterestDetectorBase<OpenCvInterestPointDetector> {
public:

  OpenCvInterestPointDetector(OpenCvIpDetectorType detector_type,
                              bool normalize,
                              bool add_description, int max_points);

  /// Detect interest points in the source image. 
  InterestPointList process_image(vw::ImageViewRef<float> const& image,
                                  int desired_num_ip=0) const;

private:
  OpenCvIpDetectorType m_detector_type;
  bool                 m_add_descriptions;
  bool                 m_normalize;
  int                  m_max_points;
  cv::Ptr<cv::FeatureDetector> m_detector;

  /// Initialize an OpenCV detector object
  cv::Ptr<cv::FeatureDetector> init_detector(OpenCvIpDetectorType detector_type,
                                             int max_points) const;

}; // End class OpenCvInterestPointDetector

}} // namespace vw::ip

#endif // End case with OpenCV installed
#endif // __VW_INTEREST_POINT_OPENCV_DETECTOR_H__
