// __BEGIN_LICENSE__
//  Copyright (c) 2006-2025, United States Government as represented by the
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
#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1

#include <vw/Core/Settings.h>
#include <vw/Core/ThreadPool.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Filter.h>
#include <vw/InterestPoint/Detector.h>

#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/Extrema.h>
#include <vw/InterestPoint/Localize.h>
#include <vw/InterestPoint/InterestOperator.h>
#include <vw/InterestPoint/WeightedHistogram.h>
#include <vw/InterestPoint/ImageOctaveHistory.h>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/features2d.hpp>
#include <vw/Image/Manipulation.h>

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
template <class ViewT>
void get_opencv_wrapper(ImageViewBase<ViewT> const& input_image,
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
  template <class ViewT>
  InterestPointList process_image(ImageViewBase<ViewT> const& image, int desired_num_ip=0) const;

private:
  OpenCvIpDetectorType m_detector_type;
  bool                 m_add_descriptions;
  bool                 m_normalize;
  int                  m_max_points;
  cv::Ptr<cv::FeatureDetector> m_detector;

  /// Initialize an OpenCV detector object
  cv::Ptr<cv::FeatureDetector> init_detector(OpenCvIpDetectorType detector_type, int max_points) const;

}; // End class OpenCvInterestPointDetector

template <class ViewT>
void get_opencv_wrapper(ImageViewBase<ViewT> const& input_image,
                        cv::Mat & cv_image,
                        ImageView<vw::uint8> &image_buffer,
                        cv::Mat & cv_mask,
                        bool normalize) {

  // Rasterize the input image so we don't suffer from slow disk access or something.
  ImageView<typename ViewT::pixel_type> input_buffer = input_image.impl();

  if (normalize) // Convert the input image to uint8 with 2%-98% intensity scaling.
    percentile_scale_convert(input_buffer, image_buffer, 0.02, 0.98);
  else {
    // Convert to uint8 using the default ranges for the input data type
    double standard_min = ChannelRange<typename ViewT::pixel_type>::min();
    double standard_max = ChannelRange<typename ViewT::pixel_type>::max();
    image_buffer = pixel_cast_rescale<vw::uint8>(clamp(input_buffer, standard_min, standard_max));
  }

  // Figure out the image buffer parameters
  int     cv_data_type = GetOpenCvPixelType<vw::uint8>::type;
  void*   raw_data_ptr = reinterpret_cast<void*>(image_buffer.data());
  size_t  pixel_size   = sizeof(vw::uint8);
  size_t  step_size    = image_buffer.cols() * pixel_size;

  // Create an OpenCV wrapper for the buffer image
  cv_image = cv::Mat(image_buffer.rows(), image_buffer.cols(),
                     cv_data_type, raw_data_ptr, step_size);

  if (input_image.channels() != 2) {
    // If there is no mask on the input image, use an empty mask.
    cv_mask = cv::Mat();
    return;
  }

  // Just make sure the masked input image pixels are zero to create the opencv mask
  ImageView<vw::uint8> mask_buffer = channel_cast_rescale<uint8>(select_channel(input_buffer, 1));

  // Wrap our mask object with OpenCV
  raw_data_ptr = reinterpret_cast<void*>(mask_buffer.data());
  cv::Mat full_cv_mask(mask_buffer.rows(), mask_buffer.cols(),
                       cv_data_type, raw_data_ptr, step_size);

  // Use OpenCV call to erode some pixels off the mask
  const int ERODE_RADIUS = 5;
  cv::Mat kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(ERODE_RADIUS, ERODE_RADIUS));
  cv::erode(full_cv_mask, cv_mask, kernel);

  return;
}

// TODO(oalexan1): Move to .cc
inline
void copy_opencv_descriptor_matrix(InterestPointList::iterator begin, 
                                   InterestPointList::iterator end,
                                   cv::Mat const& cvDescriptors, 
                                   OpenCvIpDetectorType detector_type) {

  size_t num_ip_descriptors = cvDescriptors.rows;
  size_t descriptor_length  = cvDescriptors.cols;
  // Copy the data to the output iterator
  // - Each IP needs the descriptor (a vector of floats) updated
  size_t ip_index = 0;
  for (auto iter = begin; iter != end; ++iter) {

    if (ip_index == num_ip_descriptors)
      vw_throw(LogicErr() << "copy_opencv_descriptor_matrix: IP sizes do not match!\n");

    // OpenCV descriptors can be of varying types, but InterestPoint only stores floats.
    // The workaround is to convert each element to a float here, and then convert back to
    //   the correct type when matching is performed.

    // TODO: Make sure all descriptor types work here!
    iter->descriptor.set_size(descriptor_length);
    switch (detector_type) {
      case OPENCV_IP_DETECTOR_TYPE_BRISK:
      case OPENCV_IP_DETECTOR_TYPE_ORB:
        for (size_t d=0; d<descriptor_length; ++d)
          iter->descriptor[d] 
            = static_cast<float>(cvDescriptors.at<unsigned char>(ip_index, d));
        break;
      case OPENCV_IP_DETECTOR_TYPE_SIFT:
        for (size_t d=0; d<descriptor_length; ++d)
          iter->descriptor[d] 
            = static_cast<float>(cvDescriptors.at<float>(ip_index, d))/512.0f;
        break;
      case OPENCV_IP_DETECTOR_TYPE_SURF:
        for (size_t d=0; d<descriptor_length; ++d)
          iter->descriptor[d] 
            = cvDescriptors.at<float>(ip_index, d); // TODO: May be incorrect!
        break;
      default: vw_throw(ArgumentErr() << "Unrecognized OpenCV detector type!\n");
    };
    ++ip_index;
  }
}


inline
cv::Ptr<cv::FeatureDetector>
OpenCvInterestPointDetector::init_detector(OpenCvIpDetectorType detector_type,
                                           int max_points) const {

  // Ensure the --threads option and / or the .vwrc thread setting are respected.
  cv::setNumThreads(vw_settings().default_num_threads());
 
  // Instantiate the feature detector
  // - There are a lot of detector variables that we just leave as the default here.
  switch (detector_type) {
  case OPENCV_IP_DETECTOR_TYPE_BRISK:
    vw_throw(NoImplErr() << "OpenCV BRISK option is not supported yet!\n");
    //return cv::Ptr<cv::BRISK>(new cv::BRISK());  break;
  case OPENCV_IP_DETECTOR_TYPE_ORB:
    if (max_points > 0)
      return cv::ORB::create(max_points);
    else // Unlike SIFT this will not work with a zero points input value
      return cv::ORB::create(100);
  case OPENCV_IP_DETECTOR_TYPE_SIFT:
    return cv::SIFT::create(max_points);
  case OPENCV_IP_DETECTOR_TYPE_SURF:
    vw_throw(NoImplErr() << "OpenCV SURF option is not supported yet!\n");
    //m_detector = cv::Ptr<cv::xfeatures2d::SURF >(new cv::SURF());  break;
  default:
    vw_throw(ArgumentErr() << "Unrecognized OpenCV detector type!\n");
  };
}

inline
OpenCvInterestPointDetector::OpenCvInterestPointDetector(OpenCvIpDetectorType detector_type,
                                                         bool normalize, 
                                                         bool add_descriptions, 
                                                         int max_points): 
m_detector_type(detector_type), m_add_descriptions(add_descriptions),
m_normalize(normalize), m_max_points(max_points) {

  m_detector = init_detector(detector_type, max_points);
}

/// Detect interest points in the source image.
template <class ViewT>
InterestPointList
OpenCvInterestPointDetector::process_image(ImageViewBase<ViewT> const& image,
                                           int desired_num_ip) const {
  // If the image is too small to use, don't return any interest points.
  const int MIN_DETECTOR_SIZE = 32;
  if ((image.impl().cols() < MIN_DETECTOR_SIZE) || (image.impl().rows() < MIN_DETECTOR_SIZE))
    return InterestPointList();

  // Convert the image into a plain uint8 image buffer wrapped by OpenCV
  std::cout << "--now in OpenCvInterestPointDetector::process_image\n";
  ImageView<vw::uint8> image_buffer;
  cv::Mat cv_image, cv_mask;
  try {
    get_opencv_wrapper(image, cv_image, image_buffer, cv_mask, m_normalize);
  } catch(vw::LogicErr) {
    return InterestPointList();
  }

  // Detect features
  std::vector<cv::KeyPoint> keypoints;
  cv::Mat cvDescriptors;

  // If the number of desired points has changed we need to create a temporary detector instance
  cv::Ptr<cv::FeatureDetector> detector = m_detector;
  if ((desired_num_ip > 0) && (desired_num_ip != m_max_points))
    detector = init_detector(m_detector_type, desired_num_ip);

  // Perform detection/description
  if (m_add_descriptions)
    detector->detectAndCompute (cv_image, cv_mask, keypoints, cvDescriptors);
  else // Don't add descriptions
    detector->detect(cv_image, keypoints, cv_mask);

  // Convert back to our output format
  InterestPointList ip_list;
  for (size_t i=0; i<keypoints.size(); ++i) {
    InterestPoint ip;
    ip.setFromCvKeypoint(keypoints[i]);
    ip_list.push_back(ip);
  }

  if (m_add_descriptions)
    copy_opencv_descriptor_matrix(ip_list.begin(), ip_list.end(), cvDescriptors, m_detector_type);

  return ip_list;
}


}} // namespace vw::ip

#endif // End case with OpenCV installed
#endif // __VW_INTEREST_POINT_OPENCV_DETECTOR_H__
