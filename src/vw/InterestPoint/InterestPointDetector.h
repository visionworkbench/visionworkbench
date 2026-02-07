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

/// \file InterestPointDetector.h
///
/// Built-in classes and functions for performing interest point detection.
/// The detectors in this file are not used by default in ASP but are available
/// in ipfind with the log and harris detector options.

/// The key function here is detect_interest_points()

#ifndef __VW_INTEREST_POINT_INTEREST_POINT_DETECTOR_H__
#define __VW_INTEREST_POINT_INTEREST_POINT_DETECTOR_H__

#include <vw/Core/Debugging.h>
#include <vw/Core/Settings.h>
#include <vw/Core/ThreadPool.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Filter.h>

#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/ImageInterestData.h>
#include <vw/InterestPoint/Extrema.h>
#include <vw/InterestPoint/Localize.h>
#include <vw/InterestPoint/InterestOperator.h>
#include <vw/InterestPoint/WeightedHistogram.h>
#include <vw/InterestPoint/ImageOctaveHistory.h>

namespace vw {
namespace ip {

  /// This class performs interest point detection on a source image
  /// without using scale space methods.
  /// - The input class here is an "interest" operator that returns an image with pixels
  ///   scored according to some interest metric.
  /// - After the interest image is generated, a set of fixed functions extract local maxima to
  ///   obtain the output interest points.
  template <class InterestT>
  class InterestPointDetector: public InterestDetectorBase<InterestPointDetector<InterestT>> {
  public:

    /// Setting max_points = 0 will disable interest point culling.
    /// Otherwise, the max_points most "interesting" points are returned.
    InterestPointDetector(int max_points = 1000): m_interest(InterestT()), m_max_points(max_points) {}

    /// Setting max_points = 0 will disable interest point culling.
    /// Otherwise, the max_points most "interesting" points are returned.
    InterestPointDetector(InterestT const& interest, int max_points = 1000): m_interest(interest), m_max_points(max_points) {}

    /// Detect interest points in the source image.
    InterestPointList process_image(vw::ImageViewRef<float> const& image, int desired_num_ip=0) const;

  protected:
    InterestT m_interest;
    int       m_max_points;

    // By default, use find_peaks in Extrema.h
    template <class DataT>
    inline int find_extrema(InterestPointList& points, DataT const& img_data) const {
      return find_peaks(points, img_data);
    }

    // By default, use fit_peak in Localize.h
    template <class DataT>
    inline void localize(InterestPointList& points, DataT const& img_data) const;

    template <class DataT>
    inline void threshold(InterestPointList& points, DataT const& img_data) const;

    template <class DataT>
    void assign_orientations(InterestPointList& points, DataT const& img_data) const;

    /// This method dumps the various images internal to the detector out
    /// to files for visualization and debugging.  The images written out
    /// are the x and y gradients, edge orientation and magnitude, and
    /// interest function values for the source image.
    template <class DataT>
    void write_images(DataT const& img_data) const;
  };

  /// This class performs interest point detection on a source image
  /// making use of scale space methods to achieve scale invariance.
  /// This assumes that the detector works properly over different choices of scale.
  template <class InterestT>
  class ScaledInterestPointDetector: public InterestDetectorBase<ScaledInterestPointDetector<InterestT>>,
                      private boost::noncopyable {

  public:
    // TODO: choose number of octaves based on image size
    static const int IP_DEFAULT_SCALES  = 3;
    static const int IP_DEFAULT_OCTAVES = 3;

    /// Setting max_points = 0 will disable interest point culling.
    /// Otherwise, the max_points most "interesting" points are returned.
    ScaledInterestPointDetector(int max_points = 1000):
    m_interest(InterestT()), m_scales(IP_DEFAULT_SCALES),
    m_octaves(IP_DEFAULT_OCTAVES), m_max_points(max_points) {}

    ScaledInterestPointDetector(InterestT const& interest, int max_points = 1000):
    m_interest(interest), m_scales(IP_DEFAULT_SCALES),
    m_octaves(IP_DEFAULT_OCTAVES), m_max_points(max_points) {}

    ScaledInterestPointDetector(InterestT const& interest, int scales, int octaves,
                                int max_points = 1000):
    m_interest(interest), m_scales(scales), m_octaves(octaves), m_max_points(max_points) {}

    /// Detect interest points in the source image.
    InterestPointList process_image(vw::ImageViewRef<float> const& image,
                                    int desired_num_ip=0) const;

  protected:
    InterestT m_interest;
    int m_scales, m_octaves, m_max_points;

    // By default, uses find_peaks in Extrema.h
    template <class DataT, class ViewT>
    inline int find_extrema(InterestPointList& points,
                            std::vector<DataT> const& img_data,
                            ImageOctave<ViewT> const& octave) const {
      return find_peaks(points, img_data, octave);
    }

    // By default, uses fit_peak in Localize.h
    template <class DataT, class ViewT>
    inline int localize(InterestPointList& points,
                        std::vector<DataT> const& img_data,
                        ImageOctave<ViewT> const& octave) const;

    template <class DataT, class ViewT>
    inline void threshold(InterestPointList& points,
                          std::vector<DataT> const& img_data,
                          ImageOctave<ViewT> const& octave) const;

    template <class DataT, class ViewT>
    void assign_orientations(InterestPointList& points,
                             std::vector<DataT> const& img_data,
                             ImageOctave<ViewT> const& octave) const;

    /// This method dumps the various images internal to the detector out
    /// to files for visualization and debugging.  The images written out
    /// are the x and y gradients, edge orientation and magnitude, and
    /// interest function values for all planes in the octave processed.
    template <class DataT>
    int write_images(std::vector<DataT> const& img_data) const;
  }; // End class ScaledInterestPointDetector

//-------------------------------------------------------------------
// InterestPointDetector

// Detect interest points in the source image.
template <class InterestT>
InterestPointList
InterestPointDetector<InterestT>::process_image(vw::ImageViewRef<float> const& image,
                                                int desired_num_ip) const {

  // Calculate gradients, orientations and magnitudes

  // We blur the image by a small amount to eliminate any noise
  // that might confuse the interest point detector.
  // - First rasterize the input image to a buffer and convert to required float type
  ImageView<PixelGray<float>> float_image = pixel_cast_rescale<PixelGray<float>>(image);
  typedef SeparableConvolutionView<ImageView<PixelGray<float>>,
    typename DefaultKernelT<PixelGray<float>>::type, ConstantEdgeExtension>
      blurred_image_type;
  ImageInterestData<blurred_image_type, InterestT> img_data(gaussian_filter(float_image,0.5));

  // Compute interest image
  m_interest(img_data);

  // Find extrema in interest image
  InterestPointList points;
  find_extrema(points, img_data);

  // Subpixel localization
  localize(points, img_data);

  // Threshold (after localization)
  threshold(points, img_data);

  // Handle max_points override
  int curr_max_points = m_max_points;
  if (desired_num_ip > 0)
    curr_max_points = desired_num_ip;

  // Cull (limit the number of interest points to the N "most interesting" points)
  int original_num_points = points.size();
  points.sort();
  if ((curr_max_points > 0) && (curr_max_points < original_num_points))
     points.resize(curr_max_points);

  // Assign orientations
  assign_orientations(points, img_data);

  // Return vector of interest points
  return points;
}

// By default, use fit_peak in Localize.h
template <class InterestT>
template <class DataT>
void InterestPointDetector<InterestT>::localize(InterestPointList& points, DataT const& img_data) const {
  // TODO: Remove points rejected by localizer
  for (InterestPointList::iterator i = points.begin(); i != points.end(); ++i) {
    fit_peak(img_data.interest(), *i);
  }
}

template <class InterestT>
template <class DataT>
void InterestPointDetector<InterestT>::threshold(InterestPointList& points,
                                                 DataT const& img_data) const {
  // TODO: list::remove_if
  InterestPointList::iterator pos = points.begin();
  while (pos != points.end()) {
    if (!m_interest.threshold(*pos, img_data))
      pos = points.erase(pos);
    else
      pos++;
  }
}

template <class InterestT>
template <class DataT>
void InterestPointDetector<InterestT>::assign_orientations(InterestPointList& points,
                                                           DataT const& img_data) const {

  for (InterestPointList::iterator i = points.begin(); i != points.end(); ++i) {
    i->orientation = get_orientation(pixel_cast<float>(img_data.gradient_x()), pixel_cast<float>(img_data.gradient_y()), i->x, i->y);
  }
}

/// This method dumps the various images internal to the detector out
/// to files for visualization and debugging.  The images written out
/// are the x and y gradients, edge orientation and magnitude, and
/// interest function values for the source image.
template <class InterestT>
template <class DataT>
void InterestPointDetector<InterestT>::write_images(DataT const& img_data) const {
  // Save the X gradient
  ImageView<float> grad_x_image = normalize(img_data.gradient_x());
  vw::write_image("grad_x.jpg", grad_x_image);

  // Save the Y gradient
  ImageView<float> grad_y_image = normalize(img_data.gradient_y());
  vw::write_image("grad_y.jpg", grad_y_image);

  // Save the edge orientation image
  ImageView<float> ori_image = normalize(img_data.orientation());
  vw::write_image("ori.jpg", ori_image);

  // Save the edge magnitude image
  ImageView<float> mag_image = normalize(img_data.magnitude());
  vw::write_image("mag.jpg", mag_image);

  // Save the interest function image
  ImageView<float> interest_image = normalize(img_data.interest());
  vw::write_image("interest.jpg", interest_image);
}

//-------------------------------------------------------------------
// ScaledInterestPointDetector

// Detect interest points in the source image.
template <class InterestT>
InterestPointList ScaledInterestPointDetector<InterestT>::
process_image(vw::ImageViewRef<float> const& image, int desired_num_ip) const {

  // Create scale space
  typedef ImageInterestData<ImageView<PixelGray<float>>,InterestT> DataT;

  // We blur the image by a small amount to eliminate any noise
  // that might confuse the interest point detector.
  ImageOctave<typename DataT::source_type> octave(gaussian_filter(
                    pixel_cast_rescale<PixelGray<float>>(image),0.5), m_scales);

  InterestPointList points;
  std::vector<DataT> img_data;

  for (int o = 0; o < m_octaves; ++o) {

    // Calculate intermediate data (gradients, orientations, magnitudes)
    {
      img_data.clear();
      img_data.reserve(octave.num_planes);
      for (int k = 0; k < octave.num_planes; ++k) {
        img_data.push_back(DataT(octave.scales[k]));
      }
    }

    // Compute interest images
    for (int k = 0; k < octave.num_planes; ++k) {
      m_interest(img_data[k], octave.plane_index_to_scale(k));
    }

    // Find extrema in interest image
    InterestPointList new_points;
    find_extrema(new_points, img_data, octave);

    // Subpixel localization
    localize(new_points, img_data, octave);

    // Threshold
    threshold(new_points, img_data, octave);

    // Handle max_points override
    int curr_max_points = m_max_points;
    if (desired_num_ip > 0)
      curr_max_points = desired_num_ip;

    // Cull (limit the number of interest points to the N "most interesting" points)
    int original_num_points = new_points.size();
    new_points.sort();
    if ((curr_max_points > 0) && (curr_max_points/m_octaves < (int)new_points.size()))
      new_points.resize(curr_max_points/m_octaves);

    // Assign orientations
    assign_orientations(new_points, img_data, octave);

    // Scale subpixel location to move back to original coords
    for (InterestPointList::iterator i = new_points.begin(); i != new_points.end(); ++i) {
      i->x *= octave.base_scale;
      i->y *= octave.base_scale;
      // TODO: make sure this doesn't screw up any post-processing
      i->ix = (int)(i->x + 0.5);
      i->iy = (int)(i->y + 0.5);
    }

    // Add newly found interest points
    points.splice(points.end(), new_points);

    // Build next octave of scale space
    if (o != m_octaves - 1) {
      octave.build_next();
    }
  }

  return points;
}

// By default, uses fit_peak in Localize.h
template <class InterestT>
template <class DataT, class ViewT>
int ScaledInterestPointDetector<InterestT>::localize(InterestPointList& points,
            std::vector<DataT> const& img_data,
            ImageOctave<ViewT> const& octave) const {
  for (InterestPointList::iterator i = points.begin(); i != points.end(); ++i) {
    fit_peak(img_data, *i, octave);
  }

  return 0;
}

template <class InterestT>
template <class DataT, class ViewT>
void ScaledInterestPointDetector<InterestT>::threshold(InterestPointList& points,
                                                       std::vector<DataT> const& img_data,
                                                       ImageOctave<ViewT> const& octave) const {
  // TODO: list::remove_if
  InterestPointList::iterator pos = points.begin();
  while (pos != points.end()) {
    int k = octave.scale_to_plane_index(pos->scale);
    if (!m_interest.threshold(*pos, img_data[k]))
      pos = points.erase(pos);
    else
      pos++;
  }
}

template <class InterestT>
template <class DataT, class ViewT>
void ScaledInterestPointDetector<InterestT>::assign_orientations(InterestPointList& points,
                                                                 std::vector<DataT> const& img_data,
                                                                 ImageOctave<ViewT> const& octave) const {

  for (InterestPointList::iterator i = points.begin(); i != points.end(); ++i) {
    int k = octave.scale_to_plane_index(i->scale);
    i->orientation = get_orientation(pixel_cast<float>(img_data[k].gradient_x()), pixel_cast<float>(img_data[k].gradient_y()),
                                     i->x, i->y, octave.sigma[k]/octave.sigma[1]);
  }
}

// This method dumps the various images internal to the detector out
// to files for visualization and debugging.  The images written out
// are the x and y gradients, edge orientation and magnitude, and
// interest function values for all planes in the octave processed.
template <class InterestT>
template <class DataT>
int ScaledInterestPointDetector<InterestT>::write_images(std::vector<DataT> const& img_data) const {
  for (int k = 0; k < img_data.size(); ++k) {
    int imagenum = k;
    char fname[256];

    // Save the scale
    snprintf(fname, sizeof(fname), "scale_%02d.jpg", imagenum);
    ImageView<float> scale_image = normalize(img_data[k].source());
    vw::write_image(fname, scale_image);

    // Save the X gradient
    snprintf(fname, sizeof(fname), "grad_x_%02d.jpg", imagenum);
    ImageView<float> grad_x_image = normalize(img_data[k].gradient_x());
    vw::write_image(fname, grad_x_image);

    // Save the Y gradient
    snprintf(fname, sizeof(fname), "grad_y_%02d.jpg", imagenum);
    ImageView<float> grad_y_image = normalize(img_data[k].gradient_y());
    vw::write_image(fname, grad_y_image);

    // Save the edge orientation image
    snprintf(fname, sizeof(fname), "ori_%02d.jpg", imagenum);
    ImageView<float> ori_image = normalize(img_data[k].orientation());
    vw::write_image(fname, ori_image);

    // Save the edge magnitude image
    snprintf(fname, sizeof(fname), "mag_%02d.jpg", imagenum);
    ImageView<float> mag_image = normalize(img_data[k].magnitude());
    vw::write_image(fname, mag_image);

    // Save the interest function image
    snprintf(fname, sizeof(fname), "interest_%02d.jpg", imagenum);
    ImageView<float> interest_image = normalize(img_data[k].interest());
    vw::write_image(fname, interest_image);
  }
  return 0;
}

}} // namespace vw::ip

#endif // __VW_INTEREST_POINT_INTEREST_POINT_DETECTOR_H__
