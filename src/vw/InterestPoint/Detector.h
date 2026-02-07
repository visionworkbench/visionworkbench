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

/// \file Detector.h
///
/// Built-in classes and functions for performing interest point detection.
///
/// The key function here is detect_interest_points()
///

#ifndef __VW_INTEREST_POINT_DETECTOR_H__
#define __VW_INTEREST_POINT_DETECTOR_H__

#include <vw/Core/Debugging.h>
#include <vw/Core/Settings.h>
#include <vw/Core/ThreadPool.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Filter.h>

#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/Extrema.h>
#include <vw/InterestPoint/Localize.h>
#include <vw/InterestPoint/InterestOperator.h>
#include <vw/InterestPoint/WeightedHistogram.h>
#include <vw/InterestPoint/ImageOctaveHistory.h>

namespace vw {
namespace ip {

  /// A CRTP InterestDetector base class.
  ///
  /// Subclasses of InterestDetectorBase must provide an
  /// implementation that detects interest points in an image region.
  /// The call operator() above calls this method, passing in each
  /// individual patch to be processed. That declaration is expected
  /// to be of the form:
  ///
  ///    template <class ViewT>
  ///    InterestPointList process_image(ImageViewBase<ViewT> const& image) const {
  ///
  /// This class is CRTP in order to get around the inability to have a
  ///  template <class ViewT> virtual function.
  ///
  template <class ImplT>
  struct InterestDetectorBase {

    // Methods to access the derived type
    inline ImplT      & impl() { return static_cast<ImplT      &>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }

    /// Find the interest points in an image using the provided detector.
    /// - desired_num_ip is ignored if set to zero.
    template <class ViewT>
    InterestPointList operator() (vw::ImageViewBase<ViewT> const& image,
                                  int desired_num_ip=0);
  };

  /// Get the orientation of the point at (i0,j0,k0).  This is done by
  /// computing a gaussian weighted average orientations in a region
  /// around the detected point.  This is not the most sophisticated
  /// method for determining orientations -- for example, if there is
  /// more than one dominant edge orientation in an image, this will
  /// produce a blended average of those two directions, and if those
  /// two directions are opposites, then they will cancel each other
  /// out.  However, this seems to work well enough for the time being.
  template <class OriT, class MagT>
  float get_orientation(ImageViewBase<OriT> const& x_grad,
                         ImageViewBase<MagT> const& y_grad,
                         float i0, float j0, float sigma_ratio = 1.0);

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
    template <class ViewT>
    InterestPointList process_image(ImageViewBase<ViewT> const& image, int desired_num_ip=0) const;

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
    template <class ViewT>
    InterestPointList process_image(ImageViewBase<ViewT> const& image,
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

  // The next set of classes are for performing IP detection with thread pools.
  // TODO: Maybe they should be moved to another file.

  /// This task class is used to insure that interest points are
  /// written to their list in a repeatable order that is not affected
  /// by the order in which the detection threads start and finish.
  class InterestPointWriteTask: public Task, private boost::noncopyable {
    InterestPointList  m_points;
    InterestPointList& m_global_points;
    vw::TerminalProgressCallback & m_tpc;
    double m_inc_amt;

  public:
    InterestPointWriteTask(InterestPointList  local_points,
                           InterestPointList& global_points,
                           vw::TerminalProgressCallback & tpc,
                           double inc_amt):
      m_points(local_points), m_global_points(global_points),
      m_tpc(tpc), m_inc_amt(inc_amt) {}

    virtual ~InterestPointWriteTask() {}

    virtual void operator() () {
      m_global_points.splice(m_global_points.end(), m_points);
      m_tpc.report_incremental_progress(m_inc_amt);
    }
  };

  /// IP task wrapper for use with the InterestDetectionQueue thread pool class.
  /// - After IPs are found, they are passed to an InterestPointWriteTask object.
  template <class ViewT, class DetectorT>
  class InterestPointDetectionTask: public Task, private boost::noncopyable {

    ViewT              m_view;           ///< Source image
    DetectorT        & m_detector;       ///< Interest point detection class instance (TODO: const?)
    BBox2i             m_bbox;           ///< Region of the source image to check for points
    int                m_desired_num_ip;
    int                m_job_id, m_num_jobs;
    InterestPointList& m_global_points;
    OrderedWorkQueue & m_write_queue;
    vw::TerminalProgressCallback & m_tpc;
    double m_inc_amt;

  public:
    InterestPointDetectionTask(ImageViewBase<ViewT> const& view,
                               DetectorT& detector, BBox2i const& bbox,
                               int desired_num_ip, int id, int num_jobs,
                               InterestPointList& global_list,
                               OrderedWorkQueue& write_queue,
                               vw::TerminalProgressCallback & tpc):
      m_view(view.impl()), m_detector(detector), m_bbox(bbox),
      m_desired_num_ip(desired_num_ip), m_job_id(id), m_num_jobs(num_jobs),
      m_global_points(global_list), m_write_queue(write_queue),
      m_tpc(tpc) {
        m_inc_amt = 1.0 / double(m_num_jobs);
      }

    virtual ~InterestPointDetectionTask() {}

    /// Find the IPs, then pass to an InterestPointWriteTask instance.
    void operator()();

    InterestPointList interest_point_list() { return m_global_points; }
  };

  /// A thread pool class for performing interest point detection.
  /// - There is a lot of memory allocation created on task generation. I
  ///   couldn't figure it out in a reasonable time frame. Thus now we
  ///   generate tasks on demand which should lower the instantaneous
  ///   memory requirement.
  template <class ViewT, class DetectorT>
  class InterestDetectionQueue: public WorkQueue {
    ViewT               m_view;
    DetectorT         & m_detector;
    OrderedWorkQueue  & m_write_queue;
    InterestPointList & m_ip_list;
    std::vector<BBox2i> m_bboxes;
    int                 m_tile_size;
    int                 m_desired_num_ip;
    Mutex               m_mutex;
    size_t              m_index;
    vw::TerminalProgressCallback & m_tpc;

    typedef InterestPointDetectionTask<ViewT, DetectorT> task_type;

  public:

    InterestDetectionQueue(ImageViewBase<ViewT> const& view, DetectorT& detector,
                           OrderedWorkQueue& write_queue, InterestPointList& ip_list,
                           int tile_size, int desired_num_ip,
                           vw::TerminalProgressCallback & tpc);

    size_t size() { return m_bboxes.size(); }

    virtual boost::shared_ptr<Task> get_next_task();
  };

  // End thread pool class declarations.
  // -----------------------------------------------------------------------------

  /// This function implements a multithreaded interest point detector using
  /// the passed-in detector object.
  /// - Threads are spun off to process the image in 1024x1024 pixel blocks.
  /// - Pass in desired_num_ip to enforce this limit proportional to the tile size,
  ///   otherwise each tile will use the same number regardless of size.
  template <class ViewT, class DetectorT>
  InterestPointList detect_interest_points(ImageViewBase<ViewT> const& view,
                                           DetectorT& detector,
                                           int desired_num_ip=0);

// Function definitions

//-------------------------------------------------------------------
// InterestDetectorBase

// Find the interest points in an image using the provided detector.
template <class ImplT>
template <class ViewT>
InterestPointList
InterestDetectorBase<ImplT>::operator()(vw::ImageViewBase<ViewT> const& image,
                                        int desired_num_ip) {

  InterestPointList interest_points;

  // It is up to the individual IP implementations to convert the input image
  //  to their desired input format.
  interest_points = impl().process_image(image.impl(), desired_num_ip);

  return interest_points;
}

//-------------------------------------------------------------------
// InterestPointDetectionTask

template <class ViewT, class DetectorT>
void InterestPointDetectionTask<ViewT, DetectorT>::operator()() {

  vw_out(InfoMessage, "interest_point")
    << "Locating interest points in block " << m_job_id + 1 << "/" << m_num_jobs << "   [ "
    << m_bbox << " ] with " << m_desired_num_ip << " ip.\n";

  // Use the m_detector object to find a set of image points in the cropped
  // section of the image.
  // This is a lazy crop. The actual cropping and bringing into memory will happen
  // in each detector implementation.
  InterestPointList new_ip_list = m_detector(crop(vw::pixel_cast<float>(m_view.impl()), m_bbox),
                                             m_desired_num_ip);

  for (InterestPointList::iterator pt = new_ip_list.begin(); pt != new_ip_list.end(); ++pt) {
    (*pt).x  += m_bbox.min().x();
    (*pt).ix += m_bbox.min().x();
    (*pt).y  += m_bbox.min().y();
    (*pt).iy += m_bbox.min().y();
  }

  // Append these interest points to the master list owned by the
  // detect_interest_points() function. It appears that m_write_queue
  // is accessed by only one thread at a time, so this is thread-safe.
  boost::shared_ptr<Task>
    write_task(new InterestPointWriteTask(new_ip_list, m_global_points,
                                          m_tpc, m_inc_amt));
  m_write_queue.add_task(write_task, m_job_id);

  vw_out(InfoMessage, "interest_point")
    << "Finished block " << m_job_id + 1 << "/" << m_num_jobs << std::endl;
}

//-------------------------------------------------------------------
// InterestDetectionQueue

template <class ViewT, class DetectorT>
InterestDetectionQueue<ViewT, DetectorT>::
InterestDetectionQueue(ImageViewBase<ViewT> const& view, DetectorT& detector,
                       OrderedWorkQueue& write_queue, InterestPointList& ip_list,
                       int tile_size, int desired_num_ip,
                       vw::TerminalProgressCallback & tpc):
     m_view(view.impl()), m_detector(detector),
     m_write_queue(write_queue), m_ip_list(ip_list), m_tile_size(tile_size),
     m_desired_num_ip(desired_num_ip), m_index(0), m_tpc(tpc) {

  m_bboxes = subdivide_bbox(m_view, tile_size, tile_size);
  this->notify();
}

template <class ViewT, class DetectorT>
boost::shared_ptr<Task>
InterestDetectionQueue<ViewT, DetectorT>::get_next_task() {
  Mutex::Lock lock(m_mutex);
  if (m_index == m_bboxes.size())
    return boost::shared_ptr<Task>();

  m_index++;

  int num_ip = 0; // The default value means let the detector pick the IP count.

  if (m_desired_num_ip > 0) {
    // Determine the desired number of IP for this tile based on its size
    //  relative to a full sized tile.
    const int MIN_NUM_IP = 1;
    double expected_area = m_tile_size*m_tile_size;
    double area = m_bboxes[m_index-1].area();
    double fraction = area / expected_area;
    num_ip = ceil(fraction * static_cast<double>(m_desired_num_ip));
    if (num_ip < MIN_NUM_IP)
      num_ip = MIN_NUM_IP;
    if (num_ip > m_desired_num_ip)
      num_ip = m_desired_num_ip;
  }

  return boost::shared_ptr<Task>(new task_type(m_view, m_detector,
                                               m_bboxes[m_index-1], num_ip, m_index-1,
                                               m_bboxes.size(), m_ip_list, m_write_queue,
                                               m_tpc));
}

//-------------------------------------------------------------------

// This free function implements a multithreaded interest point
// detector. Threads are spun off to process the image in 1024x1024 pixel blocks.
template <class ViewT, class DetectorT>
InterestPointList detect_interest_points(ImageViewBase<ViewT> const& view,
                                         DetectorT& detector,
                                         int desired_num_ip) {

  // Process the image in no less than 1024x1024 size pixel blocks.
  int tile_size = vw_settings().default_tile_size();
  if (tile_size < 1024)
    tile_size = 1024;

  // Progress dialog
  vw::TerminalProgressCallback tpc("asp", "Detect interest points: ");
  tpc.report_progress(0);

  // Create an ordered thread pool with one thread. Ensure that interest points
  // are written in a specific order and not by the random way threads finish.
  OrderedWorkQueue write_queue(1);
  InterestPointList ip_list;

  InterestDetectionQueue<ViewT, DetectorT>
    detect_queue(view, detector, write_queue, ip_list, tile_size, desired_num_ip, tpc);
  detect_queue.join_all();
  write_queue.join_all();
  tpc.report_finished();

  return ip_list;
}

//-------------------------------------------------------------------

// Get the orientation of the point at (i0,j0,k0).  This is done by
// computing a gaussian weighted average orientations in a region
// around the detected point.  This is not the most sophisticated
// method for determining orientations -- for example, if there is
// more than one dominant edge orientation in an image, this will
// produce a blended average of those two directions, and if those
// two directions are opposites, then they will cancel each other
// out.  However, this seems to work well enough for the time being.
template <class OriT, class MagT>
float get_orientation(ImageViewBase<OriT> const& x_grad,
                      ImageViewBase<MagT> const& y_grad,
                      float i0, float j0, float sigma_ratio) {

  // The size, in pixels, of the image patch used to compute the orientation.
  static const int IP_ORIENTATION_WIDTH = 10;

  // Nominal feature support patch is WxW at the base scale, with
  // W = IP_ORIENTATION_HALF_WIDTH * 2 + 1, and
  // we multiply by sigma[k]/sigma[1] for other planes.
  //
  // Get bounds for scaled WxW window centered at (i,j) in plane k
  int halfwidth = (int)(IP_ORIENTATION_WIDTH/2*sigma_ratio + 0.5);
  int left  = int(roundf(i0 - halfwidth));
  int top   = int(roundf(j0 - halfwidth));

  // Compute (gaussian weight)*(edge magnitude) kernel
  ImageView<float> weight(IP_ORIENTATION_WIDTH,IP_ORIENTATION_WIDTH);
  make_gaussian_kernel_2d(weight, 6 * sigma_ratio, IP_ORIENTATION_WIDTH);

  // We must compute the average orientation in quadrature.
  double weight_sum = sum_of_pixel_values(weight);

  // Compute the gaussian weighted average x_gradient
  ImageView<float> weighted_grad = weight * crop(edge_extend(x_grad.impl()),left,top,IP_ORIENTATION_WIDTH,IP_ORIENTATION_WIDTH);
  double avg_x_grad = sum_of_pixel_values(weighted_grad) / weight_sum;

  // Compute the gaussian weighted average y_gradient
  weighted_grad = weight * crop(edge_extend(y_grad.impl()),left,top,IP_ORIENTATION_WIDTH,IP_ORIENTATION_WIDTH);
  double avg_y_grad = sum_of_pixel_values(weighted_grad) / weight_sum;

  return atan2(avg_y_grad,avg_x_grad);
}

//-------------------------------------------------------------------
// InterestPointDetector

// Detect interest points in the source image.
template <class InterestT>
template <class ViewT>
InterestPointList
InterestPointDetector<InterestT>::process_image(ImageViewBase<ViewT> const& image,
                                                int desired_num_ip) const {

  // Calculate gradients, orientations and magnitudes

  // We blur the image by a small amount to eliminate any noise
  // that might confuse the interest point detector.
  // - First rasterize the input image to a buffer and convert to required float type
  ImageView<PixelGray<float>> float_image = pixel_cast_rescale<PixelGray<float>>(image);
  typedef SeparableConvolutionView<ImageView<PixelGray<float>>,
    typename DefaultKernelT<PixelGray<float> >::type, ConstantEdgeExtension>
      blurred_image_type;
  ImageInterestData<blurred_image_type, InterestT> img_data(gaussian_filter(float_image,0.5));

  // Compute interest image
  {
    m_interest(img_data);
  }

  // Find extrema in interest image
  InterestPointList points;
  {
    find_extrema(points, img_data);
  }

  // Subpixel localization
  {
    localize(points, img_data);
  }

  // Threshold (after localization)
  {
    threshold(points, img_data);
  }

  // Handle max_points override
  int curr_max_points = m_max_points;
  if (desired_num_ip > 0)
    curr_max_points = desired_num_ip;

  // Cull (limit the number of interest points to the N "most interesting" points)
  {
    int original_num_points = points.size();
    points.sort();
    if ((curr_max_points > 0) && (curr_max_points < original_num_points))
       points.resize(curr_max_points);
  }

  // Assign orientations
  {
    assign_orientations(points, img_data);
  }

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
    i->orientation = get_orientation(img_data.gradient_x(), img_data.gradient_y(), i->x, i->y);
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
template <class ViewT>
InterestPointList ScaledInterestPointDetector<InterestT>::
process_image(ImageViewBase<ViewT> const& image, int desired_num_ip) const {

  // Create scale space
  typedef ImageInterestData<ImageView<PixelGray<float> >,InterestT> DataT;

  // We blur the image by a small amount to eliminate any noise
  // that might confuse the interest point detector.
  ImageOctave<typename DataT::source_type> octave(gaussian_filter(
                    pixel_cast_rescale<PixelGray<float> >(image),0.5), m_scales);

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
    {
      for (int k = 0; k < octave.num_planes; ++k) {
        m_interest(img_data[k], octave.plane_index_to_scale(k));
      }
    }

    // Find extrema in interest image
    InterestPointList new_points;
    {
      find_extrema(new_points, img_data, octave);
    }

    // Subpixel localization
    {
      localize(new_points, img_data, octave);
    }

    // Threshold
    {
      threshold(new_points, img_data, octave);
    }

    // Handle max_points override
    int curr_max_points = m_max_points;
    if (desired_num_ip > 0)
      curr_max_points = desired_num_ip;

    // Cull (limit the number of interest points to the N "most interesting" points)
    {
      int original_num_points = new_points.size();
      new_points.sort();
      if ((curr_max_points > 0) && (curr_max_points/m_octaves < (int)new_points.size()))
        new_points.resize(curr_max_points/m_octaves);
    }

    // Assign orientations
    {
      assign_orientations(new_points, img_data, octave);
    }

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
    i->orientation = get_orientation(img_data[k].gradient_x(), img_data[k].gradient_y(),
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

#endif // __VW_INTEREST_POINT_DETECTOR_H__
