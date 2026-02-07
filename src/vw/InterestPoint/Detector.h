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
    InterestPointList operator() (vw::ImageViewRef<float> const& image,
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
  float get_orientation(vw::ImageViewRef<float> const& x_grad,
                        vw::ImageViewRef<float> const& y_grad,
                        float i0, float j0, float sigma_ratio = 1.0);

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
  template <class DetectorT>
  class InterestPointDetectionTask: public Task, private boost::noncopyable {

    vw::ImageViewRef<float> m_view;      ///< Source image
    DetectorT        & m_detector;       ///< Interest point detection class instance (TODO: const?)
    BBox2i             m_bbox;           ///< Region of the source image to check for points
    int                m_desired_num_ip;
    int                m_job_id, m_num_jobs;
    InterestPointList& m_global_points;
    OrderedWorkQueue & m_write_queue;
    vw::TerminalProgressCallback & m_tpc;
    double m_inc_amt;

  public:
    InterestPointDetectionTask(vw::ImageViewRef<float> const& view,
                               DetectorT& detector, BBox2i const& bbox,
                               int desired_num_ip, int id, int num_jobs,
                               InterestPointList& global_list,
                               OrderedWorkQueue& write_queue,
                               vw::TerminalProgressCallback & tpc):
      m_view(view), m_detector(detector), m_bbox(bbox),
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
  template <class DetectorT>
  class InterestDetectionQueue: public WorkQueue {
    vw::ImageViewRef<float> m_view;
    DetectorT         & m_detector;
    OrderedWorkQueue  & m_write_queue;
    InterestPointList & m_ip_list;
    std::vector<BBox2i> m_bboxes;
    int                 m_tile_size;
    int                 m_desired_num_ip;
    Mutex               m_mutex;
    size_t              m_index;
    vw::TerminalProgressCallback & m_tpc;

    typedef InterestPointDetectionTask<DetectorT> task_type;

  public:

    InterestDetectionQueue(vw::ImageViewRef<float> const& view, DetectorT& detector,
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
  template <class DetectorT>
  InterestPointList detect_interest_points(vw::ImageViewRef<float> const& view,
                                           DetectorT& detector,
                                           int desired_num_ip=0);

// Function definitions

//-------------------------------------------------------------------
// InterestDetectorBase

// Find the interest points in an image using the provided detector.
template <class ImplT>
InterestPointList
InterestDetectorBase<ImplT>::operator()(vw::ImageViewRef<float> const& image,
                                        int desired_num_ip) {

  InterestPointList interest_points;

  // It is up to the individual IP implementations to convert the input image
  //  to their desired input format.
  interest_points = impl().process_image(image.impl(), desired_num_ip);

  return interest_points;
}

//-------------------------------------------------------------------
// InterestPointDetectionTask

template <class DetectorT>
void InterestPointDetectionTask<DetectorT>::operator()() {

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

template <class DetectorT>
InterestDetectionQueue<DetectorT>::
InterestDetectionQueue(vw::ImageViewRef<float> const& view, DetectorT& detector,
                       OrderedWorkQueue& write_queue, InterestPointList& ip_list,
                       int tile_size, int desired_num_ip,
                       vw::TerminalProgressCallback & tpc):
     m_view(view), m_detector(detector),
     m_write_queue(write_queue), m_ip_list(ip_list), m_tile_size(tile_size),
     m_desired_num_ip(desired_num_ip), m_index(0), m_tpc(tpc) {

  m_bboxes = subdivide_bbox(m_view, tile_size, tile_size);
  this->notify();
}

template <class DetectorT>
boost::shared_ptr<Task>
InterestDetectionQueue<DetectorT>::get_next_task() {
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
template <class DetectorT>
InterestPointList detect_interest_points(vw::ImageViewRef<float> const& view,
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

  InterestDetectionQueue<DetectorT>
    detect_queue(view, detector, write_queue, ip_list, tile_size, desired_num_ip, tpc);
  detect_queue.join_all();
  write_queue.join_all();
  tpc.report_finished();

  return ip_list;
}

}} // namespace vw::ip

#endif // __VW_INTEREST_POINT_DETECTOR_H__
