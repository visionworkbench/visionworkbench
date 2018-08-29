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


/// \file Detector.h
///
/// Built-in classes and functions for performing interest point detection.
///
/// The key function here is detect_interest_points()
///

#ifndef __VW_INTERESTPOINT_DETECTOR_H__
#define __VW_INTERESTPOINT_DETECTOR_H__

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

#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"
#endif

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
    inline ImplT      & impl()       { return static_cast<ImplT      &>(*this); }
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
  /// method for determing orientations -- for example, if there is
  /// more than one dominant edge orientation in an image, this will
  /// produce a blended average of those two directions, and if those
  /// two directions are opposites, then they will cancel each other
  /// out.  However, this seems to work well enough for the time being.
  template <class OriT, class MagT>
  float get_orientation( ImageViewBase<OriT> const& x_grad,
                         ImageViewBase<MagT> const& y_grad,
                         float i0, float j0, float sigma_ratio = 1.0);


  /// This class performs interest point detection on a source image
  /// without using scale space methods.
  /// - The input class here is an "interest" operator that returns an image with pixels
  ///   scored according to some interest metric.
  /// - After the interest image is generated, a set of fixed functions extract local maxima to
  ///   obtain the output interest points.
  template <class InterestT>
  class InterestPointDetector : public InterestDetectorBase<InterestPointDetector<InterestT> > {
  public:

    /// Setting max_points = 0 will disable interest point culling.
    /// Otherwies, the max_points most "interesting" points are returned.
    InterestPointDetector(int max_points = 1000) : m_interest(InterestT()), m_max_points(max_points) {}

    /// Setting max_points = 0 will disable interest point culling.
    /// Otherwies, the max_points most "interesting" points are returned.
    InterestPointDetector(InterestT const& interest, int max_points = 1000) : m_interest(interest), m_max_points(max_points) {}

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
  class ScaledInterestPointDetector : public InterestDetectorBase<ScaledInterestPointDetector<InterestT> >,
				      private boost::noncopyable {

  public:
    // TODO: choose number of octaves based on image size
    static const int IP_DEFAULT_SCALES  = 3;
    static const int IP_DEFAULT_OCTAVES = 3;

    /// Setting max_points = 0 will disable interest point culling.
    /// Otherwies, the max_points most "interesting" points are returned.
    ScaledInterestPointDetector(int max_points = 1000)
      : m_interest(InterestT()), m_scales(IP_DEFAULT_SCALES),
        m_octaves(IP_DEFAULT_OCTAVES), m_max_points(max_points) {}

    ScaledInterestPointDetector(InterestT const& interest, int max_points = 1000)
      : m_interest(interest), m_scales(IP_DEFAULT_SCALES),
        m_octaves(IP_DEFAULT_OCTAVES), m_max_points(max_points) {}

    ScaledInterestPointDetector(InterestT const& interest, int scales, int octaves, int max_points = 1000)
      : m_interest(interest), m_scales(scales), m_octaves(octaves), m_max_points(max_points) {}

    /// Detect interest points in the source image.
    template <class ViewT>
    InterestPointList process_image(ImageViewBase<ViewT> const& image, int desired_num_ip=0) const;

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




  // TODO: Make an OpenCV interface file for this stuff

  enum OpenCvIpDetectorType {OPENCV_IP_DETECTOR_TYPE_BRISK = 0,
                             OPENCV_IP_DETECTOR_TYPE_ORB   = 1,
                             OPENCV_IP_DETECTOR_TYPE_SIFT  = 2,
                             OPENCV_IP_DETECTOR_TYPE_SURF  = 3};

#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1

  /// Struct to convert a basic type to a single channel OpenCV type
  template <typename T> struct GetOpenCvPixelType                 { static const int type=CV_8UC1;  };
  template <>           struct GetOpenCvPixelType<short         > { static const int type=CV_16SC1; };
  template <>           struct GetOpenCvPixelType<unsigned short> { static const int type=CV_16UC1; };
  template <>           struct GetOpenCvPixelType<int           > { static const int type=CV_32SC1; };
  template <>           struct GetOpenCvPixelType<float         > { static const int type=CV_32FC1; };
  template <>           struct GetOpenCvPixelType<double        > { static const int type=CV_64FC1; };

  /// Get an OpenCV wrapper, rasterizing the VW image to a provided buffer.
  template <class ViewT>
  void get_opencv_wrapper(ImageViewBase<ViewT> const& input_image,
                          cv::Mat & cv_image,
                          ImageView<vw::uint8> &image_buffer,
                          cv::Mat & cv_mask,
                          bool normalize = true);

  template <class LIST_ITER>
  void copy_opencv_descriptor_matrix(LIST_ITER begin, LIST_ITER end,
                                     cv::Mat const& cvDescriptors, OpenCvIpDetectorType detector_type);

  /// Interest point detector build using OpenCV functions
  class OpenCvInterestPointDetector : public InterestDetectorBase<OpenCvInterestPointDetector> {
  public:

    OpenCvInterestPointDetector(OpenCvIpDetectorType detector_type = OPENCV_IP_DETECTOR_TYPE_SIFT,
                                bool normalize=true,
                                bool add_descriptions=false, int max_points = 1000);

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
#else
  // If OpenCV is not defined, make a dummy version of this that just fails.

  /// Interest point detector build using OpenCV functions
  class OpenCvInterestPointDetector : public InterestDetectorBase<OpenCvInterestPointDetector> {
  public:
    OpenCvInterestPointDetector(OpenCvIpDetectorType detector_type = OPENCV_IP_DETECTOR_TYPE_BRISK,
                                bool normalize=true,
                                bool add_descriptions=false, int max_points = 1000){
      vw_throw( ArgumentErr() << "Can't use OpenCV IP detection functions if VW is not built with OpenCV!\n");
    }

    /// Detect interest points in the source image.
    template <class ViewT>
    InterestPointList process_image(ImageViewBase<ViewT> const& image, int desired_num_ip=0) const {
      vw_throw( ArgumentErr() << "Can't use OpenCV IP detection functions if VW is not built with OpenCV!\n");
      return InterestPointList();
    }
  }; // End class OpenCvInterestPointDetector

#endif // End case with no OpenCV installed



  // -----------------------------------------------------------------------------
  // The next set of classes are for performing IP detection with thread pools.
  // TODO: Maybe they should be moved to another file.

  /// This task class is used to insure that interest points are
  /// written to their list in a repeatable order that is not effected
  /// by the order in which the detection threads start and finish.
  class InterestPointWriteTask : public Task, private boost::noncopyable {
    InterestPointList  m_points;
    InterestPointList& m_global_points;

  public:
    InterestPointWriteTask( InterestPointList  local_points,
                            InterestPointList& global_points ) :
      m_points(local_points), m_global_points(global_points) {}

    virtual ~InterestPointWriteTask(){}

    virtual void operator() () {
      m_global_points.splice( m_global_points.end(), m_points );
    }
  };

  /// IP task wrapper for use with the InterestDetectionQueue thread pool class.
  /// - After IPs are found, they are passed to an InterestPointWriteTask object.
  template <class ViewT, class DetectorT>
  class InterestPointDetectionTask : public Task, private boost::noncopyable {

    ViewT              m_view;           ///< Source image
    DetectorT        & m_detector;       ///< Interest point detection class instance (TODO: const?)
    BBox2i             m_bbox;           ///< Region of the source image to check for points
    int                m_desired_num_ip; 
    int                m_id, m_max_id;
    InterestPointList& m_global_points;
    OrderedWorkQueue & m_write_queue;

  public:
    InterestPointDetectionTask(ImageViewBase<ViewT> const& view,
                               DetectorT& detector, BBox2i const& bbox, 
                               int desired_num_ip, int id, int max_id,
                               InterestPointList& global_list, OrderedWorkQueue& write_queue) :
      m_view(view.impl()), m_detector(detector), m_bbox(bbox), 
      m_desired_num_ip(desired_num_ip), m_id(id), m_max_id(max_id),
      m_global_points(global_list), m_write_queue(write_queue) {}

    virtual ~InterestPointDetectionTask(){}

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
  class InterestDetectionQueue : public WorkQueue {
    ViewT               m_view;
    DetectorT         & m_detector;
    OrderedWorkQueue  & m_write_queue;
    InterestPointList & m_ip_list;
    std::vector<BBox2i> m_bboxes;
    int                 m_tile_size;
    int                 m_desired_num_ip;
    Mutex               m_mutex;
    size_t              m_index;

    typedef InterestPointDetectionTask<ViewT, DetectorT> task_type;

  public:

    InterestDetectionQueue( ImageViewBase<ViewT> const& view, DetectorT& detector,
                            OrderedWorkQueue& write_queue, InterestPointList& ip_list,
                            int tile_size, int desired_num_ip=0 );

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
  InterestPointList detect_interest_points(ImageViewBase<ViewT> const& view, DetectorT& detector,
                                           int desired_num_ip=0);


// Include all the function definitions
#include <vw/InterestPoint/Detector.tcc>


}} // namespace vw::ip

#endif // __VW_INTERESTPOINT_DETECTOR_H__
