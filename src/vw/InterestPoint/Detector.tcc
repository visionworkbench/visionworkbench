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


// Function definitions

//-------------------------------------------------------------------
// InterestDetectorBase

// Find the interest points in an image using the provided detector.
template <class ImplT>
template <class ViewT>
InterestPointList InterestDetectorBase<ImplT>::operator() (vw::ImageViewBase<ViewT> const& image) {

  InterestPointList interest_points;
  vw_out(DebugMessage, "interest_point") << "Finding interest points in block: [ " << image.impl().cols() 
                                         << " x " << image.impl().rows() << " ]\n";

  // Note: We explicitly convert the image to PixelGray<float>
  // (rescaling as necessary) here before passing it off to the
  // rest of the interest detector code.
  interest_points = impl().process_image(pixel_cast_rescale<PixelGray<float> >(image.impl()));

  vw_out(DebugMessage, "interest_point") << "Finished processing block. Found " << interest_points.size() 
                                         << " interest points.\n";
  return interest_points;
}



//-------------------------------------------------------------------
// InterestPointDetectionTask

template <class ViewT, class DetectorT>
void InterestPointDetectionTask<ViewT, DetectorT>::operator()() {

  vw_out(InfoMessage, "interest_point") << "Locating interest points in block "
                                        << m_id + 1 << "/" << m_max_id << "   [ " << m_bbox << " ]\n";

  // Use the m_detector object to find a set of image points in the cropped section of the image.
  InterestPointList new_ip_list = m_detector(crop(m_view.impl(), m_bbox));

  for (InterestPointList::iterator pt = new_ip_list.begin(); pt != new_ip_list.end(); ++pt) {
    (*pt).x  += m_bbox.min().x();
    (*pt).ix += m_bbox.min().x();
    (*pt).y  += m_bbox.min().y();
    (*pt).iy += m_bbox.min().y();
  }

  // Append these interest points to the master list
  // owned by the detect_interest_points() function.
  boost::shared_ptr<Task> write_task( new InterestPointWriteTask(new_ip_list, m_global_points) );
  m_write_queue.add_task( write_task, m_id );

  vw_out(InfoMessage, "interest_point") << "Finished block " << m_id + 1 << "/" << m_max_id << std::endl;
}


//-------------------------------------------------------------------
// InterestDetectionQueue

template <class ViewT, class DetectorT>
InterestDetectionQueue<ViewT, DetectorT>::
InterestDetectionQueue( ImageViewBase<ViewT> const& view, DetectorT& detector,
                        OrderedWorkQueue& write_queue, InterestPointList& ip_list,
                        int tile_size ) :
     m_view(view.impl()), m_detector(detector),
     m_write_queue(write_queue), m_ip_list(ip_list), m_index(0) {
  m_bboxes = subdivide_bbox( m_view, tile_size, tile_size );
  this->notify();
}


template <class ViewT, class DetectorT>
boost::shared_ptr<Task>
InterestDetectionQueue<ViewT, DetectorT>::get_next_task() {
  Mutex::Lock lock(m_mutex);
  if ( m_index == m_bboxes.size() )
    return boost::shared_ptr<Task>();

  m_index++;
  return boost::shared_ptr<Task>( new task_type( m_view, m_detector,
                                                 m_bboxes[m_index-1], m_index-1,
                                                 m_bboxes.size(), m_ip_list, m_write_queue ) 
                                );
}

//-------------------------------------------------------------------

// This free function implements a multithreaded interest point
// detector.  Threads are spun off to process the image in 1024x1024 pixel blocks.
template <class ViewT, class DetectorT>
InterestPointList detect_interest_points (ImageViewBase<ViewT> const& view, DetectorT& detector) {

  VW_OUT(DebugMessage, "interest_point") << "Running multi-threaded interest point detector.  Input image: [ "
                                         << view.impl().cols() << " x " << view.impl().rows() << " ]\n";

  // Process the image in no less than 1024x1024 size pixel blocks.
  int tile_size = vw_settings().default_tile_size();
  if (tile_size < 1024)
    tile_size = 1024;

  // Create an ordered thread pool with one thread.
  OrderedWorkQueue write_queue(1); // Used to insure that interest points are written in a
                                   // specific order and not by the random way threads finish.
  InterestPointList ip_list;
  InterestDetectionQueue<ViewT, DetectorT>
    detect_queue( view, detector, write_queue, ip_list, tile_size );
  VW_OUT(DebugMessage, "interest_point") << "Waiting for threads to terminate.\n";
  detect_queue.join_all();
  write_queue.join_all();
  VW_OUT(DebugMessage, "interest_point") << "MT interest point detection complete.  "
                                         << ip_list.size() << " interest point detected.\n";
  return ip_list;
}

//-------------------------------------------------------------------

// Get the orientation of the point at (i0,j0,k0).  This is done by
// computing a gaussian weighted average orientations in a region
// around the detected point.  This is not the most sophisticated
// method for determing orientations -- for example, if there is
// more than one dominant edge orientation in an image, this will
// produce a blended average of those two directions, and if those
// two directions are opposites, then they will cancel each other
// out.  However, this seems to work well enough for the time being.
template <class OriT, class MagT>
float get_orientation( ImageViewBase<OriT> const& x_grad,
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
  make_gaussian_kernel_2d( weight, 6 * sigma_ratio, IP_ORIENTATION_WIDTH );

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
InterestPointList InterestPointDetector<InterestT>::process_image(ImageViewBase<ViewT> const& image) const {
  Timer total("\tTotal elapsed time", DebugMessage, "interest_point");

  // Calculate gradients, orientations and magnitudes
  vw_out(DebugMessage, "interest_point") << "\n\tCreating image data... ";
  Timer *timer = new Timer("done, elapsed time", DebugMessage, "interest_point");

  // We blur the image by a small amount to eliminate any noise
  // that might confuse the interest point detector.
  typedef SeparableConvolutionView<ViewT, typename DefaultKernelT<typename ViewT::pixel_type>::type, 
                                                                  ConstantEdgeExtension> blurred_image_type;
  ImageInterestData<blurred_image_type, InterestT> img_data(gaussian_filter(image,0.5));

  delete timer;

  // Compute interest image
  vw_out(DebugMessage, "interest_point") << "\tComputing interest image... ";
  {
    Timer t("done, elapsed time", DebugMessage, "interest_point");
    m_interest(img_data);
  }

  // Find extrema in interest image
  vw_out(DebugMessage, "interest_point") << "\tFinding extrema... ";
  InterestPointList points;
  {
    Timer t("elapsed time", DebugMessage, "interest_point");
    find_extrema(points, img_data);
    vw_out(DebugMessage, "interest_point") << "done (" << points.size() << " interest points), ";
  }

  // Subpixel localization
  vw_out(DebugMessage, "interest_point") << "\tLocalizing... ";
  {
    Timer t("elapsed time", DebugMessage, "interest_point");
    localize(points, img_data);
    vw_out(DebugMessage, "interest_point") << "done (" << points.size() << " interest points), ";
  }

  // Threshold (after localization)
  vw_out(DebugMessage, "interest_point") << "\tThresholding... ";
  {
    Timer t("elapsed time", DebugMessage, "interest_point");
    threshold(points, img_data);
    vw_out(DebugMessage, "interest_point") << "done (" << points.size() << " interest points), ";
  }

  // Cull (limit the number of interest points to the N "most interesting" points)
  vw_out(DebugMessage, "interest_point") << "\tCulling Interest Points (limit is set to " << m_max_points << " points)... ";
  {
    Timer t("elapsed time", DebugMessage, "interest_point");
    int original_num_points = points.size();
    points.sort();
    if ((m_max_points > 0) && (m_max_points < original_num_points))
       points.resize(m_max_points);
    vw_out(DebugMessage, "interest_point") << "done (removed " << original_num_points - points.size() 
                                           << " interest points, " << points.size() << " remaining.), ";
  }

  // Assign orientations
  vw_out(DebugMessage, "interest_point") << "\tAssigning orientations... ";
  {
    Timer t("elapsed time", DebugMessage, "interest_point");
    assign_orientations(points, img_data);
    vw_out(DebugMessage, "interest_point") << "done (" << points.size() << " interest points), ";
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
void InterestPointDetector<InterestT>::threshold(InterestPointList& points, DataT const& img_data) const {
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
void InterestPointDetector<InterestT>::assign_orientations(InterestPointList& points, DataT const& img_data) const {

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
process_image(ImageViewBase<ViewT> const& image) const {
  typedef ImageInterestData<ImageView<typename ViewT::pixel_type>,InterestT> DataT;

  Timer total("\t\tTotal elapsed time", DebugMessage, "interest_point");

  // Create scale space
  vw_out(DebugMessage, "interest_point") << "\n\tCreating initial image octave... ";
  Timer *t_oct = new Timer("done, elapsed time", DebugMessage, "interest_point");

  // We blur the image by a small amount to eliminate any noise
  // that might confuse the interest point detector.
  ImageOctave<typename DataT::source_type> octave(gaussian_filter(image,0.5), m_scales);

  delete t_oct;

  InterestPointList points;
  std::vector<DataT> img_data;

  for (int o = 0; o < m_octaves; ++o) {
    Timer t_loop("\t\tElapsed time for octave", DebugMessage, "interest_point");

    // Calculate intermediate data (gradients, orientations, magnitudes)
    vw_out(DebugMessage, "interest_point") << "\tCreating image data... ";
    {
      Timer t("done, elapsed time", DebugMessage, "interest_point");
      img_data.clear();
      img_data.reserve(octave.num_planes);
      for (int k = 0; k < octave.num_planes; ++k) {
        img_data.push_back(DataT(octave.scales[k]));
      }
    }

    // Compute interest images
    vw_out(DebugMessage, "interest_point") << "\tComputing interest images... ";
    {
      Timer t("done, elapsed time", DebugMessage, "interest_point");
      for (int k = 0; k < octave.num_planes; ++k) {
        m_interest(img_data[k], octave.plane_index_to_scale(k));
      }
    }

    // Find extrema in interest image
    vw_out(DebugMessage, "interest_point") << "\tFinding extrema... ";
    InterestPointList new_points;
    {
      Timer t("elapsed time", DebugMessage, "interest_point");
      find_extrema(new_points, img_data, octave);
      vw_out(DebugMessage, "interest_point") << "done (" << new_points.size() << " extrema found), ";
    }

    // Subpixel localization
    vw_out(DebugMessage, "interest_point") << "\tLocalizing... ";
    {
      Timer t("elapsed time", DebugMessage, "interest_point");
      localize(new_points, img_data, octave);
      vw_out(DebugMessage, "interest_point") << "done (" << new_points.size() << " interest points), ";
    }

    // Threshold
    vw_out(DebugMessage, "interest_point") << "\tThresholding... ";
    {
      Timer t("elapsed time", DebugMessage, "interest_point");
      threshold(new_points, img_data, octave);
      vw_out(DebugMessage, "interest_point") << "done (" << new_points.size() << " interest points), ";
    }

    // Cull (limit the number of interest points to the N "most interesting" points)
    vw_out(DebugMessage, "interest_point") << "\tCulling Interest Points (limit is set to " << m_max_points << " points)... ";
    {
      Timer t("elapsed time", DebugMessage, "interest_point");
      int original_num_points = new_points.size();
      new_points.sort();
      if ((m_max_points > 0) && (m_max_points/m_octaves < (int)new_points.size()))
        new_points.resize(m_max_points/m_octaves);
      vw_out(DebugMessage, "interest_point") << "done (removed " << original_num_points - new_points.size() 
                                             << " interest points, " << new_points.size() << " remaining.), ";
    }

    // Assign orientations
    vw_out(DebugMessage, "interest_point") << "\tAssigning orientations... ";
    {
      Timer t("elapsed time", DebugMessage, "interest_point");
      assign_orientations(new_points, img_data, octave);
      vw_out(DebugMessage, "interest_point") << "done (" << new_points.size() << " interest points), ";
    }

    // Scale subpixel location to move back to original coords
    for (InterestPointList::iterator i = new_points.begin(); i != new_points.end(); ++i) {
      i->x *= octave.base_scale;
      i->y *= octave.base_scale;
      // TODO: make sure this doesn't screw up any post-processing
      i->ix = (int)( i->x + 0.5 );
      i->iy = (int)( i->y + 0.5 );
    }

    // Add newly found interest points
    points.splice(points.end(), new_points);

    // Build next octave of scale space
    if (o != m_octaves - 1) {
      vw_out(DebugMessage, "interest_point") << "\tBuilding next octave... ";
      Timer t("done, elapsed time", DebugMessage);
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
  for (int k = 0; k < img_data.size(); ++k){
    int imagenum = k;
    char fname[256];

    // Save the scale
    sprintf( fname, "scale_%02d.jpg", imagenum );
    ImageView<float> scale_image = normalize(img_data[k].source());
    vw::write_image(fname, scale_image);

    // Save the X gradient
    sprintf( fname, "grad_x_%02d.jpg", imagenum );
    ImageView<float> grad_x_image = normalize(img_data[k].gradient_x());
    vw::write_image(fname, grad_x_image);

    // Save the Y gradient
    sprintf( fname, "grad_y_%02d.jpg", imagenum );
    ImageView<float> grad_y_image = normalize(img_data[k].gradient_y());
    vw::write_image(fname, grad_y_image);

    // Save the edge orientation image
    sprintf( fname, "ori_%02d.jpg", imagenum );
    ImageView<float> ori_image = normalize(img_data[k].orientation());
    vw::write_image(fname, ori_image);

    // Save the edge magnitude image
    sprintf( fname, "mag_%02d.jpg", imagenum );
    ImageView<float> mag_image = normalize(img_data[k].magnitude());
    vw::write_image(fname, mag_image);

    // Save the interest function image
    sprintf( fname, "interest_%02d.jpg", imagenum );
    ImageView<float> interest_image = normalize(img_data[k].interest());
    vw::write_image(fname, interest_image);
  }
  return 0;
}



#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1

template <class ViewT>
cv::Mat get_opencv_wrapper(ImageViewBase<ViewT> const& input_image,
			   ImageView<PixelGray<vw::uint8> > &image_buffer,
			   bool normalize) {

  if (normalize) // Convert the input image to uint8 with 2%-98% intensity scaling.
    percentile_scale_convert(input_image, image_buffer, 0.02, 0.98);
  else {
    // Convert to uint8 using the default ranges for the input data type
    double standard_min = ChannelRange<typename ViewT::pixel_type>::min();
    double standard_max = ChannelRange<typename ViewT::pixel_type>::max();
    image_buffer = pixel_cast_rescale<vw::uint8>(clamp(input_image, standard_min, standard_max));
  }

  // Figure out the image buffer parameters
  int     cv_data_type = GetOpenCvPixelType<vw::uint8>::type;
  void*   raw_data_ptr = reinterpret_cast<void*>(image_buffer.data());
  size_t  pixel_size   = sizeof(vw::uint8);
  size_t  step_size    = image_buffer.cols() * pixel_size;

  // Create an OpenCV wrapper for the buffer image
  cv::Mat cv_image(image_buffer.rows(), image_buffer.cols(),
		   cv_data_type, raw_data_ptr, step_size);
  return cv_image;
}

template <class LIST_ITER>
void copy_opencv_descriptor_matrix(LIST_ITER begin, LIST_ITER end,
				   cv::Mat const& cvDescriptors, OpenCvIpDetectorType detector_type) {

  size_t num_ip_descriptors = cvDescriptors.rows;
  size_t descriptor_length  = cvDescriptors.cols;
  // Copy the data to the output iterator
  // - Each IP needs the descriptor (a vector of floats) updated
  size_t ip_index = 0;
  for (LIST_ITER iter = begin; iter != end; ++iter) {

    if (ip_index == num_ip_descriptors)
      vw_throw( LogicErr() << "copy_opencv_descriptor_matrix: IP sizes do not match!\n");

    // OpenCV descriptors can be of varying types, but InterstPoint only stores floats.
    // The workaround is to convert each element to a float here, and then convert back to
    //   the correct type when matching is performed.

    // TODO: Make sure all descriptor types work here!
    iter->descriptor.set_size(descriptor_length);
    switch (detector_type)
    {
      case OPENCV_IP_DETECTOR_TYPE_BRISK:
      case OPENCV_IP_DETECTOR_TYPE_ORB:
        for (size_t d=0; d<descriptor_length; ++d)
          iter->descriptor[d] = static_cast<float>(cvDescriptors.at<unsigned char>(ip_index, d));
        break;
      case OPENCV_IP_DETECTOR_TYPE_SIFT:
        for (size_t d=0; d<descriptor_length; ++d)
          iter->descriptor[d] = static_cast<float>(cvDescriptors.at<float>(ip_index, d))/512.0f;
        break;
      case OPENCV_IP_DETECTOR_TYPE_SURF:
        for (size_t d=0; d<descriptor_length; ++d)
          iter->descriptor[d] = cvDescriptors.at<float>(ip_index, d); // TODO: May be incorrect!
        break;
      default: vw_throw( ArgumentErr() << "Unrecognized OpenCV detector type!\n");
    };
    ++ip_index;
  }
}


inline
OpenCvInterestPointDetector::OpenCvInterestPointDetector(OpenCvIpDetectorType detector_type,
			    bool normalize,
			    bool add_descriptions, int max_points)
  : m_detector_type(detector_type), m_add_descriptions(add_descriptions), m_normalize(normalize) {

  // Instantiate the feature detector
  switch (detector_type)
  {
    case OPENCV_IP_DETECTOR_TYPE_BRISK:
      vw_throw( NoImplErr() << "OpenCV BRISK option is not supported yet!\n");
      //m_detector_brisk = cv::Ptr<cv::BRISK>(new cv::BRISK());  break;
    case OPENCV_IP_DETECTOR_TYPE_ORB:
      m_detector_orb = cv::ORB::create();
      m_detector_orb->setMaxFeatures(max_points);
      break;
    case OPENCV_IP_DETECTOR_TYPE_SIFT:
      m_detector_sift = cv::xfeatures2d::SIFT::create(max_points);
      break;
    case OPENCV_IP_DETECTOR_TYPE_SURF:
      vw_throw( NoImplErr() << "OpenCV SURF option is not supported yet!\n");
      //m_detector_surf  = cv::Ptr<cv::xfeatures2d::SURF >(new cv::SURF());  break;
    default: vw_throw( ArgumentErr() << "Unrecognized OpenCV detector type!\n");
  };
  //if (!m_detector) vw_throw( LogicErr() << "OpenCvInterestPointDetector: detector failed to initialize!\n");
}

/// Detect interest points in the source image.
template <class ViewT>
InterestPointList OpenCvInterestPointDetector::process_image(ImageViewBase<ViewT> const& image) const {

  //if (!m_detector)
  //    vw_throw( LogicErr() << "OpenCvInterestPointDetector: detector is not initialized!\n");

  // If the image is too small to use, don't return any interest points.
  const int MIN_DETECTOR_SIZE = 32;
  if ( (image.impl().cols() < MIN_DETECTOR_SIZE) || (image.impl().rows() < MIN_DETECTOR_SIZE))
    return InterestPointList();

  // Convert the image into a plain uint8 image buffer wrapped by OpenCV
  ImageView<PixelGray<vw::uint8> > buffer_image;
  cv::Mat cv_image = get_opencv_wrapper(image, buffer_image, m_normalize);

  // Detect features
  std::vector<cv::KeyPoint> keypoints;
  cv::Mat cvDescriptors;

  if (m_add_descriptions) {
    switch (m_detector_type)
    {
      //case OPENCV_IP_DETECTOR_TYPE_BRISK: m_detector_brisk->operator()(cv_image, cv::Mat(), keypoints, cvDescriptors); break;
      case OPENCV_IP_DETECTOR_TYPE_ORB:   m_detector_orb->detectAndCompute (cv_image, cv::Mat(), keypoints, cvDescriptors); break;
      case OPENCV_IP_DETECTOR_TYPE_SIFT:  m_detector_sift->detectAndCompute(cv_image, cv::Mat(), keypoints, cvDescriptors); break;
      //case OPENCV_IP_DETECTOR_TYPE_SURF:  m_detector_surf->operator() (cv_image, cv::Mat(), keypoints, cvDescriptors); break;
      default: vw_throw( ArgumentErr() << "Unrecognized OpenCV detector type!\n");
    };
  } else { // Don't add descriptions
    switch (m_detector_type)
    {
      //case OPENCV_IP_DETECTOR_TYPE_BRISK: m_detector_brisk->detect(cv_image, keypoints); break;
      case OPENCV_IP_DETECTOR_TYPE_ORB:   m_detector_orb->detect (cv_image, keypoints);   break;
      case OPENCV_IP_DETECTOR_TYPE_SIFT:  m_detector_sift->detect(cv_image, keypoints);  break;
      //case OPENCV_IP_DETECTOR_TYPE_SURF:  m_detector_surf->detect(cv_image, keypoints);  break;
      default: vw_throw( ArgumentErr() << "Unrecognized OpenCV detector type!\n");
    };
  }

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

#endif // End case with OpenCV installed




