// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file Detector.h
/// 
/// Built-in classes and functions for performing interest point detection.
/// 
#ifndef __INTERESTPOINT_DETECTOR_H__
#define __INTERESTPOINT_DETECTOR_H__

#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/Extrema.h>
#include <vw/InterestPoint/Localize.h>
#include <vw/InterestPoint/Interest.h>
#include <vw/InterestPoint/Threshold.h>
#include <vw/InterestPoint/WeightedHistogram.h>
#include <vw/InterestPoint/ImageOctaveHistory.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Manipulation.h>
#include <vw/FileIO.h>

#include <vector>
#include <iostream>
#include <stdio.h>

namespace vw {
namespace ip {

// Lowe recommends 36 bins, but this could also be learned.
#define FEATURE_ORI_NBINS (36)

/// Abstract base class for all-in-one interest point detection
/// classes. The interest measure and thresholding method used
/// are passed in as parameters. Detection classes must
/// implement detect_points.
template <class T, class ThreshT = InterestThreshold<T> >
class InterestPointDetector {
public:
  InterestPointDetector(InterestBase<T> *i, KeypointThresholdBase<ThreshT> *t) :
    interest(i), thresh(t) {}
  virtual ~InterestPointDetector() {}

  /// Set the interest measure used.
  void set_interest_measure(InterestBase<T> *i) { interest = i; }

  /// Detect interest points in the source image.
  virtual std::vector<InterestPoint> detect_points(const ImageView<T> src) = 0;

  /// Write the intermediate images out to files.
  virtual int write_images() const = 0;

protected:
  InterestBase<T> *interest;
  KeypointThresholdBase<ThreshT> *thresh;

  //TODO: These abstract virtual functions are defined here to enforce
  // the decoupling of each step, but this may be too restrictive
  virtual int find_extrema(std::vector<InterestPoint>& points) = 0;

  virtual int localize(std::vector<InterestPoint>& points) = 0;

  virtual int threshold(std::vector<InterestPoint>& points) = 0;

  virtual int assign_orientations(std::vector<InterestPoint>& points) = 0;
};

//TODO: the interface for this free function needs to be standardized one
// way or another. The second version was added for compatibility with
// the Stereo Pipeline, but it breaks the first version due to
// templating ambiguities.

/// Find the keypoints in an image using the provided detector.  
/// 
/// Some images are too large to be processed for interest points all
/// at once.  If the user specifies a max_keypoint_image_dimension,
/// this value is used to segment the image into smaller images which
/// are passed individually to the keypoint detector.  This routine
/// combines the keypoints from the sub-images once detection is
/// complete.  Be aware that a few keypoints along the segment
/// borders may be lost.  A good max dimension depends on the amount
/// of RAM needed by the detector (and the total RAM available).  A
/// value of 2048 seems to work well in most cases.
template <class ViewT, class DetectorT>
std::vector<InterestPoint> interest_points(vw::ImageViewBase<ViewT> const& image, 
                                           DetectorT& detector,
                                           const int32 max_keypoint_image_dimension = 0) {
  std::vector<InterestPoint> interest_points;

  vw_out(InfoMessage) << "\tFinding interest points" << std::flush;

  // If the user has not specified a chunk size, we process the
  // entire image in one shot.
  if (!max_keypoint_image_dimension) {
    vw_out(InfoMessage) << "..." << std::flush;
    interest_points = detector.detect_points(image.impl());

  // Otherwise we segment the image and process each sub-image
  // individually.
  } else {
    
    std::vector<BBox2i> bboxes = image_blocks(image.impl(), max_keypoint_image_dimension, max_keypoint_image_dimension);
    for (int i = 0; i < bboxes.size(); ++i) {
      vw_out(InfoMessage) << "." << std::flush;
      
      std::vector<InterestPoint> new_interest_points;
      new_interest_points = detector.detect_points(crop(image.impl(), bboxes[i]));
      for (int n = 0; n < new_interest_points.size(); ++n) {
        new_interest_points[n].x += bboxes[i].min().x();
        new_interest_points[n].y += bboxes[i].min().y();
        interest_points.push_back(new_interest_points[n]);
      }
    }

  }
  vw_out(InfoMessage) << " done.";
  vw_out(InfoMessage) << "     (" << interest_points.size() << " keypoints found)\n";
  return interest_points;

}

// /// Stereo Pipeline compatibility version.
// template <class ViewT, class DetectorT>
// std::vector<InterestPoint> interest_points(vw::ImageViewBase<ViewT> const& image, 
//                                            DetectorT const& detector,
//                                            int32 max_keypoint_image_dimension = 0) {

//   std::vector<InterestPoint> interest_points;

//   vw_out(InfoMessage) << "\tFinding interest points" << std::flush;

//   // If the user has not specified a chunk size, we process the
//   // entire image in one shot.
//   if (!max_keypoint_image_dimension) {
//     vw_out(InfoMessage) << "..." << std::flush;
//     interest_points = detector(image.impl());

//   // Otherwise we segment the image and process each sub-image
//   // individually.
//   } else {
    
//     std::vector<BBox2i> bboxes = image_blocks(image.impl(), max_keypoint_image_dimension, max_keypoint_image_dimension);
//     for (int i = 0; i < bboxes.size(); ++i) {
//       vw_out(InfoMessage) << "." << std::flush;
      
//       std::vector<InterestPoint> new_interest_points;
//       new_interest_points = detector(crop(image.impl(), bboxes[i]));
//       for (int n = 0; n < new_interest_points.size(); ++n) {
//         new_interest_points[n].x += bboxes[i].min().x();
//         new_interest_points[n].y += bboxes[i].min().y();
//         interest_points.push_back(new_interest_points[n]);
//       }
//     }

//   }
//   vw_out(InfoMessage) << " done.";
//   vw_out(InfoMessage) << "     (" << interest_points.size() << " keypoints found)\n";
//   return interest_points;

// }


/// Get the orientation of the point at (i0,j0,k0).  This is done by
/// computing a weighted histogram of edge orientations in a region
/// around the detected point.  The weights for the weighted histogram
/// are computed by multiplying the edge magnitude at each point by a
/// gaussian weight.  The edge orientation histogram is then smoothed,
/// effectively computing a kernel density estimate.  This density
/// function is then searched for local peaks.
template <class T>
int get_orientation( std::vector<float>& orientation,
		     const ImageInterestData<T>& data,
		     int i0, int j0, float sigma_ratio = 1.0) {
  orientation.clear();
  // Nominal feature support patch is 41x41 at the base scale, and
  // we multiply by sigma[k]/sigma[1] for other planes.
  
  // Get bounds for scaled 41x41 window centered at (i,j) in plane k
  int halfwidth = (int)(20*sigma_ratio + 0.5);
  int left  = i0 - halfwidth;
  int top   = j0 - halfwidth;
  int width = halfwidth*2+1;
  if ( (left>=0) && (top>=0) &&
       (left+width<(int)(data.ori.cols())) &&
       (top+width<(int)(data.ori.rows()))) {
    vw::vw_out(DebugMessage) << "Computing histogram on " << width
			     << " " << width << " starting at "
			     << left << " " << top << std::endl;
    // Get cropped view of interest point support region
    ImageView<T> region_ori = vw::crop(data.ori,left,top,width,width);
    ImageView<T> region_mag = vw::crop(data.mag,left,top,width,width);
    
    // Compute (gaussian weight)*(edge magnitude) kernel
    ImageView<float> weight(width,width);
    make_gaussian_kernel_2d( weight, 6 * sigma_ratio, width );
    for (int32 j=0; j<region_mag.rows(); j++){
      for (int32 i=0; i<region_mag.cols(); i++){
        weight(i,j) *= region_mag(i,j);
      }
    }
    
    // Compute weighted histogram of edge orientations
    std::vector<double> histo;
    weighted_histogram( region_ori, weight, histo, 
                        -M_PI, M_PI, FEATURE_ORI_NBINS );
    
    // Smooth histogram
    smooth_weighted_histogram( histo, 5.0 );

    // Find modes
    std::vector<int> mode;
    find_weighted_histogram_mode( histo, mode );
    for (unsigned m=0; m<mode.size(); ++m)
      orientation.push_back( mode[m]*(2*M_PI/FEATURE_ORI_NBINS)-M_PI );
  }
  
  return 0;
}

/// This class performs interest point detection on a source image
/// without using scale space methods.
template <class T, class ThreshT = InterestThreshold<T> >
class SimpleInterestPointDetector : public InterestPointDetector<T> {
 public:
  SimpleInterestPointDetector(InterestBase<T> *i, KeypointThresholdBase<ThreshT> *t) :
    InterestPointDetector<T>(i, t) { }

  /// Detect interest points in the source image.
  virtual std::vector<InterestPoint> detect_points(const ImageView<T> src) {
    // Calculate gradients, orientations and magnitudes
    img_data.set_source(src);

    // Compute interest image
    this->interest->compute_interest(img_data);

    // Find extrema in interest image
    std::vector<InterestPoint> points;
    find_extrema(points);

    // Subpixel localization
    localize(points);

    // Threshold (after localization)
    threshold(points);

    // Assign orientations
    assign_orientations(points);

    // Return vector of interest points
    return points;
  }

  /// This method dumps the various images internal to the detector out
  /// to files for visualization and debugging.  The images written out
  /// are the x and y gradients, edge orientation and magnitude, and
  /// interest function values for the source image.
  virtual int write_images() const {
    // Save the X gradient
    ImageView<float> grad_x_image = normalize(img_data.grad_x);
    vw::write_image("grad_x.jpg", grad_x_image);

    // Save the Y gradient      
    ImageView<float> grad_y_image = normalize(img_data.grad_y);
    vw::write_image("grad_y.jpg", grad_y_image);

    // Save the edge orientation image
    ImageView<float> ori_image = normalize(img_data.ori);
    vw::write_image("ori.jpg", ori_image);

    // Save the edge magnitude image
    ImageView<float> mag_image = normalize(img_data.mag);
    vw::write_image("mag.jpg", mag_image);

    // Save the interest function image
    ImageView<float> interest_image = normalize(img_data.interest);
    vw::write_image("interest.jpg", interest_image);

    return 0;
  }

 protected:
  ImageInterestData<T> img_data;

  // By default, uses find_peaks in Extrema.h
  virtual int find_extrema(std::vector<InterestPoint>& points) {
    return find_peaks(points, img_data.interest,
		      this->interest->peak_type());
  }

  // Use fit_peak in Localize.h
  virtual int localize(std::vector<InterestPoint>& points) {
    for (int i = 0; i < points.size(); i++) {
      fit_peak(img_data.interest, points[i]);
    }

    return 0;
  }

  virtual int threshold(std::vector<InterestPoint>& points) {
    std::vector<InterestPoint>::iterator pos = points.begin();
    while (pos != points.end()) {
      if (!this->thresh->keep(*pos, img_data))
        pos = points.erase(pos);
      else
        pos++;
    }
  }

  virtual int assign_orientations(std::vector<InterestPoint>& points) {
    std::vector<float> orientation;
    std::vector<InterestPoint> tmp;
    for (int i = 0; i < points.size(); i++) {
      get_orientation(orientation, img_data, (int)(points[i].x + 0.5),
		      (int)(points[i].y + 0.5));
      for (int j = 0; j < orientation.size(); j++) {
	InterestPoint pt = points[i];
	pt.orientation = orientation[j];
	tmp.push_back(pt);
      }
    }
    points = tmp;
  }

};

/// This class performs interest point detection on a source image
/// making use of scale space methods to achieve scale invariance.
/// This assumes that the detector works properly over different
/// choices of scale.
template <class T, class ThreshT = InterestThreshold<T> >
class ScaledInterestPointDetector : public InterestPointDetector<T, ThreshT> {
 public:
  ScaledInterestPointDetector(InterestBase<T> *i, KeypointThresholdBase<ThreshT> *t) :
    InterestPointDetector<T>(i, t), num_scales(5), num_octaves(3), octave(NULL) {}

  ScaledInterestPointDetector(InterestBase<T> *i, KeypointThresholdBase<ThreshT> *t,
                              int scales, int octaves) :
    InterestPointDetector<T>(i, t), num_scales(scales), num_octaves(octaves), octave(NULL) {}

  virtual ~ScaledInterestPointDetector() {}

  /// Detect interest points in the source image.
  virtual std::vector<InterestPoint> detect_points(const ImageView<T> src) {
    //create scale space
    octave = new ImageOctave<T>(src, num_scales);
    std::vector<InterestPoint> points;
    img_data.resize(octave->num_planes);
    next_point = 0;

    for (int o = 0; o < num_octaves; o++) {
      // Calculate gradients, orientations and magnitudes
      //printf("Calculating image data\n");
      for (int k = 0; k < octave->num_planes; k++) {
        img_data[k].set_source(octave->scales[k]);
      }

      // Compute interest images
      for (int k = 0; k < octave->num_planes; k++) {
        this->interest->compute_interest(img_data[k], octave->plane_index_to_scale(k));
      }
      
      if (history) history->add_octave(img_data);

      // Find extrema in interest image
      find_extrema(points);

      // Subpixel localization
      localize(points);

      // Threshold
      threshold(points);

      // Assign orientations
      assign_orientations(points);

      // Scale subpixel location to move back to original coords
      for (; next_point < points.size(); next_point++) {
        points[next_point].x *= octave->base_scale;
        points[next_point].y *= octave->base_scale;
      }
      
      if (o != num_octaves - 1) {
        octave->build_next();
      }
    }

    delete octave;

    //return vector of interest points
    return points;
  }

  /// Record data over the entire scale space pyramid.
  void record_history(ImageOctaveHistory<ImageInterestData<T> > *h) {
    history = h;
  }

  /// This method dumps the various images internal to the detector out
  /// to files for visualization and debugging.  The images written out
  /// are the x and y gradients, edge orientation and magnitude, and
  /// interest function values for all planes in the octave processed.
  virtual int write_images() const
  {
    for (int k=0; k<img_data.size(); k++){
      //int imagenum = (int)(log((float)base_scale)/log(2.0)) * num_planes + k;
      int imagenum = k;
      char fname[256];

      // Save the scale
      sprintf( fname, "scale_%02d.jpg", imagenum );
      ImageView<float> scale_image = normalize(img_data[k].src);
      vw::write_image(fname, scale_image);
      
      // Save the X gradient
      sprintf( fname, "grad_x_%02d.jpg", imagenum );
      ImageView<float> grad_x_image = normalize(img_data[k].grad_x);
      vw::write_image(fname, grad_x_image);

      // Save the Y gradient      
      sprintf( fname, "grad_y_%02d.jpg", imagenum );
      ImageView<float> grad_y_image = normalize(img_data[k].grad_y);
      vw::write_image(fname, grad_y_image);

      // Save the edge orientation image
      sprintf( fname, "ori_%02d.jpg", imagenum );
      ImageView<float> ori_image = normalize(img_data[k].ori);
      vw::write_image(fname, ori_image);

      // Save the edge magnitude image
      sprintf( fname, "mag_%02d.jpg", imagenum );
      ImageView<float> mag_image = normalize(img_data[k].mag);
      vw::write_image(fname, mag_image);

      // Save the interest function image
      sprintf( fname, "interest_%02d.jpg", imagenum );
      ImageView<float> interest_image = normalize(img_data[k].interest);
      vw::write_image(fname, interest_image);
    }
    return 0;
  }

 protected:
  std::vector<ImageInterestData<T> > img_data;
  ImageOctave<T> *octave;
  ImageOctaveHistory<ImageInterestData<T> > *history;
  int num_scales, num_octaves;
  int next_point;

  // By default, uses find_peaks in Extrema.h
  virtual int find_extrema(std::vector<InterestPoint>& points) {
    return find_peaks(points, img_data, *octave,
		      this->interest->peak_type());
  }

  // By default, uses fit_peak in Localize.h
  virtual int localize(std::vector<InterestPoint>& points) {
    for (int i = next_point; i < points.size(); i++) {
      fit_peak(img_data, points[i], *octave);
    }

    return 0;
  }

  virtual int threshold(std::vector<InterestPoint>& points) {
    std::vector<InterestPoint>::iterator pos = points.begin();
    while (pos != points.end()) {
      int k = octave->scale_to_plane_index(pos->scale);
      if (!this->thresh->keep(*pos, img_data[k]))
        pos = points.erase(pos);
      else
        pos++;
    }
  }

  virtual int assign_orientations(std::vector<InterestPoint>& points) {
    std::vector<float> orientation;
    int size = points.size();

    for (int i = next_point; i < size; i++) {
      int k = octave->scale_to_plane_index(points[i].scale);
      get_orientation(orientation, img_data[k], (int)(points[i].x + 0.5),
		      (int)(points[i].y + 0.5), octave->sigma[k]/octave->sigma[1]);
      if (orientation.size() == 0) {
        vw::vw_out(DebugMessage) << "Cannot find orientation of interest point\n";
      } else {
        points[i].orientation = orientation[0];
        for (int j = 1; j < orientation.size(); j++) {
          InterestPoint pt = points[i];
          pt.orientation = orientation[j];
          points.push_back(pt);
        }
      }
    }
  }

};


}} // namespace vw::ip 

#endif // __INTERESTPOINT_DETECTOR_H__
