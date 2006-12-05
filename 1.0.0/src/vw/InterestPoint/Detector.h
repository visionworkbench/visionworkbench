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

/// \file Image.h
/// 
/// Basic classes and structures for storing image interest points.
/// 
#ifndef __INTERESTPOINT_DETECTOR_H__
#define __INTERESTPOINT_DETECTOR_H__

#include <vw/InterestPoint/Descriptor.h>
#include <vw/Image/Algorithms.h>

namespace vw {
namespace ip {


  // Find the keypoints in an image using the provided detector.  
  // 
  // Some images are too large to be processed for interest points all
  // at once.  If the user specifies a max_keypoint_image_dimension,
  // this value is used to segment the image into smaller images which
  // are passed individually to the keypoint detector.  This routine
  // combines the keypoints from the sub-images once detection is
  // complete.  Be aware that a few keypoints along the segment
  // borders may be lost.  A good max dimension depends on the amount
  // of RAM needed by the detector (and the total RAM available).  A
  // value of 2048 seems to work well in most cases.
  template <class ViewT, class DetectorT>
  std::vector<InterestPoint> interest_points(vw::ImageViewBase<ViewT> const& image, 
                                             DetectorT const& detector,
                                             unsigned int max_keypoint_image_dimension = 0) {

    std::vector<InterestPoint> interest_points;

    vw_out(InfoMessage) << "\tFinding interest points" << std::flush;

    // If the user has not specified a chunk size, we process the
    // entire image in one shot.
    if (!max_keypoint_image_dimension) {
      vw_out(InfoMessage) << "..." << std::flush;
      interest_points = detector(image.impl());

    // Otherwise we segment the image and process each sub-image
    // individually.
    } else {
      
      std::vector<BBox2i> bboxes = image_blocks(image, max_keypoint_image_dimension, max_keypoint_image_dimension);
      for (int i = 0; i < bboxes.size(); ++i) {
        vw_out(InfoMessage) << "." << std::flush;
        
        std::vector<InterestPoint> new_interest_points;
        new_interest_points = detector(crop(image.impl(), bboxes[i]));
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

}} // namspace vw::ip 

#endif // __INTERESTPOINT_DETECTOR_H__
