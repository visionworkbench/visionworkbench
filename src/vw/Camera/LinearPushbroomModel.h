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

/// \file LinearPushbroomModel.h
/// 
/// Linear pushbroom camera model object. 
///
#ifndef __VW_CAMERA_LINEARPUSHBROOM_MODEL_H__
#define __VW_CAMERA_LINEARPUSHBROOM_MODEL_H__

#include <vw/Camera/LinescanModel.h>
#include <vw/Camera/Extrinsics.h>
#include <vw/Math/Quaternion.h>

namespace vw { 
namespace camera {
  
  class LinearPushbroomModel : public LinescanModel<LinearPositionInterpolation,
                                                    ConstantPoseInterpolation> {

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------

    // Set up the position and pose estimator.  These are the
    // characteristics that differentiate the linear pushbroom model
    // from other linescan imagers.
    LinearPushbroomModel( double scan_duration,
                          int number_of_lines, 
                          int samples_per_line, 
                          int sample_offset, 
                          double focal_length, 
                          double along_scan_pixel_size,
                          double across_scan_pixel_size,
                          Vector3 pointing_vec,
                          Vector3 u_vec,
                          Quaternion<double> const& camera_pose,
                          Vector3 const& initial_position,
                          Vector3 const& velocity_vector) : 
      LinescanModel<LinearPositionInterpolation,
                    ConstantPoseInterpolation>::LinescanModel(number_of_lines, 
                                                              samples_per_line, 
                                                              sample_offset,
                                                              focal_length, 
                                                              along_scan_pixel_size, 
                                                              across_scan_pixel_size,
                                                              scan_duration / number_of_lines,
                                                              pointing_vec, u_vec,
                                                              LinearPositionInterpolation(initial_position, velocity_vector),
                                                              ConstantPoseInterpolation(camera_pose)) {}

    virtual ~LinearPushbroomModel() {}
    virtual std::string type() const { return "LinearPushbroom"; }
  };  

}}	// namespace vw::camera

#endif	//__VW_CAMERA_LINEARPUSHBROOM_MODEL_H__

