// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

}}      // namespace vw::camera

#endif  //__VW_CAMERA_LINEARPUSHBROOM_MODEL_H__
