// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file OrbitingPushbroomModel.h
/// 
/// Orbiting pushbroom camera model object. 
///
/// Notes:
///  
/// The orbiting pushbroom camera model takes the local curvature of a
/// spacecraft's orbit into account when computing it's position.  The
/// position estimate is based on:
/// 
/// camera_poses (std::vector of quaternions)
/// positions    (std::vector of Vector3's containing a position in 3-space)
///
/// The position is be approximated using a 2nd degree quadratic
/// curve, fit using linear least squares.  The orientation
/// (quaternion) will be extrapolated using the SLERP algorithm.  (See
/// Extrinsics.h)

#ifndef __VW_CAMERA_ORBITINGPUSHBROOM_MODEL_H__
#define __VW_CAMERA_ORBITINGPUSHBROOM_MODEL_H__

#include <vw/Camera/LinescanModel.h>
#include <vw/Camera/Extrinsics.h>

#include <vector>

namespace vw { 
namespace camera {
  
  class OrbitingPushbroomModel : public LinescanModel<Curve3DPositionInterpolation,
                                                      SLERPPoseInterpolation> {
    
  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------

    // Set up the position and pose estimators.  These are the
    // characteristics that differentiate the orbiting pushbroom model
    // from other linescan imagers.
    OrbitingPushbroomModel( int number_of_lines, 
                            int samples_per_line, 
                            int sample_offset, 
                            double focal_length, 
                            double along_scan_pixel_size,
                            double across_scan_pixel_size,
                            double line_integration_time,
                            double t0_camera_pose,
                            double dt_camera_pose,
                            double t0_position,
                            double dt_position,
                            Vector3 pointing_vec,
                            Vector3 u_vec,
                            std::vector<Quaternion<double> > const& camera_poses,
                            std::vector<Vector3> const& positions) : 
      LinescanModel<Curve3DPositionInterpolation,
                    SLERPPoseInterpolation>::LinescanModel(number_of_lines, 
                                                           samples_per_line, 
                                                           sample_offset,
                                                           focal_length, 
                                                           along_scan_pixel_size, 
                                                           across_scan_pixel_size,
                                                           line_integration_time,
                                                           pointing_vec, u_vec,
                                                           Curve3DPositionInterpolation(positions, t0_position, dt_position),
                                                           SLERPPoseInterpolation(camera_poses, t0_camera_pose, dt_camera_pose)) {}
    
    // This constructor is used when the exposure time varies from
    // scanline to scanline, e.g. in the case of HRSC.
    OrbitingPushbroomModel( int number_of_lines, 
                            int samples_per_line, 
                            int sample_offset, 
                            double focal_length, 
                            double along_scan_pixel_size,
                            double across_scan_pixel_size,
                            std::vector<double> line_times,
                            double t0_camera_pose,
                            double dt_camera_pose,
                            double t0_position,
                            double dt_position,
                            Vector3 pointing_vec,
                            Vector3 u_vec,
                            std::vector<Quaternion<double> > const& camera_poses,
                            std::vector<Vector3> const& positions) :
      LinescanModel<Curve3DPositionInterpolation,
                    SLERPPoseInterpolation>::LinescanModel(number_of_lines, 
                                                           samples_per_line, 
                                                           sample_offset,
                                                           focal_length, 
                                                           along_scan_pixel_size, 
                                                           across_scan_pixel_size,
                                                           line_times,
                                                           pointing_vec, u_vec,
                                                           Curve3DPositionInterpolation(positions, t0_position, dt_position),
                                                           SLERPPoseInterpolation(camera_poses, t0_camera_pose, dt_camera_pose)) {}

    virtual ~OrbitingPushbroomModel() {}
    virtual std::string type() const { return "OrbitingPushbroom"; }

  };  
}}	// namespace vw::camera

#endif	//__VW_CAMERA_ORBITINGPUSHBROOM_MODEL_H__

