// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

/// \file OrbitingPusbroomModel.h
/// 
/// Orbiting pushbroom camera model object. 
///
/// You can then access the camera matrices using calls to
/// get_camera_matrix(), etc. Due to numerical instability in the
/// linear pushbroom model, it is also sometimes necessary to access
/// the individual components of the total linear pushbroom camera
/// matrix.  These are the
///
///  - Rigid transformation matrix
///  - Velocity matrix
///  - Perspective projection matrix
///
/// These are composed into the complete LP matrix as follows:
///
///  LP = Perpective * Velocity * RigidTransformation
///
/// The camera model will assume a straight flight path based on the 
/// position and velocity vectors.  Orientation is constant throughout
/// the flight path as well.
///
/// Notes:
///  
/// The orbiting pushbroom camera model is based on the linear
/// pushbroom model, but it takes the local curvature of a spacecraft's
/// orbit into account when computing the camera matrix.  These
/// positions are based on a time series of data that is supplied in
/// the variables
/// 
/// camera_poses (std::vector of quaternions)
/// positions    (std::vector of Vector3's containing a position in 3-space)
/// velocity_vectors  (std::vector of Vector3's, velocity vector in 3-space)
///
/// The velocity, and orientation will be approximated using a 2nd
/// degree quadratic curve, fit using linear least squares.  The
/// orientation (quaternion) will be extrapolated using the SLERP
/// algorithm.  (See CurveFitting.cc and SLERP.cc)

#ifndef _VW_CAMERA_ORBITINGPUSHBROOM_MODEL_H_
#define _VW_CAMERA_ORBITINGPUSHBROOM_MODEL_H_

#include <vw/Camera/LinearPushbroomModel.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>

#include <vector>

namespace vw { 
namespace camera {
  
  class OrbitingPushbroomModel : public LinearPushbroomModel {

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    OrbitingPushbroomModel() {}
    OrbitingPushbroomModel( double scan_duration,
                            double number_of_lines, 
                            double samples_per_line, 
                            double sample_offset, 
                            double focal_length, 
                            double along_scan_pixel_size,
                            double across_scan_pixel_size,
                            double t0_camera_pose,
                            double dt_camera_pose,
                            double t0_position,
                            double dt_position,
                            std::vector<Quaternion<double> > const& camera_poses,
                            std::vector<Vector3> const& positions,
                            std::vector<Vector3> const& velocity_vectors);

    virtual ~OrbitingPushbroomModel() {}

    //------------------------------------------------------------------
    // Interface
    //------------------------------------------------------------------

    virtual Vector2 vector_to_pixel(Vector3 const& vec) const;

    /// Given a pixel in image coordinates, what is the pointing
    /// vector in 3-space if you apply the camera model.
    ///
    /// Important Note: For linear pushbroom sensors, the orientation
    /// of the image is important.  In the Camera module we adopt the
    /// following convention for pushbroom imagers:
    ///
    /// u is the crosstrack (i.e., perspective) direction
    /// v is the downtrack (i.e., orthographic) direction
    virtual Vector3 pixel_to_vector(Vector2 const& pix) const;
    virtual Vector3 camera_center(Vector2 const& pix = Vector2() ) const;

    //------------------------------------------------------------------
    // Public Methods
    //------------------------------------------------------------------
    
    /// Return the camera matrix that describes the complete set of
    /// intrinsic and extrinsic parameters.
    ///
    /// The camera matrix is recomputed only when parameters change,
    /// so you can call get_camera_matrix() repeatedly with little penalty.
    /// 
    /// IMPORTANT NOTE: 
    ///
    /// This camera model can be applied to a 3D point to project it into
    /// the camera image space with coordinates [u v 1].  HOWEVER, this
    /// camera model assumes that the u variable is the is the
    /// ORTHOGRAPHIC projection axis (the time axis of the linear
    /// pusbroom image) and v is the PERSPECTIVE projection (the axis
    /// perpendicular to the time axis).  Be sure that you have your u
    /// and v coordinates straight before you may need to switch your use
    /// this camera model.
    const Matrix<double>& get_camera_matrix(Vector2 const& pix = Vector2() );

    /// Alternatively, you can extract each of the component of the
    ///full camera matrix (projection, velocity, and rigid
    ///transformation. This is handy in some cases if you are having
    ///numerical stability problems, and you want to rearrange the
    ///order that use to multipy points by the camera matrix.
    void get_decomposed_camera_matrix(Matrix<double> &projection,
                                              Matrix<double> &velocity_matrix,
                                              Matrix<double> &similarity, 
                                              Vector2 const& pix = Vector2() );
    // Accessors 
    inline double t0_camera_pose() const { return m_t0_camera_pose; }
    inline double dt_camera_pose() const { return m_dt_camera_pose; }
    inline double t0_position() const { return m_t0_position; }
    inline double dt_position() const { return m_dt_position; }
    inline std::vector<Quaternion<double> > camera_poses() const { return m_camera_poses; } 
    inline std::vector<Vector3> positions() const { return m_positions; } 
    inline std::vector<Vector3> velocity_vectors() const { return m_velocity_vectors; }
    inline Vector3 const initial_position() const { return m_positions[0]; }
    inline Vector3 const velocity_vector() const { return m_velocity_vectors[0]; }
    inline Quaternion<double> const camera_pose() const { return m_camera_poses[0]; }
    

    void set_camera_poses(std::vector<Quaternion<double> > const& val, double t0, double dt) {  
      m_t0_camera_pose = t0;
      m_dt_camera_pose = dt;
      m_camera_poses = val; 
    }
    void set_positions( std::vector<Vector3> const& val, double t0, double dt);
    void set_velocity_vectors( std::vector<Vector3> const& val, double t0, double dt);

  private:
    
    //------------------------------------------------------------------
    // Private Variables
    //------------------------------------------------------------------      
    std::vector<Quaternion<double> > m_camera_poses;
    std::vector<Vector3> m_positions;
    std::vector<Vector3> m_velocity_vectors;

    // Time dependent camera model parameters
    double m_t0_camera_pose;
    double m_dt_camera_pose;
    double m_t0_position;
    double m_dt_position;

    // Time dependent camera model parameters (these are coefficients to 
    // curves that are fit the the spacecraft motion 
    Matrix<double,3,3> m_position_coeff;
    Matrix<double,3,3> m_velocity_vector_coeff;
  };  

  /// Output stream method for printing a summary of the linear
  /// pushbroom camera model parameters.
  std::ostream& operator<<( std::ostream& os, OrbitingPushbroomModel const& camera_model);

}}	// namespace vw::camera

#endif	//_VW_CAMERA_ORBITINGPUSHBROOM_MODEL_H_

