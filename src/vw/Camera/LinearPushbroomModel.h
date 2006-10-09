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

/// \file LinearPusbroomModel.h
/// 
/// Linear pushbroom camera model object. 
///
/// You can then access the camera matrices using calls to
/// get_camera_matrix(), etc. Due to numerical instability in the
/// Linear pushbroom model, it is also sometimes necessary to access
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
#ifndef _VW_CAMERA_LINEARPUSHBROOM_MODEL_H_
#define _VW_CAMERA_LINEARPUSHBROOM_MODEL_H_

#include <vw/Camera/CameraModel.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Quaternion.h>

namespace vw { 
namespace camera {
  
  class LinearPushbroomModel : public CameraModel { 
  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    LinearPushbroomModel() {}
    LinearPushbroomModel( double scan_duration,
                          double number_of_lines, 
                          double samples_per_line, 
                          double sample_offset, 
                          double focal_length, 
                          double along_scan_pixel_size,
                          double across_scan_pixel_size,
                          Quaternion<double> const& camera_pose,
                          Vector3 const& initial_position,
                          Vector3 const& velocity_vector);

    virtual ~LinearPushbroomModel() {}
    
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
    const Matrix<double>& get_camera_matrix(Vector2 const& pix = Vector2() ) {
      return m_camera_matrix;
    }

    /// Alternatively, you can extract each of the component of the
    ///full camera matrix (projection, velocity, and rigid
    ///transformation. This is handy in some cases if you are having
    ///numerical stability problems, and you want to rearrange the
    ///order that use to multipy points by the camera matrix.
    virtual void get_decomposed_camera_matrix(Matrix<double> &projection,
                                              Matrix<double> &velocity_matrix,
                                              Matrix<double> &similarity, 
                                              Vector2 const& pix = Vector2() ) {
      projection = m_projection_matrix;
      velocity_matrix = m_velocity_matrix;
      similarity = m_rigid_transform_matrix;
    }

    // Accessors 
    inline double scan_duration() const { return m_scan_duration; }
    inline double number_of_lines() const { return m_number_of_lines; } 
    inline double samples_per_line() const { return m_samples_per_line; } 
    inline double sample_offset() const { return m_sample_offset; } 
    inline double focal_length() const { return m_focal_length; } 
    inline double along_scan_pixel_size() const { return m_along_scan_pixel_size; } 
    inline double across_scan_pixel_size() const { return m_across_scan_pixel_size; } 
    inline Quaternion<double> camera_pose() const { return m_camera_pose; } 
    inline Vector3 initial_position() const { return m_initial_position; } 
    inline Vector3 velocity_vector() const { return m_velocity_vector; } 
    
    inline void set_scan_duration(double val) { m_scan_duration = val; recompute_camera_matrix(); }
    inline void set_number_of_lines(double val) { m_number_of_lines = val; recompute_camera_matrix(); }
    inline void set_samples_per_line(double val) { m_samples_per_line = val; recompute_camera_matrix(); }
    inline void set_sample_offset(double val) { m_sample_offset = val; recompute_camera_matrix(); }
    inline void set_along_scan_pixel_size(double val) { m_along_scan_pixel_size = val; recompute_camera_matrix(); }
    inline void set_across_scan_pixel_size(double val) {m_across_scan_pixel_size = val; recompute_camera_matrix(); }
    inline void set_focal_length(double val) {m_focal_length = val; recompute_camera_matrix(); }
    inline void set_camera_pose(Quaternion<double> const& val) { m_camera_pose = val; recompute_camera_matrix(); }
    inline void set_initial_position(Vector3 const& val) { m_initial_position = val; recompute_camera_matrix(); }
    inline void set_velocity_vector(Vector3 const& val) { m_velocity_vector = val; recompute_camera_matrix(); }
    
  protected:
    //------------------------------------------------------------------
    // Protected Methods
    //------------------------------------------------------------------

    /// Take the parameters that are stored in the class, and use them
    /// to compute the camera model matrix.  This is the
    /// representation that is most useful when manipulating data
    /// using this camera model.
    /// 
    /// This camera model is based on Gupta and Hartley's paper entitled
    /// "Linear Pushbroom Cameras."
    void recompute_camera_matrix();

    Quaternion<double> m_camera_pose;
    Vector3 m_initial_position;
    Vector3 m_velocity_vector;
    
    Matrix<double> m_camera_matrix;
    Matrix<double> m_projection_matrix;
    Matrix<double> m_velocity_matrix;
    Matrix<double> m_rigid_transform_matrix;
    
    double m_scan_duration;
    double m_number_of_lines;
    double m_samples_per_line;
    double m_focal_length;
    double m_across_scan_pixel_size;
    double m_along_scan_pixel_size;
    double m_sample_offset;
  };  

  /// Output stream method for printing a summary of the linear
  /// pushbroom camera model parameters.
  std::ostream& operator<<( std::ostream& os, LinearPushbroomModel const& camera_model);

}}	// namespace vw::camera

#endif	//_VW_CAMERA_LINEARPUSHBROOM_MODEL_H_

