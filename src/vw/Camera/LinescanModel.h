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

/// \file LinescanModel.h
/// 
/// A generic linescan camera model object
///
/// 
#ifndef _VW_CAMERA_LINESCAN_MODEL_H_
#define _VW_CAMERA_LINESCAN_MODEL_H_

#include <vw/Camera/CameraModel.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Quaternion.h>

namespace vw { 
namespace camera {
  
  class LinescanModel : public CameraModel { 
  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    LinearPushbroomModel() {}
    LinearPushbroomModel( double number_of_lines, 
                          double samples_per_line, 
                          double sample_offset, 
                          double focal_length, 
                          double along_scan_pixel_size,
                          double across_scan_pixel_size,
                          std::vector<double> const& line_integration_times,
                          std::vector<Quaternion<double> > const& camera_poses,
                          std::vector<Vector3> const& positions);

    virtual ~LinearPushbroomModel() {}
    
    //------------------------------------------------------------------
    // Interface
    //------------------------------------------------------------------
    virtual Vector2 point_to_pixel(Vector3 const& vec) const;

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
    
    // Accessors 
    inline double number_of_lines() const { return m_number_of_lines; } 
    inline double samples_per_line() const { return m_samples_per_line; } 
    inline double sample_offset() const { return m_sample_offset; } 
    inline double focal_length() const { return m_focal_length; } 
    inline double along_scan_pixel_size() const { return m_along_scan_pixel_size; } 
    inline double across_scan_pixel_size() const { return m_across_scan_pixel_size; } 
    inline std::vector<Quaternion<double> > camera_poses() const { return m_camera_poses; } 
    inline std::vector<Vector3> positions() const { return m_positions; } 
    inline std::vector<float> line_integration_times() const { return m_line_integration_times; }
    
    inline void set_number_of_lines(double val) { m_number_of_lines = val; }
    inline void set_samples_per_line(double val) { m_samples_per_line = val; }
    inline void set_sample_offset(double val) { m_sample_offset = val; }
    inline void set_along_scan_pixel_size(double val) { m_along_scan_pixel_size = val; }
    inline void set_across_scan_pixel_size(double val) {m_across_scan_pixel_size = val; }
    inline void set_focal_length(double val) {m_focal_length = val; }
    inline void set_camera_poses(std::vector<Quaternion<double> > const& val) { m_camera_pose = val;}
    inline void set_positions(std::vector<Vector3> const& val) { m_positions = val;}
    inline void set_line_integration_times(std::vector<float> val) const { m_line_integration_times = val; }

  protected:
    //------------------------------------------------------------------
    // Protected Methods
    //------------------------------------------------------------------

    std::vector<Quaternion<double> > m_camera_poses;
    std::vector<Vector3> m_positions;
    
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

