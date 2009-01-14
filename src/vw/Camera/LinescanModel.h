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

#include <vw/Math/Quaternion.h>
#include <vw/Camera/CameraModel.h>

namespace vw { 
namespace camera {

  template <class PositionFuncT, class PoseFuncT>
  class LinescanModel : public CameraModel { 

    // Extrinsics
    PositionFuncT m_position_func;
    PoseFuncT m_pose_func;
    
    // Intrinsics
    int m_number_of_lines;
    int m_samples_per_line;
    int m_sample_offset;
    double m_focal_length;
    double m_across_scan_pixel_size;
    double m_along_scan_pixel_size;
    std::vector<double> m_line_times;
    Vector3 m_pointing_vec;
    Vector3 m_u_vec;

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------

    /// This version of the constructor assumes that the line
    /// integration varies one each scanline of the image.  The
    /// line_integration_times vector must be equal to the
    /// number_of_lines
    LinescanModel( int number_of_lines, 
                   int samples_per_line, 
                   int sample_offset, 
                   double focal_length, 
                   double along_scan_pixel_size,
                   double across_scan_pixel_size,
                   std::vector<double> const& line_times,
                   Vector3 pointing_vec,
                   Vector3 u_vec,
                   PositionFuncT const& position_func,
                   PoseFuncT const& pose_func) : m_position_func(position_func),
                                                 m_pose_func(pose_func) {
      
      VW_ASSERT(int(line_times.size()) == number_of_lines,
                ArgumentErr() << "LinescanModel: number of line integration times does not match the number of scanlines.\n");
      
      // Intrinsics
      m_number_of_lines = number_of_lines;
      m_samples_per_line = samples_per_line;
      m_sample_offset = sample_offset;
      m_focal_length = focal_length;
      m_along_scan_pixel_size = along_scan_pixel_size;
      m_across_scan_pixel_size = across_scan_pixel_size;

      m_line_times = line_times;

      m_pointing_vec = normalize(pointing_vec);
      m_u_vec = normalize(u_vec);
    }

    /// This version of the constructor assumes that the line
    /// integration time is consistent throughout the image.
    LinescanModel( int number_of_lines, 
                   int samples_per_line, 
                   int sample_offset, 
                   double focal_length, 
                   double along_scan_pixel_size,
                   double across_scan_pixel_size,
                   double line_integration_time,
                   Vector3 pointing_vec,
                   Vector3 u_vec,
                   PositionFuncT const& position_func,
                   PoseFuncT const& pose_func) : m_position_func(position_func),
                                                 m_pose_func(pose_func) {
      
      // Intrinsics
      m_number_of_lines = number_of_lines;
      m_samples_per_line = samples_per_line;
      m_sample_offset = sample_offset;
      m_focal_length = focal_length;
      m_along_scan_pixel_size = along_scan_pixel_size;
      m_across_scan_pixel_size = across_scan_pixel_size;
      
      m_line_times.resize(number_of_lines);
      double sum = 0;
      for (int i = 0; i < number_of_lines; ++i) {
        m_line_times[i] = sum;
        sum += line_integration_time;
      }

      m_pointing_vec = normalize(pointing_vec);
      m_u_vec = normalize(u_vec);
    }
    
    virtual ~LinescanModel() {}
    virtual std::string type() const { return "Linescan"; }

    //------------------------------------------------------------------
    // Interface
    //------------------------------------------------------------------
    virtual Vector2 point_to_pixel(Vector3 const& point) const {
      vw_throw( vw::NoImplErr() << "LinescanModel::point_to_pixel is not yet implemented." );
      return Vector2(); // never reached
    }

    /// Given a pixel in image coordinates, what is the pointing
    /// vector in 3-space if you apply the camera model.
    ///
    /// Important Note: For linear pushbroom sensors, the orientation
    /// of the image is important.  In the Camera module we adopt the
    /// following convention for pushbroom imagers:
    ///
    /// u is the crosstrack (i.e., perspective) direction
    /// v is the downtrack (i.e., orthographic) direction
    virtual Vector3 pixel_to_vector(Vector2 const& pix) const {
      double u = pix[0], v = pix[1];
      
      // Check to make sure that this is a valid pixel
      if (int(round(v)) < 0 || int(round(v)) >= int(m_line_times.size()))
        vw_throw( PixelToRayErr() << "LinescanModel: requested pixel " << pix << " is not on a valid scanline." );

      // The view_matrix takes vectors from the camera (extrinsic)
      // coordinate system to the world frame
      //
      // The position and veloctiy are not actually needed, since we are
      // purely interested in returning the direction of the ray at this
      // point and not its origin.
      //
      // The v pixel need not be an integer in every case, therefore
      // we need to linearly interpolate line times that fall in
      // between pixels.
      int y = int(floor(pix[1]));
      double normy = pix[1] - y;
      double approx_line_time = double( m_line_times[y] + (m_line_times[y+1] - m_line_times[y]) * normy );

      Quaternion<double> pose = m_pose_func(approx_line_time);
      Matrix<double,3,3> rotation_matrix = transpose(pose.rotation_matrix());
      
      // The viewplane is the [pointing_vec cross u_vec] plane of the
      // camera coordinate system.  Assuming the origin of the
      // coordinate system is at the center of projection, the image
      // plane is at pointing_vec = +f, and the pixel position in
      // camera coordinates is:
      double pixel_pos_u = (u + m_sample_offset) * m_across_scan_pixel_size;
      Vector<double, 3> pixel_direction = pixel_pos_u * m_u_vec + m_focal_length * m_pointing_vec;
      
      // Transform to world coordinates using the rigid rotation
      return normalize(rotation_matrix * pixel_direction);
    }
    
    virtual Vector3 camera_center(Vector2 const& pix = Vector2() ) const {
      // Check to make sure that this is a valid pixel
      if (int(round(pix[1])) < 0 || int(round(pix[1])) >= int(m_line_times.size()))
        vw_throw( PixelToRayErr() << "LinescanModel: requested pixel " << pix << " is not on a valid scanline." );

      // The v pixel need not be an integer in every case, therefore
      // we need to linearly interpolate line times that fall in
      // between pixels.
      int y = int(floor(pix[1]));
      double normy = pix[1] - y;
      double approx_line_time = double( m_line_times[y] + (m_line_times[y+1] - m_line_times[y]) * normy );

      return m_position_func(approx_line_time);
    }

    /// Returns the pose (as a quaternion) of the camera for a given
    /// pixel.
    virtual Quaternion<double> camera_pose(Vector2 const& pix) const {
      // Check to make sure that this is a valid pixel
      if (int(round(pix[1])) < 0 || int(round(pix[1])) >= int(m_line_times.size()))
        vw_throw( PixelToRayErr() << "LinescanModel::camera_pose(): requested pixel " << pix << " is not on a valid scanline." );

      // The v pixel need not be an integer in every case, therefore
      // we need to linearly interpolate line times that fall in
      // between pixels.
      int y = int(floor(pix[1]));
      double normy = pix[1] - y;
      double approx_line_time = double( m_line_times[y] + (m_line_times[y+1] - m_line_times[y]) * normy );

      return m_pose_func(approx_line_time);
    }
    
    //------------------------------------------------------------------
    // Public Methods
    //------------------------------------------------------------------
    
    // Accessors 
    virtual double number_of_lines() const { return m_number_of_lines; } 
    virtual double samples_per_line() const { return m_samples_per_line; } 
    virtual double sample_offset() const { return m_sample_offset; } 
    virtual double focal_length() const { return m_focal_length; } 
    virtual double along_scan_pixel_size() const { return m_along_scan_pixel_size; } 
    virtual double across_scan_pixel_size() const { return m_across_scan_pixel_size; } 
    virtual std::vector<double> line_times() const { return m_line_times; }
    virtual Vector3 camera_position(double t) const { return m_position_func(t); }
    virtual Quaternion<double> camera_pose(double t) const { return m_pose_func(t); }
    
    virtual void set_number_of_lines(int val) { m_number_of_lines = val; }
    virtual void set_samples_per_line(int val) { m_samples_per_line = val; }
    virtual void set_sample_offset(int val) { m_sample_offset = val; }
    virtual void set_along_scan_pixel_size(double val) { m_along_scan_pixel_size = val; }
    virtual void set_across_scan_pixel_size(double val) {m_across_scan_pixel_size = val; }
    virtual void set_focal_length(double val) {m_focal_length = val; }
    virtual void set_line_times(std::vector<double> val) { m_line_times = val; }
  };  

  /// Output stream method for printing a summary of the linear
  /// pushbroom camera model parameters.
  template <class PositionFuncT, class PoseFuncT>
  std::ostream& operator<<( std::ostream& os, LinescanModel<PositionFuncT, PoseFuncT> const& camera_model) {
    os << "\n-------------------- Linescan Camera Model -------------------\n\n";
    os << " Camera center @ origin :   " << camera_model.camera_center(Vector2(0,0)) << "\n";
    os << " Number of Lines        :   " << camera_model.number_of_lines() << "\n";
    os << " Samples per Line       :   " << camera_model.samples_per_line() << "\n";
    os << " Sample Offset          :   " << camera_model.sample_offset() << "\n";
    os << " Focal Length           :   " << camera_model.focal_length() << "\n";
    os << " Across Scan Pixel Size :   " << camera_model.across_scan_pixel_size() << "\n";
    os << " Along Scan Pixel Size  :   " << camera_model.along_scan_pixel_size() << "\n";
    os << "\n------------------------------------------------------------------------\n\n";
    return os;
  }

}}	// namespace vw::camera

#endif	//_VW_CAMERA_LINEARPUSHBROOM_MODEL_H_

