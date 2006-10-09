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

#include <vw/Camera/LinearPushbroomModel.h>

namespace vw {
namespace camera {

  // ------------------------------------------------
  //           Constructors / Destructors            
  // ------------------------------------------------

  LinearPushbroomModel::LinearPushbroomModel( double scan_duration,
                                              double number_of_lines, 
                                              double samples_per_line, 
                                              double sample_offset, 
                                              double focal_length, 
                                              double along_scan_pixel_size,
                                              double across_scan_pixel_size,
                                              Quaternion<double> const& camera_pose,
                                              Vector3 const& initial_position,
                                              Vector3 const& velocity_vector) {
    m_scan_duration = scan_duration;
    m_number_of_lines = number_of_lines;
    m_samples_per_line = samples_per_line;
    m_sample_offset = sample_offset;
    m_focal_length = focal_length;
    m_along_scan_pixel_size = along_scan_pixel_size;
    m_across_scan_pixel_size = across_scan_pixel_size;
    m_camera_pose = camera_pose;
    m_initial_position = initial_position;
    m_velocity_vector = velocity_vector;
    recompute_camera_matrix();
  }
  
  // ------------------------------------------------
  //           Interface
  // ------------------------------------------------
  Vector2 LinearPushbroomModel::vector_to_pixel(Vector3 const& vec) const {
    throw vw::NoImplErr() << "LinearPushbroomModel::vector_to_pixel is not yet implemented.";
  }
  
  Vector3 LinearPushbroomModel::pixel_to_vector(Vector2 const& pix) const {
    // u is the crosstrack (i.e., perspective) direction
    // v is the downtrack (i.e., orthographic) direction
    double u = pix[0], v = pix[1];

    // Compute the time for the given pixel
    double time_per_line = m_scan_duration / m_number_of_lines;
    double pixel_time = v * time_per_line;
    
    // The view_matrix takes vectors from the camera (extrinsic)
    // coordinate system to the world frame
    //
    // The position and veloctiy are not actually needed, since we are
    // purely interested in returning the direction of the ray at this
    // point and not its origin.
    Matrix<double,3,3> rot = transpose(m_camera_pose.rotation_matrix());
    Matrix<double,4,4> view_matrix;
    view_matrix.set_identity();
    submatrix(view_matrix, 0, 0, 3, 3) = rot;

    double pixel_size_u = m_across_scan_pixel_size;
    double center_offset = m_sample_offset * m_across_scan_pixel_size;
    
    // The viewplane is the y-z plane of the camera coordinate system.
    // Assuming the origin of the coordinate system is at the center
    // of projection, the image plane is z = +f, and the pixel
    // position in camera coordinates is:
    //
    // The following assumes u, v coords have origin at lower left hand
    // corner of image...
    double pixel_pos_y = (u * pixel_size_u) + center_offset;
    double f = m_focal_length;
    Vector<double, 4> pixel_pos(0.0, pixel_pos_y, f, 1.0);

    // Transform to world coordinates using the rigid rotation
    Vector<double, 4> direction_vec = view_matrix * pixel_pos;
    return normalize(Vector3 ( direction_vec[0],
                               direction_vec[1], 
                               direction_vec[2] ));
  }
  
  Vector3 LinearPushbroomModel::camera_center(Vector2 const& pix ) const {
    double t = pix[1] / m_number_of_lines * m_scan_duration;
    return m_initial_position + m_velocity_vector * t;
  }
  
  // ------------------------------------------------
  //                 Protected Methods                   
  // ------------------------------------------------]
  void LinearPushbroomModel::recompute_camera_matrix() {
  
   // Start by determining the rotation and translation (extrinsic
   // parameters).  This takes us from Mars Centered coordinates to the
   // coordinates of the spacecraft.
   //
   // We must also permute the entries of the quaternion, since our
   // quaternion to rotation matrix function expects the quaternion to
   // appear as [q1 q2 q3 w].  Note that we also include the R_moc
   // correction to move us from the spacecraft frame to the MOC frame.
    Matrix<double,3,3> R = m_camera_pose.rotation_matrix();

    /* For debugging
     *
     std::cout << "R is: \n " << "\n";
     vnl_matrix<double> T = (Q.rotation_matrix_transpose()).transpose();
     printf("\t  %1.15f \t %1.15f \t %1.15f \n", T(0,0), T(0,1), T(0,2));
     printf("\t  %1.15f \t %1.15f \t %1.15f \n", T(1,0), T(1,1), T(1,2));
     printf("\t  %1.15f \t %1.15f \t %1.15f \n", T(2,0), T(2,1), T(2,2));
    */

    Matrix<double,3,4> T1;
    submatrix(T1, 0, 0, 3, 3) = R;
    select_col(T1,3) = R * m_initial_position;
    
    /* For debugging 
     *
     std::cout << "T1 is: \n " << "\n";
     printf("\t  %1.15f \t %1.15f \t %1.15f \t %1.15f \n", T1(0,0), T1(0,1), T1(0,2), T1(0,3));
     printf("\t  %1.15f \t %1.15f \t %1.15f \t %1.15f \n", T1(1,0), T1(1,1), T1(1,2), T1(1,3));
     printf("\t  %1.15f \t %1.15f \t %1.15f \t %1.15f \n", T1(2,0), T1(2,1), T1(2,2), T1(2,3));
    */
    
    // Next, we compute the velocity transform.
    // 
    // We must first rotate the velocity vectors into the frame of reference of the MOC 
    // instrument.  Then we build the second transformation matrix, which takes us into 
    // coordinates of the form [t_im, y_im, z_im]. 
    Vector3 V = R * m_velocity_vector;
    Matrix<double,3,3> T2;
    T2.set_identity();
    T2(0,0) = 1/V(0);    T2(1,0) = -V(1)/V(0);   T2(2,0) = -V(2)/V(0);

    /* For debugging... 
     *
     std::cout << "T2 is: \n " << "\n";
     printf("\t  %1.15f \t %1.15f \t %1.15f \n", T2(0,0), T2(0,1), T2(0,2));
     printf("\t  %1.15f \t %1.15f \t %1.15f \n", T2(1,0), T2(1,1), T2(1,2));
     printf("\t  %1.15f \t %1.15f \t %1.15f \n", T2(2,0), T2(2,1), T2(2,2));
    */
    
    
    // Finally, we build the instrinsic transformation matrix from the ephemeris data.
    // This is the matrix that takes us from units of meters and seconds to units of pixels 
    // in image space.
    double line_integration_time  = m_scan_duration / m_number_of_lines;   // Line integration time
    
    Matrix<double> T3(3,3);
    T3.set_identity();
    T3(0,0) = 1 / line_integration_time;
    T3(1,1) = m_focal_length / m_across_scan_pixel_size;
    T3(1,2) = -1.0 * m_sample_offset;
    
    /* For debugging....
     *
     std::cout << "T3 is: \n " << "\n";
     printf("\t  %1.15f \t %1.15f \t %1.15f \n", T3(0,0), T3(0,1), T3(0,2));
     printf("\t  %1.15f \t %1.15f \t %1.15f \n", T3(1,0), T3(1,1), T3(1,2));
     printf("\t  %1.15f \t %1.15f \t %1.15f \n", T3(2,0), T3(2,1), T3(2,2));
    */
    
    // Finally, compose the complete camera model from T1, T2, and T3.
    m_camera_matrix = T3 * T2 * T1;
    m_projection_matrix = T3;
    m_velocity_matrix = T2;
    m_rigid_transform_matrix = T1;
    
    /*For debugging 
      std::cout << "CameraMatrix is: \n " << cameraMatrix << "\n";
    */
  }

  std::ostream& operator<<( std::ostream& os, LinearPushbroomModel const& camera_model) {
    os << "\n-------------------- Linear Pushbroom Camera Model -------------------\n\n";
    os << " Camera Pose      :   " << camera_model.camera_pose() << "\n";
    os << " Initial Position :   " << camera_model.initial_position() << "\n";
    os << " Velocity Vector  :   " << camera_model.velocity_vector() << "\n";
    os << " Scan Duration    :   " << camera_model.scan_duration() << "\n";
    os << " Number of Lines  :   " << camera_model.number_of_lines() << "\n";
    os << " Samples per Line :   " << camera_model.samples_per_line() << "\n";
    os << " Focal Length     :   " << camera_model.focal_length() << "\n";
    os << " Sample Offset    :   " << camera_model.sample_offset() << "\n";
    os << " Across Scan Pixel Size :   " << camera_model.across_scan_pixel_size() << "\n";
    os << " Along Scan Pixel Size  :   " << camera_model.along_scan_pixel_size() << "\n";
    os << "\n------------------------------------------------------------------------\n\n";
  }

}} // namespace vw::camera

