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

/// \file PinholeModel.h
/// 
/// This file contains the pinhole camera model.
///
#ifndef __VW_CAMERAMODEL_PINHOLE_H__
#define __VW_CAMERAMODEL_PINHOLE_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

#include <vw/Camera/CameraModel.h>
#include <vw/Camera/LensDistortion.h>

#include <iostream>

namespace vw { 
namespace camera {

  /// This is a simple "generic" pinhole camera model.
  ///
  /// The extrinsic parameters are the camera center (translation) and
  /// rotation matrix (pose) of the camera. This transforms from world
  /// frame to the frame of reference of the camera.  In the frame of
  /// the camera, the camera's pointing vector is the +z unit vector,
  /// and the image plane is aligned such that the positive x-pixel
  /// direction (increasing image columns) is the camera frame's +x
  /// vector, and the positive y-pixel direction (increasing image
  /// rows) is the frame's -y vector.  Note that this discrepancy in y
  /// frames is due to the fact that images stored in memory are most
  /// naturally indexed starting in the upper left hand corner.  The
  /// intrinsic portion of the camera matrix is nominally stored as
  ///
  ///    [  fx   0   cx  ]
  /// K= [  0   -fy  cy  ]
  ///    [  0    0   1   ]
  ///
  /// with fx, fy the focal length of the system (in horizontal and
  /// vertical pixels), and (cx, cy) the pixel coordinates of the
  /// central pixel (the principal point on the image
  /// plane). Combining both the intrinsic camera matrix K with the
  /// extrinsic matrices, R and C, we see that a real-world point (x,
  /// y, z), to pixel p in an image by:
  ///
  ///     [  row  ]                     [ x ]
  /// p = [  col  ]  =  K * [R | -R C]  [ y ]
  ///     [   w   ]                     [ z ]
  ///
  /// p is then in homogenous coordinates, so the w has to be divided
  /// out so that w=1. Here R and C are the extrinsic parameters
  /// (rotation and position of the camera center in world
  /// coordinates, respectively).
  /// 
  /// XXX: Add some comments here about the distortion model...
  ///
  class PinholeModel : public CameraModel {
    boost::shared_ptr<LensDistortion> m_distortion_model_ptr;
    Matrix<double,3,4> m_camera_matrix;

    // Stored for easy access.
    Matrix<double,3,3> m_intrinsics;
    Matrix<double,3,3> m_rotation;
    Vector3 m_camera_center;
    
    // Cached values for pixel_to_vector
    Vector3 m_camera_matrix_rightmost_row;
    Matrix<double,3,3> m_inv_camera_matrix_block;

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    
    /// Initialize an empty camera model.
    PinholeModel() {
      m_intrinsics.set_identity();
      fill(m_camera_center,0.0);
      m_rotation.set_identity();
      this->rebuild_camera_matrix();

      m_distortion_model_ptr = boost::shared_ptr<LensDistortion>(new NullLensDistortion());
    }
    
    /// Initialize the pinhole model with explicit parameters.
    /// 
    /// This limits you to cameras where the intrinsic matrix is of
    /// the form:
    ///
    ///     | fx        cx |
    /// K = |     -fy   cy |
    ///     |           1  |
    ///
    /// Here, the focal length f is given in units of pixels, as are
    /// the principal point offset values, px and py.  If you are
    /// starting from a focal length in a physical unit (e.g. meters),
    /// you can find the focal length in pixels by dividing by the
    /// pixel scale (usually in meters/pixel).  Note the negative sign
    /// in the (2,2) entry.  This reflects the difference between the
    /// camera coordinate system (+y increases as you move up the
    /// image) and the VW standard frame of reference wherein (0,0) is
    /// the upper left hand corner of the image and the coordinates
    /// increase as you move down the image.
    ///
    PinholeModel(Vector3 camera_center, 
                 Matrix<double,3,3> rotation,
                 double fx, double fy, 
                 double cx, double cy,
                 LensDistortion const& distortion_model) : m_camera_center(camera_center),
                                                           m_rotation(rotation) {
      m_intrinsics.set_identity();
      m_intrinsics(0,0) = fx;
      m_intrinsics(1,1) = -fy;
      m_intrinsics(0,2) = cx;
      m_intrinsics(1,2) = cy;
      rebuild_camera_matrix();

      m_distortion_model_ptr = distortion_model.copy();
      m_distortion_model_ptr->set_parent_camera_model(this);
    }

    /// Construct a basic pinhole model with no lens distortion
    PinholeModel(Vector3 camera_center, 
                 Matrix<double,3,3> rotation,
                 double fx, double fy, 
                 double cx, double cy) : m_camera_center(camera_center),
                                         m_rotation(rotation) {
      m_intrinsics.set_identity();
      m_intrinsics(0,0) = fx;
      m_intrinsics(1,1) = -fy;
      m_intrinsics(0,2) = cx;
      m_intrinsics(1,2) = cy;
      rebuild_camera_matrix();
      m_distortion_model_ptr = boost::shared_ptr<LensDistortion>(new NullLensDistortion());
    }


    /// Initialize the pinhole model.  The user can specify an
    /// arbitrary intrinsic camera parameter matrix.
    PinholeModel(Vector3 camera_center, 
                 Matrix<double,3,3> rotation,
                 Matrix<double,3,3> intrinsics,
                 LensDistortion const& distortion_model) : m_camera_center(camera_center),
                                                           m_rotation(rotation),
                                                           m_intrinsics(intrinsics) {
      rebuild_camera_matrix();

      m_distortion_model_ptr = distortion_model.copy();
      m_distortion_model_ptr->set_parent_camera_model(this);
    }
    
    /// Construct a basic pinhole model with no lens distortion
    PinholeModel(Vector3 camera_center, 
                 Matrix<double,3,3> rotation,
                 Matrix<double,3,3> intrinsics) : m_camera_center(camera_center),
                                                  m_rotation(rotation),
                                                  m_intrinsics(intrinsics) {
      rebuild_camera_matrix();
      m_distortion_model_ptr = boost::shared_ptr<LensDistortion>(new NullLensDistortion());
    }

    virtual ~PinholeModel() {}
    
    //------------------------------------------------------------------
    // Methods
    //------------------------------------------------------------------
    virtual Vector2 point_to_pixel(Vector3 const& point) const {

      // Multiply the pixel location by the camera matrix.
      double denominator = m_camera_matrix(2,0)*point(0) + m_camera_matrix(2,1)*point(1) + m_camera_matrix(2,2)*point(2) + m_camera_matrix(2,3);
      Vector2 pixel = Vector2( (m_camera_matrix(0,0)*point(0) + m_camera_matrix(0,1)*point(1) + m_camera_matrix(0,2)*point(2) + m_camera_matrix(0,3)) / denominator,
                               (m_camera_matrix(1,0)*point(0) + m_camera_matrix(1,1)*point(1) + m_camera_matrix(1,2)*point(2) + m_camera_matrix(1,3)) / denominator);

      // Apply the lens distortion model
      return m_distortion_model_ptr->get_distorted_coordinates(pixel);
    }

    virtual Vector3 pixel_to_vector (Vector2 const& pix) const {

      // Apply the inverse lens distortion model
      Vector2 undistorted_pix = m_distortion_model_ptr->get_undistorted_coordinates(pix);
      
      // Compute the direction of the ray emanating from the camera center.
      Vector3 p(0,0,1);
      subvector(p,0,2) = undistorted_pix;
      return normalize(m_inv_camera_matrix_block * (p-m_camera_matrix_rightmost_row));
    }

    virtual Vector3 camera_center(Vector2 const& pix = Vector2() ) const { return m_camera_center; };
    void set_camera_center(Vector3 const& position) { m_camera_center = position; rebuild_camera_matrix(); }
    
    Matrix<double,3,3> camera_pose() const { return m_rotation; };
    void set_camera_pose(Matrix<double,3,3> const& pose) { m_rotation = pose; rebuild_camera_matrix(); }

    Matrix<double,3,3> intrinsic_matrix() const { return m_intrinsics; };
    void set_intrinsic_matrix(Matrix<double,3,3> const& intrinsics) { m_intrinsics = intrinsics; rebuild_camera_matrix(); }

    boost::shared_ptr<LensDistortion> lens_distortion() const { return m_distortion_model_ptr; };
    void set_lens_distortion(LensDistortion const& distortion) { m_distortion_model_ptr = distortion.copy(); }
    
  private:
    void rebuild_camera_matrix() {
      Matrix<double,3,4> extrinsics;
      submatrix(extrinsics,0,0,3,3) = m_rotation;
      select_col(extrinsics,3) = -m_rotation * m_camera_center;
      m_camera_matrix = m_intrinsics * extrinsics;
      
      m_inv_camera_matrix_block = inverse(submatrix(m_camera_matrix,0,0,3,3));
      m_camera_matrix_rightmost_row = select_col(m_camera_matrix,3);
    }

  };

  /// TSAI Lens Distortion Model
  /// 
  /// For a given set of observed (distorted) pixel coordinates, return the 
  /// location where the pixel would have appeared if there were no lens distortion.
  /// 
  /// The equations which produce these are:
  /// (u, v) = undistorted coordinates
  /// (u', v') = observed (distorted) coordinates
  /// (x, y) = object coordinates of projected point
  /// r2 = x * x + y * y   -- square of distance from object to primary vector
  /// k1, k2 are radial distortion parameters; p1, p2 are tangential distortion
  /// parameters. principal point is at (cx, cy).
  ///
  /// u' = u + (u - cx) * (k1 * r2 + k2 * r4 + 2 * p1 * y + p2 * (r2/x + 2x))
  /// v' = v + (v - cy) * (k1 * r2 + k2 * r4 + 2 * p2 * x + p1 * (r2/y + 2y))
  ///
  /// k1 is distortion[0], k2 is distortion[1],  p1 is distortion[2], p2 is distortion[3]
  class TsaiLensDistortion : public LensDistortionBase<TsaiLensDistortion, PinholeModel> {
    Vector4 m_distortion;
  public:
    TsaiLensDistortion(Vector4 params) : m_distortion(params) {}
    virtual ~TsaiLensDistortion() {}

    Vector4 distortion_parameters() { return m_distortion; }

    virtual Vector2 get_distorted_coordinates(Vector2 const& p) const {
      double du = p[0] - this->camera_model().intrinsic_matrix()(0, 2);
      double dv = p[1] - this->camera_model().intrinsic_matrix()(1, 2);
      
      double x = du / this->camera_model().intrinsic_matrix()(0, 0);  // find (x, y) using similar triangles;
      double y = dv / this->camera_model().intrinsic_matrix()(1, 1);  // assumed z=1.
      
      double x1 = m_distortion[3] / x;
      double y1 = m_distortion[2] / y;
      
      double r2 = x * x + y * y;
      
      double x3 = 2.0 * m_distortion[3] * x;
      double y3 = 2.0 * m_distortion[2] * y;
      
      double bx = r2 * (m_distortion[0] + r2 * m_distortion[1]) + x3 + y3;
      double by = bx + r2 * y1;
      bx += r2 * x1;
      
      return Vector2( p[0] + bx * du,
                      p[1] + by * dv);
    }
  };

//   /// Given two pinhole camera models, this method returns two new camera
//   /// models that have been epipolar rectified.
//   template <>
//   void epipolar(PinholeModel<NoLensDistortion> const& src_camera0, 
//                 PinholeModel<NoLensDistortion> const& src_camera1, 
//                 PinholeModel<NoLensDistortion> &dst_camera0, 
//                 PinholeModel<NoLensDistortion> &dst_camera1);

  /// Function to remove lens distortion from a pinhole camera model.
  inline PinholeModel linearize_camera(PinholeModel const& camera_model) {
    return PinholeModel(camera_model.camera_center(),
                        camera_model.camera_pose(),
                        camera_model.intrinsic_matrix());
  }

  std::ostream& operator<<(std::ostream& str, PinholeModel const& model) {
    str << "Pinhole camera: \n";
    str << "\tCamera Center: " << model.camera_center() << "\n";
    str << "\tRotation Matrix: " << model.camera_pose() << "\n";
    str << "\tIntrinsic Matrix: " << model.intrinsic_matrix() << "\n";
  }

}}	// namespace vw::camera

#endif	//__CAMERAMODEL_CAHV_H__
