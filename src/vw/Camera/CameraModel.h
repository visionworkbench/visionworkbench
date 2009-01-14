// __BEGN_LICENSE__
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

/// \file CameraModel.h
/// 
/// This file contains the abstract base class from which all camera
/// models must derive their interface.
/// 
#ifndef __VW_CAMERA_CAMERAMODEL_H__
#define __VW_CAMERA_CAMERAMODEL_H__

#include <fstream>
#include <vw/Math/Quaternion.h>

namespace vw { 
namespace camera {

  /// This is the abstract base class for a camera model object.  You
  /// should create a child that adheres to the interface set forth
  /// below.  If your camera model subclass does not implement one of
  /// these methods, it should vw_throw a vw::NotImplErr() exception.
  ///
  /// The most fundamental camera model operation is the forward
  /// projection method ( point_to_pixel() ) that "images" a point in 3D
  /// space; this returns the pixel location where the projection
  /// intersects the image plane.
  ///
  /// It is often necessary to perform the inverse operation:
  /// projecting a ray from the camera center through the pixel in an
  /// image plane to determine the direction of the original 3D point.
  /// The 3D point corresponding to pixel (u,v) lies along a ray whose
  /// origin can be determined by calling camera_center(Vector2(u,v))
  /// and whose direction can be determined by calling
  /// pixel_to_vector(Vector2(u,v)).  Note that the absolute position
  /// of the 3D point cannot be computed from the pixel location in a
  /// single image; this information is lost when the point is forward
  /// projected.  However stereo vision techniques can be used to
  /// determine the location of the original point by intersecting two
  /// rays from two distinct cameras.
  class CameraModel {

  public:

    virtual ~CameraModel() {}

    //------------------------------------------------------------------
    // Generic Camera Model Interface
    //------------------------------------------------------------------

    /// Computes the image of the point 'point' in 3D space on the
    /// image plane.  Returns a pixel location (col, row) where the
    /// point appears in the image.  It is possible that the selected
    /// point will not be imaged by the camera (e.g. if it lies behind
    /// the camera).  In this case the method should vw_throw a
    /// vw::camera::PointToPixelErr()
    virtual Vector2 point_to_pixel (Vector3 const& point) const = 0;

    /// Returns a pointing vector from the camera center through the
    /// position of the pixel 'pix' on the image plane.  For
    /// consistency, the pointing vector should generally be
    /// normalized.
    virtual Vector3 pixel_to_vector(Vector2 const& pix) const = 0;

    /// Returns the camera center in the frame of reference of the
    /// camera model.  This method is often used to find the origin of
    /// a ray emanating from the focal point of the camera through a
    /// pixel on the image plane (e.g. for computing ray-ray
    /// intersection in a stereo vision algorithm).
    virtual Vector3 camera_center(Vector2 const& pix) const = 0;

    /// Subclasses must define a method that return the camera type as a string.
    virtual std::string type() const = 0;

    /// Returns the pose (as a quaternion) of the camera for a given
    /// pixel.
    virtual Quaternion<double> camera_pose(Vector2 const& /*pix*/) const {
      vw_throw( NoImplErr() << "CameraModel: this camera model has not implemented camera_pose()" );
      return Quaternion<double>();
    }

  };  


  /// This class is useful if you have an existing camera model, and
  /// you want to systematically "tweak" its extrinsic parameters
  /// (position and pose).  This is particularly useful in Bundle
  /// Adjustment.
  class AdjustedCameraModel : public CameraModel {
    
    boost::shared_ptr<CameraModel> m_camera;
    Vector3 m_translation;
    Quaternion<double> m_rotation;
    Quaternion<double> m_rotation_inverse;
    
  public:
    AdjustedCameraModel(boost::shared_ptr<CameraModel> camera_model) : m_camera(camera_model) {
      m_rotation = math::Quaternion<double>(math::identity_matrix<3>());
      m_rotation_inverse = math::Quaternion<double>(math::identity_matrix<3>());
    }

    AdjustedCameraModel(boost::shared_ptr<CameraModel> camera_model,
                        Vector3 const& translation, math::Quaternion<double> const& rotation) : 
      m_camera(camera_model), m_translation(translation), m_rotation(rotation), m_rotation_inverse(inverse(rotation)) {}

    virtual ~AdjustedCameraModel() {}
    virtual std::string type() const { return "Adjusted"; }

    Vector3 translation() const { return m_translation; }
    Quaternion<double> rotation() const { return m_rotation; }
    Matrix<double,3,3> rotation_matrix() const { return m_rotation.rotation_matrix(); }
    Vector3 axis_angle_rotation() const {
      Quaternion<double> quat = this->rotation();
      Vector4 rot;
      rot(0) = quat.w();
      rot(1) = quat.x();
      rot(2) = quat.y();
      rot(3) = quat.z();
      if (rot(0) > 1)
        rot = normalize(rot);
      double angle = 2 * acos(rot[0]);
      double s = sqrt(1-rot[0]*rot[0]); // assuming quaternion normalised then w is less than 1, so term always positive.
      Vector3 result;
      // if s is close to zero, then direction of axis is not important
      if (s < 0.001) {
        result[0] = rot.x();   // if it is important that axis is normalised then replace with x=1; y=z=0
        result[1] = rot.y();
        result[2] = rot.z();
      } else {
        result[0] = rot.x()/s; // normalize axis
        result[1] = rot.y()/s;
        result[2] = rot.z()/s;
      }

      result *= angle;
      return result;
    }

    void set_translation(Vector3 const& translation) { m_translation = translation; }
    void set_rotation(Quaternion<double> const& rotation) { 
      m_rotation = rotation;
      m_rotation_inverse = inverse(m_rotation);
    }
    void set_rotation(Matrix<double,3,3> const& rotation) {
      m_rotation = Quaternion<double>(rotation);
      m_rotation_inverse = inverse(m_rotation);
    }

    void set_axis_angle_rotation(Vector3 const& axis_angle) {
      Quaternion<double> rot; 
      double angle = norm_2(axis_angle);
      Vector3 temp = normalize(axis_angle);
      double s = sin(angle/2);
      if (angle == 0) {
        rot.x() = 0;
        rot.y() = 0;
        rot.z() = 0;
      } else {
        rot.x() = temp[0] * s;
        rot.y() = temp[1] * s;
        rot.z() = temp[2] * s;
      }
      rot.w() = cos(angle/2);
      this->set_rotation(rot);
    }

    virtual Vector2 point_to_pixel (Vector3 const& point) const {
      Vector3 offset_pt = point-m_camera->camera_center(Vector2(0,0))-m_translation;
      Vector3 new_pt = m_rotation_inverse.rotate(offset_pt) + m_camera->camera_center(Vector2(0,0));
      return m_camera->point_to_pixel(new_pt);
    }

    virtual Vector3 pixel_to_vector (Vector2 const& pix) const {
      return m_rotation.rotate(m_camera->pixel_to_vector(pix));
    }

    virtual Vector3 camera_center (Vector2 const& pix) const {
      return m_camera->camera_center(pix) + m_translation;
    }

    virtual Quaternion<double> camera_pose(Vector2 const& pix) const {
      return m_camera->camera_pose(pix)*m_rotation_inverse;      
    }

    void write(std::string filename) {
      std::ofstream ostr(filename.c_str());
      ostr << m_translation[0] << " " << m_translation[1] << " " << m_translation[2] << "\n";
      ostr << m_rotation.w() << " " << m_rotation.x() << " " << m_rotation.y() << " " << m_rotation.z() << "\n";
    }

    void read(std::string filename) {
      Quaternion<double> rot;
      Vector3 pos;
      std::ifstream istr(filename.c_str());
      istr >> pos[0] >> pos[1] >> pos[2];
      istr >> rot.w() >> rot.x() >> rot.y() >> rot.z();
      this->set_translation(pos);
      this->set_rotation(rot);
    }

  };


  /// Error during projection of a 3D point onto the image plane.
  VW_DEFINE_EXCEPTION(PointToPixelErr, vw::Exception);

  /// Error during reverse projection of a pixel to a pointing vector
  /// from the camera center.
  VW_DEFINE_EXCEPTION(PixelToRayErr, vw::Exception);

}}	// namespace vw::camera

#endif // __VW_CAMERA_CAMERAMODEL_H__
