// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
    /// pixel. It represents the rotation from the camera frame to world frame.
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
    Quat m_rotation;
    Quat m_rotation_inverse;

  public:
    AdjustedCameraModel(boost::shared_ptr<CameraModel> camera_model) : m_camera(camera_model) {
      m_rotation = Quat(math::identity_matrix<3>());
      m_rotation_inverse = Quat(math::identity_matrix<3>());
    }

    AdjustedCameraModel(boost::shared_ptr<CameraModel> camera_model,
                        Vector3 const& translation, Quat const& rotation) :
      m_camera(camera_model), m_translation(translation), m_rotation(rotation), m_rotation_inverse(inverse(rotation)) {}

    virtual ~AdjustedCameraModel() {}
    virtual std::string type() const { return "Adjusted"; }

    Vector3 translation() const { return m_translation; }
    Quat rotation() const { return m_rotation; }
    Matrix<double,3,3> rotation_matrix() const { return m_rotation.rotation_matrix(); }
    Vector3 axis_angle_rotation() const;
    void set_rotation(Quat const&);

    template <class MatrixT>
    void set_rotation(MatrixBase<MatrixT> const& m) {
      m_rotation = Quat(m.impl());
      m_rotation_inverse = inverse(m_rotation);
    }
    template <class VectorT>
    void set_translation(VectorBase<VectorT> const& v) {
      m_translation = v.impl();
    }
    template <class VectorT>
    void set_axis_angle_rotation(VectorBase<VectorT> const& v) {
      this->set_rotation( axis_angle_to_quaternion(v.impl()) );
    }

    virtual Vector2 point_to_pixel (Vector3 const&) const;
    virtual Vector3 pixel_to_vector (Vector2 const&) const;
    virtual Vector3 camera_center (Vector2 const&) const;
    virtual Quat camera_pose(Vector2 const&) const;

    void write(std::string const&);
    void read(std::string const&);

    friend std::ostream& operator<<(std::ostream&, AdjustedCameraModel const&);
  };

  std::ostream& operator<<(std::ostream&, AdjustedCameraModel const&);

  /// Error during projection of a 3D point onto the image plane.
  VW_DEFINE_EXCEPTION(PointToPixelErr, vw::Exception);

  /// Error during reverse projection of a pixel to a pointing vector
  /// from the camera center.
  VW_DEFINE_EXCEPTION(PixelToRayErr, vw::Exception);

}} // namespace vw::camera

#endif // __VW_CAMERA_CAMERAMODEL_H__
