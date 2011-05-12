// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file PinholeModel.h
///
/// This file contains the pinhole camera model.
///
#ifndef __VW_CAMERAMODEL_PINHOLE_H__
#define __VW_CAMERAMODEL_PINHOLE_H__

#include <vw/Math/Quaternion.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Camera/LensDistortion.h>

#include <iostream>
#include <fstream>

namespace vw {
namespace camera {

  /// This is a simple "generic" pinhole camera model.
  ///
  /// To specify the EXTRINSIC paramters of the camera, we specify the
  /// position of the camera center in the world frame (m_camera_center)
  /// and the pose (or orientation) of the camera in world frame (m_rotation)
  /// (which is the transformation from the camera's frame to the world frame).
  /// In the default Vision Workbench camera frame, the camera's pointing vector
  /// is the +z unit vector, and the image plane is aligned such that the
  /// positive x-pixel direction (increasing image columns) is the camera frame's
  /// +x vector, and the positive y-pixel direction (increasing image
  /// rows) is the frame's -y vector.  Note that this discrepancy in y
  /// frames is due to the fact that images stored in memory are most
  /// naturally indexed starting in the upper left hand corner.
  ///
  /// --->The user can re-define the direction of increasing x-pixel,
  ///     increasing y-pixel, and pointing vector by specifying
  ///     orthonormal vectors u,v,w. These are intended to simplify
  ///     movement between different camera coordinate conventions,
  ///     rather than encoding the complete rotation between world
  ///     and camera coordinate frames.
  ///
  /// The INTRINSIC portion of the camera matrix is nominally stored as
  ///
  ///    [  fx   0   cx  ]
  /// K= [  0   -fy  cy  ]
  ///    [  0    0   1   ]
  ///
  /// with fx, fy the focal length of the system (in horizontal and
  /// vertical pixels), and (cx, cy) the pixel offset of the
  /// principal point of the camera on the image plane. --Note that
  /// the default v direction is <0,-1,0>, so
  /// K will be create with a POSITIVE fy term in the center; it
  /// becomes negative when multiplied with the v_direction vector).
  ///
  /// Combining both the intrinsic camera matrix K with the
  /// extrinsic matrices, (u,v,w rotation, R and C) we see that a real-world point (x,
  /// y, z), to pixel p in an image by:
  ///
  ///     [  row  ]         [ -u- ]               [ x ]
  /// p = [  col  ]  =  K * [ -v- ] * [R | -R C]  [ y ]
  ///     [   w   ]         [ -w- ]               [ z ]
  ///
  /// p is then in homogenous coordinates, so the w has to be divided
  /// out so that w=1. Here R and C are the extrinsic parameters; R and -R*C
  /// rotate and translate a vector in world coordinates into camera coordinates.
  ///
  ///
  ///  The Tsai distortion model describes radial and tangential lens distortion. See below.
  ///

  class PinholeModel : public CameraModel {
    typedef boost::shared_ptr<const LensDistortion> DistortPtr;

    DistortPtr m_distortion;
    Matrix<double,3,4> m_camera_matrix;

    // Stored for easy access.
    Vector3 m_camera_center;
    Matrix<double,3,3> m_rotation;
    Matrix<double,3,3> m_intrinsics;
    Matrix<double,3,4> m_extrinsics;

    // Intrinsic parameters, in pixel units
    double m_fu, m_fv, m_cu, m_cv;

    // Vectors that define how the coordinate system of the camera
    // relate to the directions: +u (increasing image columns), +v
    // (increasing image rows), and +w (out along optical axis)
    Vector3 m_u_direction;
    Vector3 m_v_direction;
    Vector3 m_w_direction;

    // Pixel Pitch, if the above units were not in pixels this should
    // convert it to that. For example, if distortion and focal length
    // have been described in mm. Pixel pitch would then be described
    // in mm/px.
    double m_pixel_pitch;

    // Cached values for pixel_to_vector
    Matrix<double,3,3> m_inv_camera_transform;

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------

    /// Initialize an empty camera model.
    PinholeModel() : m_distortion(DistortPtr(new NullLensDistortion)),
                     m_camera_center(Vector3(0,0,0)),
                     m_fu(1), m_fv(1), m_cu(0), m_cv(0),
                     m_u_direction(Vector3(1,0,0)),
                     m_v_direction(Vector3(0,-1,0)),
                     m_w_direction(Vector3(0,0,1)), m_pixel_pitch(1) {

      m_rotation.set_identity();
      this->rebuild_camera_matrix();
    }

    /// Initialize from a file on disk.
    PinholeModel(std::string const& filename) : m_distortion(DistortPtr(new NullLensDistortion)) {
      read(filename);
    }

    /// Initialize the pinhole model with explicit parameters.
    ///
    /// The user supplies the basic intrinsic camera parameters:
    ///
    /// f_u : focal length (in units of pixels) in the u direction
    /// f_v : focal length (in units of pixels) in the v direction
    /// c_u : principal point offset (in pixels) in the u direction
    /// c_v : principal point offset (in pixels) in the v direction
    ///
    /// The direction vectors define how the coordinate system of the
    /// camera relate to the directions: +u (increasing image
    /// columns), +v (increasing image rows), and +w (complete the RH
    /// coordinate system with u and v -- points into the image)
    ///
    /// If you start from a focal length in a physical unit
    /// (e.g. meters), you can find the focal length in pixels by
    /// dividing by the pixel scale (usually in meters/pixel).
    ///
    /// Remember that the VW standard frame of reference is such that
    /// (0,0) is the upper left hand corner of the image and the v
    /// coordinates increase as you move down the image. There is an
    /// illustration in the VisionWorkbook.
    ///

    PinholeModel(Vector3 camera_center, Matrix<double,3,3> rotation,
                 double f_u, double f_v, double c_u, double c_v,
                 Vector3 u_direction, Vector3 v_direction,
                 Vector3 w_direction,
                 LensDistortion const& distortion_model) : m_distortion(DistortPtr(distortion_model.copy())),
                                                           m_camera_center(camera_center),
                                                           m_rotation(rotation),
                                                           m_fu(f_u), m_fv(f_v), m_cu(c_u), m_cv(c_v),
                                                           m_u_direction(u_direction),
                                                           m_v_direction(v_direction),
                                                           m_w_direction(w_direction),
                                                           m_pixel_pitch(1) {
      this->rebuild_camera_matrix();
    }

    /// Initialize the pinhole model with explicit parameters.
    ///
    /// The user supplies the basic intrinsic camera parameters:
    ///
    /// f_u : focal length (in units of pixels) in the u direction
    /// f_v : focal length (in units of pixels) in the v direction
    /// c_u : principal point offset (in pixels) in the u direction
    /// c_v : principal point offset (in pixels) in the v direction
    ///
    /// The direction vectors defining the coordinate system of the
    /// camera are set to default values in this version of the
    /// constructor:
    ///
    ///   +u (increasing image columns)                     =  +X   [1 0 0]
    ///   +v (increasing image rows)                        =  +Y   [0 -1 0]
    ///   +w (complete the RH coordinate system with
    ///       u and v -- points into the image)             =  +Z   [0 0 1]
    ///
    /// If you start from a focal length in a physical unit
    /// (e.g. meters), you can find the focal length in pixels by
    /// dividing by the pixel scale (usually in meters/pixel).
    ///
    /// Remember that the VW standard frame of reference is such that
    /// (0,0) is the upper left hand corner of the image and the v
    /// coordinates increase as you move down the image.
    ///
    PinholeModel(Vector3 camera_center, Matrix<double,3,3> rotation,
                 double f_u, double f_v, double c_u, double c_v,
                 LensDistortion const& distortion_model) : m_distortion(DistortPtr(distortion_model.copy())),
                                                           m_camera_center(camera_center),
                                                           m_rotation(rotation),
                                                           m_fu(f_u), m_fv(f_v), m_cu(c_u), m_cv(c_v),
                                                           m_u_direction(Vector3(1,0,0)),
                                                           m_v_direction(Vector3(0,-1,0)),
                                                           m_w_direction(Vector3(0,0,1)),
                                                           m_pixel_pitch(1) {
      rebuild_camera_matrix();
    }


    /// Construct a basic pinhole model with no lens distortion
    PinholeModel(Vector3 camera_center, Matrix<double,3,3> rotation,
                 double f_u, double f_v,
                 double c_u, double c_v) : m_distortion(DistortPtr(new NullLensDistortion)),
                                           m_camera_center(camera_center),
                                           m_rotation(rotation),
                                           m_fu(f_u), m_fv(f_v), m_cu(c_u), m_cv(c_v),
                                           m_u_direction(Vector3(1,0,0)),
                                           m_v_direction(Vector3(0,-1,0)),
                                           m_w_direction(Vector3(0,0,1)),
                                           m_pixel_pitch(1) {
      rebuild_camera_matrix();
    }

    virtual std::string type() const { return "Pinhole"; }
    virtual ~PinholeModel() {}

    /// Read / Write a pinhole model from a file on disk.
    /// Files will end in format .pinhole
    void read(std::string const& filename);
    void write(std::string const& filename) const;

    /// DEPRECATED FILE IO
    void read_file(std::string const& filename) VW_DEPRECATED;
    void write_file(std::string const& filename) const VW_DEPRECATED;
    void read_old_file(std::string const& filename) VW_DEPRECATED;

    //------------------------------------------------------------------
    // Methods
    //------------------------------------------------------------------

    //  Computes the image of the point 'point' in 3D space on the
    //  image plane.  Returns a pixel location (col, row) where the
    //  point appears in the image.
    virtual Vector2 point_to_pixel(Vector3 const& point) const;

    // Is a valid projection of point is possible?
    // This is equal to: Is the point in front of the camera (z > 0)
    // after extinsic transformation?
    virtual bool projection_valid(Vector3 const& point) const;

    // Returns a (normalized) pointing vector from the camera center
    //  through the position of the pixel 'pix' on the image plane.
    virtual Vector3 pixel_to_vector (Vector2 const& pix) const;

    virtual Vector3 camera_center(Vector2 const& /*pix*/ = Vector2() ) const {
      return m_camera_center;
    };
    void set_camera_center(Vector3 const& position) {
      m_camera_center = position; rebuild_camera_matrix();
    }

    // Pose is a rotation which moves a vector in camera coordinates
    // into world coordinates.
    virtual Quaternion<double> camera_pose(Vector2 const& /*pix*/ = Vector2() ) const {
      return Quaternion<double>(m_rotation);
    }
    void set_camera_pose(Quaternion<double> const& pose) {
      m_rotation = pose.rotation_matrix(); rebuild_camera_matrix();
    }
    void set_camera_pose(Matrix<double,3,3> const& pose) {
      m_rotation = pose; rebuild_camera_matrix();
    }

    //  u_direction, v_direction, and w_direction define how the coordinate
    //  system of the camera relate to the directions in the image:
    //  +u (increasing image columns),
    //  +v (increasing image rows), and
    //  +w (pointing away from the focal point in the direction of the imaged object).
    //
    //  All three vectors must be of unit length.

    void coordinate_frame(Vector3 &u_vec, Vector3 &v_vec, Vector3 &w_vec) const {
      u_vec = m_u_direction;
      v_vec = m_v_direction;
      w_vec = m_w_direction;
    }

    void set_coordinate_frame(Vector3 u_vec, Vector3 v_vec, Vector3 w_vec) {
      m_u_direction = u_vec;
      m_v_direction = v_vec;
      m_w_direction = w_vec;

      rebuild_camera_matrix();
    }

    // Redundant...
    Vector3 coordinate_frame_u_direction() const { return m_u_direction; }
    Vector3 coordinate_frame_v_direction() const { return m_v_direction; }
    Vector3 coordinate_frame_w_direction() const { return m_w_direction; }

    const LensDistortion* lens_distortion() const { return m_distortion.get(); };
    void set_lens_distortion(LensDistortion const& distortion) {
      m_distortion = distortion.copy();
    }

    //  f_u and f_v :  focal length in horiz and vert. pixel units
    //  c_u and c_v :  principal point in pixel units
    void intrinsic_parameters(double& f_u, double& f_v,
                              double& c_u, double& c_v) const VW_DEPRECATED;
    void set_intrinsic_parameters(double f_u, double f_v,
                                  double c_u, double c_v) VW_DEPRECATED;

    Vector2 focal_length() const { return Vector2(m_fu,m_fv); }
    void set_focal_length(Vector2 const& f, bool rebuild=true ) {
      m_fu = f[0]; m_fv = f[1];
      if (rebuild) rebuild_camera_matrix();
    }
    Vector2 point_offset() const { return Vector2(m_cu,m_cv); }
    void set_point_offset(Vector2 const& c, bool rebuild=true ) {
      m_cu = c[0]; m_cv = c[1];
      if (rebuild) rebuild_camera_matrix();
    }
    double pixel_pitch() const { return m_pixel_pitch; }
    void set_pixel_pitch( double pitch ) { m_pixel_pitch = pitch; }

    // Ingest camera matrix
    // This performs a camera matrix decomposition and rewrites most variables
    void set_camera_matrix( Matrix<double,3,4> const& p );

    Matrix<double,3,4> camera_matrix() const {
      return m_camera_matrix;
    }

  private:
    /// This must be called whenever camera parameters are modified.
    void rebuild_camera_matrix();
  };

  //   /// Given two pinhole camera models, this method returns two new camera
  //   /// models that have been epipolar rectified.
  //   template <>
  //   void epipolar(PinholeModel<NoLensDistortion> const& src_camera0,
  //                 PinholeModel<NoLensDistortion> const& src_camera1,
  //                 PinholeModel<NoLensDistortion> &dst_camera0,
  //                 PinholeModel<NoLensDistortion> &dst_camera1);

  PinholeModel scale_camera(PinholeModel const& camera_model,
                            float const& scale);
  PinholeModel linearize_camera(PinholeModel const& camera_model);

  std::ostream& operator<<(std::ostream& str, PinholeModel const& model);

}}      // namespace vw::camera

#endif  //__CAMERAMODEL_CAHV_H__
