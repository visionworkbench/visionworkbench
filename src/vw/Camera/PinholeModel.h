// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file PinholeModel.h
///
/// This file contains the pinhole camera model.
///
#ifndef __VW_CAMERAMODEL_PINHOLE_H__
#define __VW_CAMERAMODEL_PINHOLE_H__

#include <vw/Math/Quaternion.h>
#include <vw/Camera/CameraModel.h>

#include <iostream>
#include <fstream>

namespace vw {
namespace camera {

  class LensDistortion;

  /// This is a simple "generic" pinhole camera model.
  ///
  /// To specify the EXTRINSIC parameters of the camera, we specify the
  /// position of the camera center in the world frame (m_camera_center)
  /// and the pose (or orientation) of the camera in world frame (m_rotation)
  /// (which is the transformation from the camera's frame to the world frame).
  /// In the default Vision Workbench camera frame, the camera's pointing vector
  /// is the +z unit vector, and the image plane is aligned such that the
  /// positive x-pixel direction (increasing image columns) is the camera frame's
  /// +x vector, and the positive y-pixel direction (increasing image
  /// rows) is the frame's y vector.
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
  ///    [  fx   0    cx  ]
  /// K= [  0    fy   cy  ]
  ///    [  0    0    1   ]
  ///
  /// with fx, fy the focal length of the system (in horizontal and
  /// vertical pixels), and (cx, cy) the pixel offset of the
  /// principal point of the camera on the image plane.
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
  ///  Lens distortion models can be found in the file LensDistortion.h
  ///

  class PinholeModel : public CameraModel {
    typedef boost::shared_ptr<LensDistortion> DistortPtr;

    DistortPtr m_distortion;
    Matrix<double,3,4> m_camera_matrix;

    // Stored for easy access.
    Vector3 m_camera_center;
    Matrix<double,3,3> m_rotation;
    Matrix<double,3,3> m_intrinsics;
    Matrix<double,3,4> m_extrinsics;

    /// Intrinsic parameters, in pixel units
    /// Focal length in u and v, 
    double m_fu, m_fv, m_cu, m_cv;

    // Vectors that define how the coordinate system of the camera
    // relate to the directions: +u (increasing image columns), +v
    // (increasing image rows), and +w (out along optical axis)
    Vector3 m_u_direction;
    Vector3 m_v_direction;
    Vector3 m_w_direction;

    // Pixel Pitch, if the above units were not in pixels this should
    // convert it to that. For example, if distortion and focal length
    // have been described in mm. Pixel pitch would then be described in mm/px.
    // - To clarify, if the above units are given in pixels this should equal to 1.0.
    // - It is kind of annoying that m_cu and m_cv must be in the same units as 
    //   all the numbers in the lens distortion parameters.
    double m_pixel_pitch;

    /// If true, perform an additional check to make sure point_to_pixel
    /// does not return erroneous results (but this can be slow).
    /// - Default is true.
    bool m_do_point_to_pixel_check;

    /// Cached values for pixel_to_vector
    Matrix<double,3,3> m_inv_camera_transform;

  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------

    /// Initialize an empty camera model.
    PinholeModel();

    /// Initialize from a file on disk.
    PinholeModel(std::string const& filename);

    /// Copy constructor. A deep copy is made of the distortion model held by a pointer.
    PinholeModel(PinholeModel const& other);

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
                 LensDistortion const* distortion_model = 0,
                 double pixel_pitch = 1.0);

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
    ///   +v (increasing image rows)                        =  +Y   [0 1 0]
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
                 LensDistortion const* distortion_model = 0,
                 double pixel_pitch = 1.0);

    virtual std::string type() const { return "Pinhole"; }
    virtual ~PinholeModel() {}

    /// Read / Write a pinhole model from a file on disk.
    /// - Supported file formats are .pinhole, .tsai
    void read (std::string const& filename);
    void write(std::string const& filename) const;

    //------------------------------------------------------------------
    // Methods
    //------------------------------------------------------------------

    //  Computes the image of the point 'point' in 3D space on the
    //  image plane.  Returns a pixel location (col, row) where the
    //  point appears in the image.
    virtual Vector2 point_to_pixel(Vector3 const& point) const;

    /// Skips the pixel_to_vector call used for a sanity check in point_to_pixel.
    Vector2 point_to_pixel_no_check(Vector3 const& point) const;

    /// Option to turn off and on the point to pixel check (default is on).
    void set_do_point_to_pixel_check(bool value) {m_do_point_to_pixel_check = value;}

    /// As point_to_pixel, but ignoring any lens distortion.
    Vector2 point_to_pixel_no_distortion(Vector3 const& point) const;

    // Is a valid projection of point is possible?
    // This is equal to: Is the point in front of the camera (z > 0)
    // after extrinsic transformation?
    virtual bool projection_valid(Vector3 const& point) const;

    // Returns a (normalized) pointing vector from the camera center
    //  through the position of the pixel 'pix' on the image plane.
    virtual Vector3 pixel_to_vector (Vector2 const& pix) const;

    // The pinhole camera position does not vary by pixel so the input pixel is ignored.
    virtual Vector3 camera_center(Vector2 const& /*pix*/ = Vector2() ) const;
    void set_camera_center(Vector3 const& position);

    // Pose is a rotation which moves a vector in camera coordinates
    // into world coordinates.
    // - The pinhole camera position does not vary by pixel so the input pixel is ignored.
    virtual Quaternion<double> camera_pose(Vector2 const& /*pix*/ = Vector2() ) const;
    void set_camera_pose(Quaternion<double> const& pose);
    void set_camera_pose(Matrix<double,3,3> const& pose);
    
    // WARNING: There may be a bug copying pose between cameras so use this function
    //          instead of camera_pose() for copies until it is resolved!
    Matrix<double,3,3> const& get_rotation_matrix() const {return m_rotation;}

    //  u_direction, v_direction, and w_direction define how the coordinate
    //  system of the camera relate to the directions in the image:
    //  +u (increasing image columns),
    //  +v (increasing image rows), and
    //  +w (pointing away from the focal point in the direction of the imaged object).
    //
    //  All three vectors must be of unit length.
    void coordinate_frame(Vector3 &u_vec, Vector3 &v_vec, Vector3 &w_vec) const;
    void set_coordinate_frame(Vector3 u_vec, Vector3 v_vec, Vector3 w_vec);

    // Redundant...
    Vector3 coordinate_frame_u_direction() const;
    Vector3 coordinate_frame_v_direction() const;
    Vector3 coordinate_frame_w_direction() const;

    LensDistortion const* lens_distortion() const;
    void set_lens_distortion(LensDistortion const* distortion); // Makes a copy

    //  f_u and f_v :  focal length in horiz and vert. pixel units
    //  c_u and c_v :  principal point in pixel units
    void intrinsic_parameters(double& f_u, double& f_v,
                              double& c_u, double& c_v) const;
    void set_intrinsic_parameters(double f_u, double f_v,
                                  double c_u, double c_v);

    Vector2 focal_length() const;
    Vector2 point_offset() const;
    double  pixel_pitch () const;
    
    void set_focal_length(Vector2 const& f, bool rebuild=true );
    void set_point_offset(Vector2 const& c, bool rebuild=true );
    void set_pixel_pitch (double pitch);

    // Ingest camera matrix
    // This performs a camera matrix decomposition and rewrites most variables
    void set_camera_matrix( Matrix<double,3,4> const& p );

    Matrix<double,3,4> camera_matrix() const;

    // Apply a given rotation + translation + scale transform to a pinhole camera
    void apply_transform(vw::Matrix4x4 const & transform);
    void apply_transform(vw::Matrix3x3 const & rotation,
                         vw::Vector3   const & translation,
                         double                scale);

  private:
    /// This must be called whenever camera parameters are modified.
    void rebuild_camera_matrix();
    
    /// Initialize m_distortion with the correct type of lens distortion
    ///  model depending on a string from an input .tsai file.
    /// - Returns true if it found a name or false if it did not and
    ///   just created the default TSAI distortion model.
    bool construct_lens_distortion(std::string const& config_line, int camera_version);
  };

  // TODO: Any use for an epipolar alignment function that operates on pinhole cameras?

  /// Used to modify camera in the event to user resizes the image
  /// - Under the hood all this does is change the pixel pitch.
  PinholeModel scale_camera(PinholeModel const& camera_model, double scale);
                            
  /// Returns a copy of the camera model with no lens distortion.
  /// - This does not account for the distortion in any way, 
  ///   it just removes the distortion model!
  PinholeModel strip_lens_distortion(PinholeModel const& camera_model);

  /// From a pair of input pinhole cameras, produce a pair of associated
  ///  zero-distortion epipolar aligned cameras which can be used for stereo.
  void epipolar(PinholeModel const &src_camera0, PinholeModel const &src_camera1,
                PinholeModel       &dst_camera0, PinholeModel       &dst_camera1);

  /// Write a description of a PinholeModel to the stream.
  std::ostream& operator<<(std::ostream& str, PinholeModel const& model);

}}      // namespace vw::camera

#endif  //__CAMERAMODEL_PINHOLE_H__
