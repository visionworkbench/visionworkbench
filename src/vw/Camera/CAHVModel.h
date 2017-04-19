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


/// \file CAHVModel.h
///
/// This file contains the CAHV pinhole camera model.
///
#ifndef __VW_CAMERAMODEL_CAHV_H__
#define __VW_CAMERAMODEL_CAHV_H__

#include <vw/Camera/CameraModel.h>

namespace vw {
namespace camera {

  class PinholeModel;

  /// The CAHV pinhole camera model has been widely used in NASA
  /// planetary mission for rover navigation and scientific camera
  /// systems.  In the CAHV model, camera intrinsics and extrincsics
  /// are jointly parameterized by four 3-Vectors C,A,H, and V.  The
  /// strength of this model is the mathematical ease with which 3D
  /// points can be projected onto the image plane.  The CAHV model is
  /// computationally efficient, however it does not model lens
  /// distortion.  For this, its close relative the CAHVOR model
  /// should be used.  Refer to vw::camera::CAHVORModel for more
  /// information.
  ///
  /// References:
  ///
  /// Yakimovsky, Y. and Cunningham R., "A System for Extracting
  /// Three-Dimensional Measurements from a Stereo Pair of TV
  /// Cameras". Computer Graphics and Image Processing 7,
  /// pp. 195-210. (1978)
  ///
  class CAHVModel : public CameraModel {
  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------

    /// Initialize an empty camera model where the 4 CAHV vectors are
    /// set to zero.
    CAHVModel() {}

    /// This constructor takes a filename and reads in a camera model
    /// from the file.  The file may contain either CAHV parameters or
    /// pinhole camera parameters.
    CAHVModel(std::string const& filename);

    /// Initialize the CAHV vectors directly in the native CAHV format.
    CAHVModel(Vector3 const& C_vec,
              Vector3 const& A_vec,
              Vector3 const& H_vec,
              Vector3 const& V_vec);

    /// Initialize a CAHV model from a pinhole model
    /// - Since the CAHV model does not support lens distortion, any lens
    ///   distortion information is simply lost.
    CAHVModel(PinholeModel const& pin_model);

    CAHVModel operator= (PinholeModel const& pin_model);

    virtual std::string type() const;

    /// Initialize the CAHV vectors indirectly using pinhole camera
    /// parameters.  In this variant, the view matrix is supplied directly.
    ///
    /// f         - focal length in world units
    /// pixel_size - (width, height) of a pixel on the image plane in world units
    /// xMin, xMax, yMin, yMax - image plane boundaries in units of pixels
    ///
    /// View matrix transformation:
    ///
    /// |R00 R01 R02 Tx|
    /// |R10 R11 R12 Ty|
    /// |R20 R21 R22 Tz|
    /// | 0   0   0   1|
    ///
    CAHVModel(double f, Vector2 const& pixel_size,
              double xmin, double /*xmax*/, double ymin, double /*ymax*/,
              math::Matrix<double, 4, 4> const& view_matrix);

    /// Initialize the CAHV vectors indirectly using pinhole camera parameters:
    ///
    /// f         - focal length in world units
    /// pixel_size - (width, height) of a pixel on the image plane in world units
    /// Hc, Vc    - location of projection of C on the image plane in units of
    ///             pixels from lower left hand corner
    ///
    /// C         - center of projection in world coordinates
    /// A         - unit vector perpendicular to image plane in world coords
    ///             pointing out the "front" of the camera
    /// Hvec      - unit horizontal image plane vector in world coords
    /// Vvec      - unit vertical image plane vector in world coords
    ///
    CAHVModel(double f, Vector2 const& pixel_size, double Hc, double Vc,
              Vector3 const& Cinit, Vector3 const& Ainit,
              Vector3 const& Hvec,  Vector3 const& Vvec);

    virtual ~CAHVModel() {}

    //------------------------------------------------------------------
    // Methods
    //------------------------------------------------------------------
    virtual Vector2 point_to_pixel (Vector3 const& point) const;
    virtual Vector3 pixel_to_vector(Vector2 const& pix  ) const;
    virtual Vector3 camera_center  (Vector2 const& /*pix*/ = Vector2() ) const;

    /// Write CAHV model to file
    void write(std::string const& filename);

    //------------------------------------------------------------------
    // Exposed Variables
    //------------------------------------------------------------------
    Vector3   C; //< The pinhole center coordinate.
    Vector3   A; //< A unit vector normal to the sensor plane.
    Vector3   H; //< Horizontal information vector
    Vector3   V; //< Vertical information vector

  private:
    void read_cahv   (std::string const& filename);
    void read_pinhole(std::string const& filename);
  };

  /// Given two CAHV camera models, this method returns two new camera
  /// models that have been epipolar rectified.
  void epipolar(CAHVModel const &src_camera0, CAHVModel const &src_camera1,
                CAHVModel       &dst_camera0, CAHVModel       &dst_camera1);

  std::ostream& operator<<(std::ostream& str, CAHVModel const& model);

}} // namespace vw::camera

#endif //__CAMERAMODEL_CAHV_H__
