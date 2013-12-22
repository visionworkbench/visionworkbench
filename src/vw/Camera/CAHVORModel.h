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


/// \file CAHVORModel.h
///
/// This file contains the CAHVOR pinhole camera model.  This camera
/// model is a refinement of the \ref vw::camera::CAHVModel: it adds extra terms
/// to model basic radial camera distortion.
///
#ifndef __VW_CAMERAMODEL_CAHVOR_H__
#define __VW_CAMERAMODEL_CAHVOR_H__

#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Camera/CameraModel.h>

#include <string>

namespace vw {
namespace camera {

  class CAHVModel;

  class CAHVORModel : public CameraModel {
  public:
    typedef CAHVModel linearized_type;

    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------

    /// Default constructor creates CAHVOR vectors equal to 0.
    CAHVORModel() {}

    /// Read a CAHVOR file from disk.
    CAHVORModel(std::string const& filename);

    /// Initialize the CAHVOR vectors directly in the native CAHVOR format.
    CAHVORModel(Vector3 const& C_vec, Vector3 const& A_vec,
                Vector3 const& H_vec, Vector3 const& V_vec,
                Vector3 const& O_vec, Vector3 const& R_vec);

    virtual ~CAHVORModel();
    virtual std::string type() const;

    //------------------------------------------------------------------
    // Methods
    //------------------------------------------------------------------
    virtual Vector2 point_to_pixel(Vector3 const& point) const;
    virtual Vector3 pixel_to_vector(Vector2 const& pix) const;
    virtual Vector3 camera_center(Vector2 const& /*pix*/ = Vector2() ) const;

    // Overloaded versions also return partial derviatives in a Matrix.
    Vector2 point_to_pixel(Vector3 const& point, Matrix<double> &partial_derivatives) const;
    Vector3 pixel_to_vector(Vector2 const& pix, Matrix<double> &partial_derivatives) const;

    void get_point_derivatives( Vector3 const& P, double& u, double& v,
                                Vector3& grad_u, Vector3& grad_v,
                                Matrix3x3& hess_u, Matrix3x3& hess_v ) const;

    /// Write CAHVOR model to file.
    void write(std::string const& filename);

    //------------------------------------------------------------------
    // Exposed Variables
    //------------------------------------------------------------------
    Vector3   C;
    Vector3   A;
    Vector3   H;
    Vector3   V;
    Vector3   O;
    Vector3   R;
  };

  /// Function to "map" the CAHVOR parameters into CAHV parameters:
  /// requires dimensions of input image and output image (usually the
  /// same) You must supply the dimensions of the CAHVOR image that
  /// you are dewarping, as well as the dimensions of the dewarped
  /// cahv image (if they differ).
  CAHVModel linearize_camera( CAHVORModel const& camera_model,
                              Vector2i const& cavhor_image_size,
                              Vector2i const& cahv_image_size );

  CAHVModel linearize_camera( CAHVORModel const& camera_model,
                              int32 ix, int32 iy, int32 ox, int32 oy ) VW_DEPRECATED;

}}      // namespace vw::camera

#endif  //__CAMERAMODEL_CAHVOR_H__
