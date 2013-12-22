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


/// \file CAHVOREModel.h
///
/// This file contains the CAHVORE camera model.  This camera
/// model is a refinement of the \ref vw::camera::CAHVORModel: it adds extra terms
/// to model fisheye lenses
///
#ifndef __VW_CAMERAMODEL_CAHVORE_H__
#define __VW_CAMERAMODEL_CAHVORE_H__

#include <vw/Core/Log.h>
#include <vw/Camera/CameraModel.h>

namespace vw {
namespace camera {

  class CAHVModel;

  class CAHVOREModel : public CameraModel {
  public:
    typedef CAHVModel linearized_type;

    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------

    /// Default constructor creates CAHVORE vectors equal to 0.
    CAHVOREModel() {}

    /// Read a CAHVORE file from disk.
    CAHVOREModel(std::string const& filename);

    /// Initialize the CAHVORE vectors directly in the native CAHVORE format.
    CAHVOREModel(Vector3 const& C_vec, Vector3 const& A_vec,
                 Vector3 const& H_vec, Vector3 const& V_vec,
                 Vector3 const& O_vec, Vector3 const& R_vec,
                 Vector3 const& E_Vec);
    CAHVOREModel(Vector3 const& C_vec, Vector3 const& A_vec,
                 Vector3 const& H_vec, Vector3 const& V_vec,
                 Vector3 const& O_vec, Vector3 const& R_vec,
                 Vector3 const& E_Vec, int T, double P_v);
    CAHVOREModel(Vector3 const& C_vec, Vector3 const& A_vec,
                 Vector3 const& H_vec, Vector3 const& V_vec,
                 Vector3 const& O_vec, Vector3 const& R_vec,
                 Vector3 const& E_Vec, double P_v);

    virtual ~CAHVOREModel();
    virtual std::string type() const;

    //------------------------------------------------------------------
    // Methods
    //------------------------------------------------------------------
    virtual Vector2 point_to_pixel(Vector3 const& point) const;
    virtual Vector3 pixel_to_vector(Vector2 const& pix) const;
    virtual Vector3 camera_center(Vector2 const& /*pix*/ = Vector2() ) const;

    /// Write CAHVORE model to file.
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
    Vector3   E;
    double    P; // We don't have T as it is redundant information
  private:
    bool check_line( std::istream& istream, char letter );
  };

  // Function to "map" the CAHVORE parameters into CAHV:
  CAHVModel linearize_camera( CAHVOREModel const& camera_model,
                              Vector2i const& cahvore_image_size,
                              Vector2i const& cahv_image_size );

}} // namespace vw::camera

#endif//__VW_CAMERAMODEL_CAHVORE_H__
