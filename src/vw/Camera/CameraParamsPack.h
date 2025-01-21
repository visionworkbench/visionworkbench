// __BEGIN_LICENSE__
//  Copyright (c) 2006-2025, United States Government as represented by the
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

/// \file CameraParamsPack.h
///
/// Minor functions for packing camera params to a vector and vice versa.
/// Some of these are templated so need to be in the .h file.

#ifndef __VW_CAMERA_PARAMS_PACK_H__
#define __VW_CAMERA_PARAMS_PACK_H__

#include <vw/Camera/CameraModel.h>

namespace vw {
namespace camera {

// Pack a pinhole or optical bar model's rotation and camera center into a vector
template <class CAM>
void camera_to_vector(CAM const& P, vw::Vector<double> & C) {

  vw::Vector3 ctr = P.camera_center();
  vw::Vector3 axis_angle = P.camera_pose().axis_angle();
  C.set_size(6);
  vw::math::subvector(C, 0, 3) = ctr;
  vw::math::subvector(C, 3, 3) = axis_angle;
}

// Pack a vector into a pinhole or optical bar model. We assume the model
// already has set its optical center, focal length, etc.
template <class CAM>
void vector_to_camera(CAM & P, vw::Vector<double> const& C) {

  if (C.size() != 6)
    vw::vw_throw(vw::LogicErr() << "Expecting a vector of size 6.\n");

  vw::Vector3 ctr = subvector(C, 0, 3);
  vw::Quat rotation = vw::math::axis_angle_to_quaternion(subvector(C, 3, 3));
  P.set_camera_center(ctr);
  P.set_camera_pose(rotation);
}

// Unpack a vector into a rotation + translation + scale
void vector_to_transform(vw::Vector<double> const & C,
                         vw::Matrix3x3            & rotation,
                         vw::Vector3              & translation,
                         double                   & scale);

// Pack a rotation + translation + scale into a vector
void transform_to_vector(vw::Vector<double>  & C,
                         vw::Matrix3x3 const & rotation,
                         vw::Vector3   const & translation,
                         double    const & scale);

}} // namespace vw::camera

#endif  //__CAMERA_PARAMS_PACK_H__
