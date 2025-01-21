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

/// \file CameraParamsPack.cc
///
/// Minor functions for packing camera params to a vector and vice versa.
/// Some of these are templated so need to be in the .h file.
///

#include <vw/Camera/CameraParamsPack.h>

namespace vw {
namespace camera {

// Unpack a vector into a rotation + translation + scale
void vector_to_transform(vw::Vector<double> const & C,
                         vw::Matrix3x3            & rotation,
                         vw::Vector3              & translation,
                         double                   & scale) {

  if (C.size() != 7)
    vw_throw(LogicErr() << "Expecting a vector of size 7.\n");

  translation = vw::math::subvector(C, 0, 3);
  rotation    = vw::math::axis_angle_to_quaternion(subvector(C, 3, 3)).rotation_matrix();
  scale       = C[6];

  return;
}

// Pack a rotation + translation + scale into a vector
void transform_to_vector(vw::Vector<double>  & C,
                         vw::Matrix3x3 const & rotation,
                         vw::Vector3   const & translation,
                         double    const & scale) {

  vw::Vector3 axis_angle = vw::Quat(rotation).axis_angle();

  C.set_size(7);
  vw::math::subvector(C, 0, 3) = translation;
  vw::math::subvector(C, 3, 3) = axis_angle;
  C[6] = scale;

  return;
}

}} // end namespace vw::camera
