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

/// \file LinescanErr.h
///
#ifndef __CAMERA_LINESCAN_ERR_H__
#define __CAMERA_LINESCAN_ERR_H__

#include <vw/Camera/CameraModel.h>

namespace vw {

namespace camera {

  // Find the projection of a ground point into a linescan camera by 
  // employing the 2D to 2D Newton-Raphson method in a plane that is
  // approximately at the ground level.

class LinescanErr {

  const CameraModel* m_model;
  Vector3 m_point;
  
  // Two vectors that define the ground plane
  vw::Vector3 m_perp1, m_perp2;

public:

  // Constructor
  LinescanErr(const CameraModel* model, const vw::Vector3& pt, vw::Vector2 const& guess);

  // This must have the signature expected by Newton's method. Can throw
  // exceptions. The math is described above.
  vw::Vector2 operator()(vw::Vector2 const& pix) const;
  
}; // End class LinescanErr

}} // end vw::camera

#endif // __CAMERA_LINESCAN_ERR_H__
