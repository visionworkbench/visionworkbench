// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

#include <vw/Math/VectorUtils.h>
#include <vw/Core/Exception.h>
#include <cmath>

namespace vw {
namespace math {

// Form an orthonormal basis from a vector
// Given input vector V, returns orthonormal basis vectors X, Y, Z where:
// - Z is aligned with V (normalized)
// - X and Y are perpendicular to Z and to each other
// - All three form a right-handed orthonormal coordinate system
void formBasis(Vector3 const& V, Vector3& X, Vector3& Y, Vector3& Z) {
  // Normalize V to get Z
  double norm = vw::math::norm_2(V);
  if (norm == 0.0)
    vw::vw_throw(vw::ArgumentErr() << "Cannot form basis from zero vector");

  Z = V / norm;

  // Find largest component by magnitude to ensure numerical stability
  int max_idx = 0;
  double max_abs = std::abs(Z[0]);
  for (int i = 1; i < 3; i++) {
    if (std::abs(Z[i]) > max_abs) {
      max_abs = std::abs(Z[i]);
      max_idx = i;
    }
  }

  // Create perpendicular vector by swapping components
  int idx1 = (max_idx + 1) % 3;
  int idx2 = (max_idx + 2) % 3;
  X = Vector3(0, 0, 0);
  X[idx1] = -Z[idx2];
  X[idx2] = Z[idx1];

  // Normalize X
  X = X / vw::math::norm_2(X);

  // Cross product to get Y, ensuring right-handed system
  Y = cross_prod(Z, X);
}

// Generate two orthogonal tangent vectors that lie in a plane defined by its normal.
// The vectors are scaled by the given offset distance.
// The tangent vectors are perpendicular to the normal and to each other.
void computePlaneTangents(Vector3 const& normal, double offset,
                          Vector3& tangent1, Vector3& tangent2) {
  
  if (std::abs(normal[2]) > 0.1) {
    // If plane isn't too vertical, use X and Y directions projected onto plane
    tangent1 = Vector3(1, 0, -normal[0]/normal[2]);
    tangent2 = Vector3(0, 1, -normal[1]/normal[2]);
  } else if (std::abs(normal[1]) > 0.1) {
    // Plane is nearly vertical in Z, use X and Z
    tangent1 = Vector3(1, -normal[0]/normal[1], 0);
    tangent2 = Vector3(0, -normal[2]/normal[1], 1);
  } else {
    // Plane is nearly vertical in Y and Z, use Y and Z
    tangent1 = Vector3(-normal[1]/normal[0], 1, 0);
    tangent2 = Vector3(-normal[2]/normal[0], 0, 1);
  }

  // Normalize tangent vectors and scale by offset
  tangent1 = normalize(tangent1) * offset;
  tangent2 = normalize(tangent2) * offset;
}

}} // namespace vw::math
