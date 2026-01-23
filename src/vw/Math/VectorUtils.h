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

/// \file VectorUtils.h
///
/// Vector utility functions.
///
#ifndef __VW_MATH_VECTOR_UTILS_H__
#define __VW_MATH_VECTOR_UTILS_H__

#include <vw/Math/Vector.h>

namespace vw {
namespace math {

  // Form an orthonormal basis from a vector
  // Given input vector V, returns orthonormal basis vectors X, Y, Z where:
  // - Z is aligned with V (normalized)
  // - X and Y are perpendicular to Z and to each other
  // - All three form a right-handed orthonormal coordinate system
  void formBasis(Vector3 const& V, Vector3& X, Vector3& Y, Vector3& Z);

  // Generate two orthogonal tangent vectors that lie in a plane defined by its normal.
  // The vectors are scaled by the given offset distance.
  // The tangent vectors are perpendicular to the normal and to each other.
  void computePlaneTangents(Vector3 const& normal, double offset,
                            Vector3& tangent1, Vector3& tangent2);

}} // namespace vw::math

#endif // __VW_MATH_VECTOR_UTILS_H__
