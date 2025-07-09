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

#include <vw/Math/Quaternion.h>

namespace vw {
namespace math {

Quat slerp(double alpha, Quat const& a, Quat const& b, int spin) {
  const double SLERP_EPSILON = 1.0E-6;              // a tiny number
  double beta;                      // complementary interp parameter
  double theta;                     // angle between A and B
  double sin_t, cos_t;              // sine, cosine of theta
  double phi;                       // theta plus spins
  int bflip;                        // use negation of B?

  // cosine theta = dot product of A and B
  cos_t = a(1)*b(1) + a(2)*b(2) + a(3)*b(3) + a(0)*b(0);

  // if B is on opposite hemisphere from A, use -B instead
  if (cos_t < 0.0) {
    cos_t = -cos_t;
    bflip = true;
  } else {
    bflip = false;
  }

  // if B is (within precision limits) the same as A,
  // just linear interpolate between A and B.
  // Can't do spins, since we don't know what direction to spin.
  if (1.0 - cos_t < SLERP_EPSILON) {
    beta = 1.0 - alpha;
  } else {                          /* normal case */
    theta = acos(cos_t);
    phi = theta + spin * M_PI;
    sin_t = sin(theta);
    beta = sin(theta - alpha*phi) / sin_t;
    alpha = sin(alpha*phi) / sin_t;
  }

  if (bflip)
    alpha = -alpha;

  // interpolate
  return Quat(beta*a(0) + alpha*b(0),
              beta*a(1) + alpha*b(1),
              beta*a(2) + alpha*b(2),
              beta*a(3) + alpha*b(3));
}

// Compute the n-weight slerp, analogous to the linear combination
// w[0]*Q[0] + ... + w[n-1]*Q[n-1]. This is experimental.
// We assume the sum of weights is 1.
Quat slerp_n(std::vector<double> const& w, std::vector<Quat> const& Q, int spin) {
  VW_ASSERT(w.size() == Q.size(),
	    ArgumentErr() << "Expecting as many quaternions as weights.\n");

  VW_ASSERT(Q.size() >= 1,
	    ArgumentErr() << "Expecting at least one quaternion and weight.\n");

  if (Q.size() == 1)
    return Q[0];

  if (Q.size() == 2) {
    VW_ASSERT(std::abs(w[0] + w[1] - 1.0) < 1e-6 && w[0] >= 0 && w[1] >= 0,
	    ArgumentErr() << "Expecting the weights to be >= 0 and sum up to 1.\n");
    return slerp(w[1], Q[0], Q[1], spin);
  }

  // Call recursively this function with fewer terms

  double sum = w[0] + w[1];
  if (sum == 0) sum = 1.0;

  Quat q = slerp(w[1]/sum, Q[0], Q[1], spin);
  std::vector<double> w2 = w;
  std::vector<Quat>   Q2 = Q;
  w2.erase(w2.begin());
  Q2.erase(Q2.begin());
  w2[0] = sum;
  Q2[0] = q;
  return slerp_n(w2, Q2, spin);
}

}} // namespace vw::math
