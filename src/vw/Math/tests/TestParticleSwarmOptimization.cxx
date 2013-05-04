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


#include <test/Helpers.h>
#include <vw/Math/ParticleSwarmOptimization.h>

using namespace vw;
using namespace vw::math;

// This quadratic function has a single minimum at [0.1962, 0.4846].
struct QuadraticFunction {
  typedef double result_type;
  typedef Vector2 domain_type;
  typedef Vector2 gradient_type;

  result_type operator()( domain_type const& x ) const {
    return 1.2 * pow(x[0] - 0.6, 2) + 1.7 * pow(x[1] - 0.6, 2) + 2 * x[0] * x[1];
  }
  gradient_type gradient( domain_type const& x ) const {
    return Vector2( 2.4*x[0]-1.44+2*x[1],
                    3.4*x[1]-2.04+2*x[0]);
  }

  unsigned dimension() const { return 2; }
};

struct SinFunction {
  typedef double result_type;
  typedef Vector4 domain_type;
  typedef Vector4 gradient_type;

  result_type operator()( domain_type const& x ) const {
    return std::sin(x[0])*std::sin(x[1])*std::sin(x[2])*std::sin(x[3]);
  }

  gradient_type gradient( domain_type const& x ) const {
    return Vector4( std::cos(x[0])*std::sin(x[1])*std::sin(x[2])*std::sin(x[3]),
                    std::sin(x[0])*std::cos(x[1])*std::sin(x[2])*std::sin(x[3]),
                    std::sin(x[0])*std::sin(x[1])*std::cos(x[2])*std::sin(x[3]),
                    std::sin(x[0])*std::sin(x[1])*std::sin(x[2])*std::cos(x[3]) );
  }

  unsigned dimension() const { return 4; }
};

inline double modulo(double v, double d) {
  return (v/d - static_cast<int>(v/d))*d;
}

TEST(ParticleSwarmOptimization, quadratic) {
  Vector2 min(-2, -2);
  Vector2 max(2, 2);

  QuadraticFunction::domain_type result = particle_swarm_optimization( QuadraticFunction(), min, max );
  EXPECT_VECTOR_NEAR( Vector2(0.1962, 0.4846), result, 1e-2 );
}

TEST(ParticleSwarmOptimization, quadratic_large_search) {
  Vector2 min(-200, -200);
  Vector2 max(200, 200);

  QuadraticFunction::domain_type result = particle_swarm_optimization( QuadraticFunction(), min, max, false, 5, 1000, 2000);
  EXPECT_VECTOR_NEAR( Vector2(0.1962, 0.4846), result, 1e-3 );
}

TEST(ParticleSwarmOptimization, sinus_function) {
  Vector4 min(-200, -200, -200, -200);
  Vector4 max(200, 200, 200, 200);

  SinFunction cost_functor;
  SinFunction::domain_type result = particle_swarm_optimization( cost_functor, min, max, false, 2, 100, 2000);

  EXPECT_NEAR( -1.0, cost_functor(result), 2.7e-3 );

  double pi2 = static_cast<double>(M_PI)/2;

  // TODO: Trig precision is so bad...
  EXPECT_NEAR( pi2, std::fabs(modulo(result(0), M_PI)), 1e-1 );
  EXPECT_NEAR( pi2, std::fabs(modulo(result(1), M_PI)), 1e-1 );
  EXPECT_NEAR( pi2, std::fabs(modulo(result(2), M_PI)), 1e-1 );
  EXPECT_NEAR( pi2, std::fabs(modulo(result(3), M_PI)), 1e-1 );
}
