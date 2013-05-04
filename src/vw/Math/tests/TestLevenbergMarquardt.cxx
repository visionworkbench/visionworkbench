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
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/LevenbergMarquardt.h>

using namespace vw;
using namespace vw::math;

struct TestLeastSquaresModel : public LeastSquaresModelBase<TestLeastSquaresModel> {

  typedef Vector<double> result_type;
  typedef Vector<double> domain_type;
  typedef Matrix<double> jacobian_type;

  TestLeastSquaresModel() {}

  /// Evaluate h(x)
  inline result_type operator()( domain_type const& x ) const {

    // For now make up a function to get started with
    Vector<double> h(5);
    h(0) = sin(x(0)+0.1);
    h(1) = cos(x(1) * x(2));
    h(2) = x(1) * cos(x(2));
    h(3) = atan2(x(0),x(3));
    h(4) = atan2(x(2),x(1));
    return h;
  }
};

TEST(LevenbergMarquardt, least_squares_model) {

  typedef Vector<double,5> Vector5;
  TestLeastSquaresModel model;

  Vector5 z;
  Vector4 x( 0.2, 0.3, 0.4, 0.5 );
  Vector5 expected_first(
      .29552020666133957510,
      .99280863585386625224,
      .27631829820086552483,
      .38050637711236488630,
      .92729521800161223242);

  Vector<double> first = model(x);

  EXPECT_VECTOR_DOUBLE_EQ( expected_first, first );

  Matrix<double> jacobian = model.jacobian(x);

  //TODO: Vet this.
  EXPECT_NEAR( -0.0478849, jacobian(1,1), 1e-5 );
  EXPECT_NEAR( -0.116826,  jacobian(2,2), 1e-5 );
  EXPECT_NEAR(  1.72414,   jacobian(3,0), 1e-5 );
}

TEST(LevenbergMarquardt, levenberg_marquardt) {
  TestLeastSquaresModel model;
  Vector<double> target(5);
  target(0) = 0.2;
  target(1) = 0.3;
  target(2) = 0.4;
  target(3) = 0.5;
  target(4) = 0.6;

  Vector<double> seed(4);
  seed[0] = 1.0;
  seed[1] = 1.0;
  seed[2] = 1.0;
  seed[3] = 1.0;

  int status;
  Vector<double> best = levenberg_marquardt( model, seed, target, status );

  //TODO: Vet this.
  Vector4 expected_best( 0.101358, 1.15485, 1.12093, 0.185534 );

  EXPECT_EQ(vw::math::optimization::eConvergedRelTolerance, status);
  EXPECT_VECTOR_NEAR( expected_best, best, 1e-5 );
}
