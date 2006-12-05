// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

// TestOptimization.h
#include <cxxtest/TestSuite.h>
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

class TestOptimization : public CxxTest::TestSuite
{
public:

  void test_least_squares_model()
  {
    TestLeastSquaresModel model;
    Vector<double> z(5);
    Vector<double> x(4);
    x(0) = 0.2;
    x(1) = 0.3;
    x(2) = 0.4;
    x(3) = 0.5;
    
    TS_ASSERT_DELTA(model(x)(0), 0.29552, 0.0001);
    TS_ASSERT_DELTA(model(x)(1), 0.992809, 0.0001);
    TS_ASSERT_DELTA(model(x)(2), 0.276318, 0.0001);
    TS_ASSERT_DELTA(model(x)(3), 0.380506, 0.0001);
    TS_ASSERT_DELTA(model(x)(4), 0.927295, 0.0001);

    TS_ASSERT_DELTA(model.jacobian(x)(1,1),-0.0478849,0.0001); 
    TS_ASSERT_DELTA(model.jacobian(x)(2,2),-0.116826,0.0001); 
    TS_ASSERT_DELTA(model.jacobian(x)(3,0),1.72414,0.0001); 
   }

  void test_levenberg_marquardt()
  {
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
    TS_ASSERT_EQUALS(status, vw::math::optimization::eConvergedRelTolerance);
    TS_ASSERT_DELTA(best(0),0.101358,0.0001);
    TS_ASSERT_DELTA(best(1),1.15485,0.0001);
    TS_ASSERT_DELTA(best(2),1.12093,0.0001);
    TS_ASSERT_DELTA(best(3),0.185534,0.0001);
  }

}; // class TestOptimization
