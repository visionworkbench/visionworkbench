// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestConjugateGradient.h
#include <cxxtest/TestSuite.h>
#include <vw/Math/Vector.h>
#include <vw/Math/ConjugateGradient.h>

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


class TestConjugateGradient : public CxxTest::TestSuite
{
public:

  void test_steepest_descent()
  {
    //    set_debug_level(VerboseDebugMessage);
    Vector2 initial_guess(2,2);
    int numiters = 100;
    int max_stepsize = 1;
    QuadraticFunction cost_functor;
    QuadraticFunction::domain_type result = steepest_descent( cost_functor, initial_guess, ArmijoStepSize(max_stepsize), numiters);
//    TS_TRACE(stringify(result));
//    TS_TRACE(stringify(cost_functor(result)))
//    TS_TRACE(stringify(cost_functor.gradient(result)))

    TS_ASSERT_DELTA(result[0], 0.1962, 0.001);
    TS_ASSERT_DELTA(result[1], 0.4846, 0.001);
  }

  void test_conjugate_gradient()
  {
    //    set_debug_level(VerboseDebugMessage);
    Vector2 initial_guess(2,2);
    int numiters = 100;
    int max_stepsize = 1;
    QuadraticFunction cost_functor;
    QuadraticFunction::domain_type result = conjugate_gradient( cost_functor, initial_guess, ArmijoStepSize(max_stepsize), numiters);
//     TS_TRACE(stringify(result));
//     TS_TRACE(stringify(cost_functor(result)));
//     TS_TRACE(stringify(cost_functor.gradient(result)));

    TS_ASSERT_DELTA(result[0], 0.1962, 0.001);
    TS_ASSERT_DELTA(result[1], 0.4846, 0.001);
  }

}; // class TestVector
