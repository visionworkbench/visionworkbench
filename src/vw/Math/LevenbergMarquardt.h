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


/// \file LevenbergMarquardt.h
///
/// Optimization classes for carrying out nonlinear optimization.
///
/// This code provides some optimization routines for doing estimation.
/// The goal is to provide some sufficiently general support to that
/// anything from Levenberg-Marquardt for ICP to Kalman filters for
/// navigation can use the same toolbox.

#ifndef __VW_MATH_OPTIMIZATION_H__
#define __VW_MATH_OPTIMIZATION_H__

// Vision Workbench
#include <vw/Core/Log.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/LinearAlgebra.h>

// Boost
#include <boost/concept_check.hpp>

namespace vw {
namespace math {

  /// First thing we need is a generic idea of a measurement function
  /// or model function.  The model function needs to provide a way to
  /// evaluate h(x) as well as a way to differentiate h() to get its
  /// Jacobian.  In the past we've implemented a fallback option by
  /// computing the Jacobian numerically which seems to work pretty
  /// well.  This LeastSquaresModelFunction shoulde be implemented as
  /// an inherited class using the curiously recurring template
  /// mechanism (CRTP) using the following syntax:
  ///
  /// class MyLeastSquaresModel : public LeastSquaresModelBase<MyLeastSquaresModel> {
  /// ...
  ///
  /// You will (at the very least) need to provide the following
  /// things in your base class:
  ///
  /// * Defines a result_type that is the type returned by
  ///   evaluating the functor.  Typically Vector<float> or
  ///   Vector<double>
  /// * Defines a domain_type that is the type of the search
  ///   space.  Often a Vector<foo>, but can reflect other
  ///   topologies if needed.
  /// * Defines a jacobian_type corresponding to the space of
  ///   jacobian matrices.  Typically Matrix<foo>.
  /// * Defines a method: result_type operator()( domain_type const& x ) const;
  ///   that evaluates the model function at the given point.
  /// * The domain_type must implement a method: domain_type domain_type::operator+( gradient_type const& g ) const;
  ///   that adds a tangent vector to a position.  You get this for
  ///   free if both domain_type and gradient_type are Vector<foo>.
  ///   This is where you do most of the hard work if domain_type
  ///   represents some non-trivial topological space.
  /// * The result_type must implement a method: double result_type::norm_2( result_type const& g ) const;
  ///   that is used in some optimizers to compute the error.  You get this for
  ///   free if result_type is a Vector<foo>.
  /// * The jacobian_type must implement several matrix-like methods such as
  ///   scalar multiplication on the left.  You get this for free in the usual case when
  ///   jacobian_type is just Matrix<foo>.
  ///
  /// In addition, depending on the application, you may want to also define:
  ///
  /// * Defines a method: jacobian_type jacobian( domain_type const& x ) const;
  ///   that evaluates the jacobian of the model function at the
  ///   given point.  This will override the default implementation in this base
  ///   class which computes the derivatives numerically.
  ///

  template <class ImplT>
  struct LeastSquaresModelBase {

    /// \cond INTERNAL
    // Methods to access the derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
    /// \endcond

    /// This default implementation evaluates dh(x;p)/dx at the
    /// supplied location using numerical differentiation.  It highly
    /// recommended that you override this method in your sub-class and
    /// provide an exact computation of the gradient.
    ///
    /// Note that this numerical derivative assumes that the
    /// result_type and jacobian_type types behave like vector and
    /// matrix objects, respectively.
    template <class DomainT>
    inline Matrix<double> jacobian( DomainT const& x ) const {

      // Get nominal function value
      Vector<double> h0 = impl().operator()(x);

      // Jacobian is #params x #outputs
      Matrix<double> H(h0.size(), x.size());

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( unsigned i=0; i<x.size(); ++i ){
        DomainT xi = x;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(xi(i)*1e-7);
        xi(i) += epsilon;

        // Evaluate function with this step and compute the derivative w.r.t. parameter i
        Vector<double> hi = impl().operator()(xi);
        select_col(H,i) = this->difference(hi,h0)/epsilon;
      }
      return H;
    }

    /// The utility of a differencing function is that it allows
    /// looping topologies like angles to be handled without putting
    /// the logic in the L-M implementation.  You probably won't need
    /// to override this method in your base class implementation
    /// unless you are working with unusual result_space topologies.
    template <class T>
    inline T difference( T const& a, T const& b ) const {
      return (a-b);
    }
  };

  namespace optimization {
    enum LM_STATUS_CODES { eDidNotConverge = -1,
                           eStatusUnknown = 0,
                           eConvergedAbsTolerance = 1,
                           eConvergedRelTolerance = 2 };
  }

  /// Levenberg-Marquardt is an algorithm for solving problems of the
  /// form
  ///
  /// J(p) = sum_i ( z_i - h(x_i) )^2
  ///
  /// That is, a least squares problem where the objective is to find a
  /// parameter vector x such that the model function h(x), evaluates
  /// as closely as possible to the observations z in a least squares
  /// sense, that is, the cost function J(p) is minimized.
  ///
  /// Requires:
  /// - a seed parameter vector x
  /// - an observation vector z
  /// - a data model derived from LeastSquaresModelBase, which includes
  ///     - the model function,
  ///     - its Jacobian
  ///
  /// The cost function in L-M is always the inner product of the
  /// difference between an observation and the expected observation
  /// given the model parameters.  This means we can compute the cost
  /// function and its derivatives if we know the measurement function
  /// and its derivatives.
  ///

#define VW_MATH_LM_ABS_TOL (1e-16)
#define VW_MATH_LM_REL_TOL (1e-16)
#define VW_MATH_LM_MAX_ITER (100)

  template <class ImplT>
  typename ImplT::domain_type levenberg_marquardt( LeastSquaresModelBase<ImplT> const& least_squares_model,
                                                   typename ImplT::domain_type const& seed,
                                                   typename ImplT::result_type const& observation,
                                                   int &status,
                                                   double abs_tolerance = VW_MATH_LM_ABS_TOL,
                                                   double rel_tolerance = VW_MATH_LM_REL_TOL,
                                                   double max_iterations = VW_MATH_LM_MAX_ITER) {

    status = optimization::eDidNotConverge;

    const ImplT& model = least_squares_model.impl();
    bool done = false;
    double Rinv = 10;
    double lambda = 0.1;

    typename ImplT::domain_type x_try, x = seed;
    typename ImplT::result_type h = model(x);
    typename ImplT::result_type error = model.difference(observation, h);
    double norm_start = norm_2(error);

    VW_OUT(DebugMessage, "math") << "LM: initial guess for the model is " << seed << std::endl;
    VW_OUT(VerboseDebugMessage, "math") << "LM: starting error " << error << std::endl;
    VW_OUT(DebugMessage, "math") << "LM: starting norm is: " << norm_start << std::endl;

    // Solution may already be good enough
    if (norm_start < abs_tolerance) {
      status = optimization::eConvergedAbsTolerance;
      VW_OUT(DebugMessage, "math") << "CONVERGED TO ABSOLUTE TOLERANCE\n";
      done = true;
    }

    int outer_iter = 0;
    while (!done){

      bool shortCircuit = false;
      outer_iter++;
      VW_OUT(DebugMessage, "math") << "LM: outer iteration " << outer_iter << "   x = " << x << std::endl;

      // Compute the value, derivative, and hessian of the cost function
      // at the current point.  These remain valid until the parameter
      // vector changes.

      // expected measurement with new x
      h = model(x);

      // Difference between observed and predicted and error (2-norm of difference)
      error = model.difference(observation, h);
      norm_start = norm_2(error);
      //VW_OUT(DebugMessage, "math") << "LM: outer iteration starting robust norm: " << norm_start << std::endl;

      // Measurement Jacobian
      typename ImplT::jacobian_type J = model.jacobian(x);

      Vector<double> del_J = -1.0 * Rinv * (transpose(J) * error);

      // Hessian of cost function (using Gauss-Newton approximation)
      Matrix<double> hessian = Rinv * (transpose(J) * J);

      int iterations = 0;
      double norm_try = norm_start+1.0;
      while (norm_try > norm_start){

        // Increase diagonal elements to dynamically mix gradient
        // descent and Gauss-Newton.
        Matrix<double> hessian_lm = hessian;
        for ( unsigned i=0; i < hessian_lm.rows(); ++i ){
          hessian_lm(i,i) += hessian_lm(i,i)*lambda + lambda;
        }

        // Solve for update
        typename ImplT::domain_type delta_x;
        if (hessian_lm.rows() <= 2 && det(hessian_lm) > 0.0){
          // Direct method is more efficient for small matrices, also
          // here we avoid calling LAPACK which we've seen misbehave
          // in this situation in a multi-threaded environment.
          delta_x = inverse(hessian_lm)*del_J;
        }else{
          try{
            // By construction, hessian_lm is symmetric and
            // positive-definite.
            delta_x = solve_symmetric(hessian_lm, del_J);
          }catch ( const ArgumentErr& e ) {
            // If lambda is very small, the matrix becomes numerically
            // singular. In that case use the more general
            // least_squares solver.
            delta_x = least_squares(hessian_lm, del_J);
          }

        }

        // update parameter vector
        x_try = x - delta_x;

        typename ImplT::result_type h_try = model(x_try);

        typename ImplT::result_type error_try = model.difference(observation, h_try);
        norm_try = norm_2(error_try);

        //VW_OUT(VerboseDebugMessage, "math") << "LM: inner iteration " << iterations << " error is " << error_try << std::endl;
        //VW_OUT(DebugMessage, "math") << "\tLM: inner iteration " << iterations << " norm is " << norm_try << std::endl;

        if (norm_try > norm_start)
          // Increase lambda and try again
          lambda *= 10;

        ++iterations; // Sanity check on iterations in this loop
        if (iterations > 5) {
          //VW_OUT(DebugMessage, "math") << "\n****LM: too many inner iterations - short circuiting\n" << std::endl;
          shortCircuit = true;
          norm_try = norm_start;
        }
        //VW_OUT(DebugMessage, "math") << "\tlambda = " << lambda << std::endl;
      }

      // Percentage change convergence criterion
      if (((norm_start-norm_try)/norm_start) < rel_tolerance) {
        status = optimization::eConvergedRelTolerance;
        VW_OUT(DebugMessage, "math") << "CONVERGED TO RELATIVE TOLERANCE\n";
        done = true;
      }

      // Absolute error convergence criterion
      if (norm_try < abs_tolerance) {
        status = optimization::eConvergedAbsTolerance;
        VW_OUT(DebugMessage, "math") << "CONVERGED TO ABSOLUTE TOLERANCE\n";
        done = true;
      }

      // Max iterations convergence criterion
      if (outer_iter >= max_iterations) {
        VW_OUT(DebugMessage, "math") << "REACHED MAX ITERATIONS!";
        done = true;
      }

      // Take trial parameters as new parameters
      // If we short-circuited the inner loop, then we didn't actually find a
      // better p, so don't update it.
      if (!shortCircuit)
        x = x_try;

      // Take trial error as new error
      norm_start = norm_try;

      // Decrease lambda
      lambda /= 10;
      //VW_OUT(DebugMessage, "math") << "lambda = " << lambda << std::endl;
      //VW_OUT(DebugMessage, "math") << "LM: end of outer iteration " << outer_iter << " with error " << norm_try << std::endl;
    }
    VW_OUT(DebugMessage, "math") << "LM: finished with: " << outer_iter << "\n";
    return x;
  }

}} // namespace vw::math

#endif // __VW_OPTIMIZATION_H__
