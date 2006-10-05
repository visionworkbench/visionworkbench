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

/// \file Transform.h
/// 
/// Optimization classes for carrying out nonlinear optimization.
/// 
/// This code provides some optimization routines for doing estimation.
/// The goal is to provide some sufficiently general support to that
/// anything from Levenberg-Marquardt for ICP to Kalman filters for
/// navigation can use the same toolbox.

#ifndef __VW_OPTIMIZATION_H__
#define __VW_IMAGE_TRANSFORM_H__

// Vision Workbench
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/LinearAlgebra.h>

// Boost 
#include <boost/concept_check.hpp>

#define ABS_TOL (0.001)
#define REL_TOL (0.001)

namespace vw {

  // Optimization

  // Levenberg-Marquardt is an algorithm for solving problems of the
  // form
  //
  // J(p) = sum_i ( z_i - h(x_i;p) )^2
  //
  // That is, a least squares problem where the objective is to find a
  // parameter vector p such that the model function h(x;p) for given
  // data x, evaluates as closely as possible to the observations z in
  // a least squares sense.

  // First thing we need is a generic idea of a measurement function
  // or model function.  The model function needs to provide a way to
  // evaluate h(x;p) as well as a way to differentiate h() to get its
  // Jacobian.  In the past we've implemented a fallback option by
  // computing the Jacobian numerically which seems to work pretty
  // well.  This ModelFunction might be implemented as an inherited
  // class, or something like the transform concepts.  For now I'm
  // just going to focus on providing a quick port of
  // Levenberg-Marquardt to get started so I'll make up a bogus model
  // function and use it directly.
  template <class T>
  class ModelFunction {
  private:
    Vector<T> m_x; // data
  public:
    ModelFunction( const Vector<T>& x )
    {
      m_x = x;
    }

    inline Vector<T> operator()( const Vector<T>& p ) const {
      // Evaluate h(x;p) provided the parameter p

      // For now make up a function to get started with
      Vector<T> h(5);
      h(0) = sin(p(0)+0.1);
      h(1) = cos(p(1) * p(2));
      h(2) = p(1) * cos(p(2));
      h(3) = atan2(p(0),p(3));
      h(4) = atan2(p(2),p(1));
      return h;
    }
    
    // Less obvious is the utility of a differencing function, so that
    // things like angles can be handled without putting the logic in
    // the L-M implementation
    inline Vector<T> diff( const Vector<T>& a,
			   const Vector<T>& b ) const
    {
      return (a-b);
    }
      
    inline Matrix<T> Jacobian( const Vector<T>& p ) const {
      // Evaluate dh(x;p)/dp provided the parameter p.  One option is
      // to put a numerical derivative in the base class and use it
      // unless overloaded.
      
      // Get nominal function value
      Vector<T> h0 = this->operator()(p);
      //std::cout << "h0 = " << h0 << std::endl;
      
      // Jacobian is #params x #outputs
      Matrix<T> H(h0.size(), p.size());
      
      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for (int i=0; i<p.size(); i++){
	// Take a delta step in each dimension
	Vector<T> pi = p;
	// Variable step size, depending on parameter value
	double epsilon = 1e-7 + pi(i)*1e-7;
	pi(i) += epsilon;
	// Evaluate function with this step
	Vector<T> hi = this->operator()(pi);
	//std::cout << "hi = " << hi << std::endl;
	//std::cout << "H = " << H << std::endl;
	// Compute derivative w.r.t. parameter i
	select_col(H,i) = (hi-h0)/epsilon;
      }
      return H;
    }
  };

  // Levenberg-Marquardt optimization
  template <class T>
  void levenberg_marquardt( Vector<T> &p,
			    Vector<T> &z,
			    const Vector<T> &x, // do we need this?
			    const ModelFunction<T>& model )
  {
    // Implements Levenberg Marquardt optimization.  
    // Requires:
    // - an initial parameter vector p
    // - an observation z
    // - some data x
    // - a sensor model, which includes 
    //     - the sensor model function, 
    //     - its Jacobian, and 
    //     - the sensor noise covariance.
    
    // The cost function in L-M is always the inner product of the
    // difference between an observation and the expected observation
    // given the model parameters.  This means we can compute the cost
    // function and its derivatives if we know the measurement function
    // and its derivatives.
    
    //  double errorold;
    // double error, errortry;
    bool converged = false;
    
    double Rinv = 10; // inverse(model.R(p));
    
    Matrix<double> A, Alm;
    Vector<double> e, h, ptry, htry;
    // Vector<double> weight, norm;
    //Huber_robust_estimator estimator;
    //L2_robust_estimator estimator;
    Vector<double> delta_p;
    double lambda = 0.1;
    
    h = model(p);
    e = model.diff(z,h);
    
    //Matrix<double> em(e.get_size(), 1, e.begin());
    //  error = (em.transpose_times(Rinv)*em)(0, 0);
    //  error = Rinv * ((em.transpose_times(em)))(0, 0);
    // error = Rinv * e.dot(e);
    
    std::cout << "LM: initial guess for p is " << p << std::endl;
    // std::cout << "LM: starting error " << error << std::endl;
    
    // Compute robust weights
    // compute_robust_weight( Rinv*e, weight, estimator );
    // compute_robust_norm( Rinv*e, norm, estimator );
    
    double norm_start = 0.0;
    norm_start = sum( elem_prod(e,e) );
    //for (int i=0; i<e.rows(); i++){
    //  norm_start += norm(i);
    //}
    std::cout << "LM: starting norm is: " << norm_start << std::endl;
    
    // Solution may already be good enough
    //if (error<ABS_TOL)
    if (norm_start<ABS_TOL)
      converged = true;
    
    int outer_iter = 0;
    
    while (!converged){
      
      bool shortCircuit = false;
      
      std::cout << "LM: start of outer iteration " << ++outer_iter << std::endl;
      std::cout << "LM: p at start of outer iteration is " << p << std::endl;
      
      // Compute the value, derivative, and hessian of the cost function
      // at the current point.  These remain valid until the parameter
      // vector changes.
      
      // expected measurement with new p
      h = model(p);
      
      // Difference between observed and predicted
      e = model.diff(z, h);
      //    std::cout << "z - h = " << e << std::endl;
      
      // Shouldn't have to compute this, but let's debug it anyway
      // error = Rinv * e.dot(e);
      // std::cout << "LM: error is " << error << std::endl;
      
      // For robust estimation:
      // Use robust norm rather than least squares
      // Compute robust weights
      // compute_robust_weight( Rinv*e, weight, estimator );
      // compute_robust_norm( Rinv*e, norm, estimator );
      norm_start = 0.0;
      norm_start = sum( elem_prod( e,e ) );
      // For robust estimation:
      //for (int i=0; i<norm.get_num_of_cols(); i++){
      //norm_start += norm(i);
      //}
      std::cout << "LM: outer iteration starting robust norm: " << norm_start << std::endl;

      // For robust estimation:
      //Matrix<double> ee(e.get_size(), 1, e.begin());
      
      // Measurement Jacobian
      //std::cout << "computing jacobian" << std::endl;
      Matrix<double> H = model.Jacobian(p);
      //Matrix<double> hterinv = transpose(-1.0*H);
      
      //    std::cout << "-H' --> " << std::endl << hterinv << std::endl;
      
      // For robust estimation:
      // Multiply each ee(i) by robust weight
      // for (int i=0; i<ee.get_num_of_rows(); i++){
      // ee(i,0) *= weight(i);
      // }
      
      // For robust estimation:
      //hterinv = hterinv * ee;
      //hterinv = hterinv * e;
      //hterinv *= Rinv;
      //    std::cout << hterinv << std::endl;
      Vector<double> del_J;
      del_J = -1.0 * transpose(H) * Rinv * e;
      //    del_J = (H * -1.0).transpose_times(ee) * Rinv;
      std::cout << "LM: del_J is " << std::endl << del_J << std::endl;
      
      // Hessian of cost function (using Gauss-Newton approximation)
      A = transpose(H) * Rinv * H;
      
      // For robust estimation:
      // Hessian of robust cost function
      //Matrix<double> Hweighted;
      //Hweighted.resize(H.get_num_of_rows(), H.get_num_of_cols());
      //for (int i=0; i<H.get_num_of_rows(); i++){
      //for (int j=0; j<H.get_num_of_cols(); j++){
      //  Hweighted(i,j) = H(i,j) * sqrt(weight(i));
      //}
      //}
      //A = Hweighted.transpose_times(Hweighted) * Rinv;
      //std::cout << " LM: A is " << std::endl << A << std::endl;
      
      int iterations = 0;
      // errortry = error+1.0;
      // while (errortry>error){
      double norm_try = norm_start+1.0;
      while (norm_try>norm_start){
	
	// Increase diagonal elements to dynamically mix gradient
	// descent and Gauss-Newton.
	Alm = A;
	for (int i=0; i<Alm.rows(); i++){
	  Alm(i,i) += Alm(i,i)*lambda + lambda;
	}
	
	// Solve for update
	Matrix<double> Ainv = pseudoinverse( Alm );
	delta_p = Ainv * del_J;
	
	// update parameter vector
	std::cout << " *** p is " << p << std::endl
	     << " *** delta_p is " << delta_p << std::endl;
	ptry = p - delta_p;
	std::cout << "p is " << p << std::endl;
	std::cout << "delta p is " << delta_p << std::endl;
	std::cout << "p try is " << std::endl << ptry << std::endl;
	
	htry = model(ptry);
	
	Vector<double> difftry = model.diff(z, htry);
	//      Matrix<double> etmp(difftry.get_size(), 1, difftry.begin());
	//      errortry = (etmp.transpose_times(etmp))(0, 0) * Rinv;
	
	// Use robust norm here too:
	// errortry = difftry.dot(difftry) * Rinv;
	// std::cout << "LM: inner iteration " << iterations << " error is " 
	//   << errortry << std::endl;
	
	// For robust estimation:
	// compute_robust_norm( Rinv*difftry, norm, estimator );
	// Compute sum of norms for this trial solution
	// norm_try = 0.0;
	//for (int i=0; i<norm.get_num_of_cols(); i++){
	//  norm_try += norm(i);
	//}
	norm_try = sum( elem_prod( difftry, difftry ) );
	
	//std::cout << "LM: inner iteration " << iterations << " error is " 
	//   << errortry << std::endl;
	std::cout << "LM: inner iteration " << iterations << " norm is " 
	     << norm_try << std::endl;
	
	if (norm_try>norm_start)
	  // Increase lambda and try again
	  lambda *= 10;
	
	++iterations; // Sanity check on iterations in this loop
	if (iterations > 5) {
	  std::cerr << "LM: too many inner iterations - short circuiting" << std::endl;
	  shortCircuit = true;
	  norm_try = norm_start;
	}
	std::cout << "lambda = " << lambda << std::endl;
      }
      /* 
      // Percentage change convergence criterion
      if (((error-errortry)/errortry)<REL_TOL)
      converged = true;
      // Absolute error convergence criterion
      if (errortry<ABS_TOL)
      converged = true;
      // Take trial parameters as new parameters
      p = ptry;
      // Take trial error as new error
      error = errortry;
      // Decrease lambda
      lambda /= 10;
      std::cout << "LM: end of outer iteration " << outer_iter << " with error "
      << error << std::endl;
      */
      
      // Percentage change convergence criterion
      if (((norm_start-norm_try)/norm_start) < REL_TOL)
	converged = true;
      // Absolute error convergence criterion
      if (norm_try < ABS_TOL)
	converged = true;
      
      // Take trial parameters as new parameters
      // If we short-circuited the inner loop, then we didn't actually find a
      // better p, so don't update it.
      if (!shortCircuit)
	p = ptry;
      // Take trial error as new error
      norm_start = norm_try;
      // Decrease lambda
      lambda /= 10;
      std::cout << "lambda = " << lambda << std::endl;
      //std::cout << "LM: end of outer iteration " << outer_iter << " with error "
      //	 << error << std::endl;
      std::cout << "LM: end of outer iteration " << outer_iter << " with error "
	   << norm_try << std::endl;
    }
  }

} // namespace vw

#endif // __VW_OPTIMIZATION_H__
