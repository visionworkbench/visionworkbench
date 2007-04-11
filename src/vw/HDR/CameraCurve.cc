// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
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

/// \file CameraCurve.cc
/// 
/// Functions for deducing the camera response curve by comparing 
/// relative brightness values from several images of the same scene.

#include <vw/HDR/CameraCurve.h>
#include <boost/algorithm/string.hpp>

namespace vw {
namespace hdr {

  void write_curves_file(std::string const& generated_curves_file, 
                         std::vector<vw::Vector<double> > const &curves) {
    
    FILE* output_file = fopen(generated_curves_file.c_str(), "w");
    if ( !output_file ) vw_throw( IOErr() << "write_curves_file: failed to open file for writing." );
    for ( unsigned i = 0; i < curves.size(); ++i ) {
      for ( unsigned j = 0; j < curves[i].size(); ++j ) {
        fprintf(output_file, "%f ", curves[i][j]);
      }
      fprintf(output_file, "\n");
    }    
    fclose(output_file);
  }


  void read_curves_file(std::vector<vw::Vector<double> > &curves, 
                        std::string const& curves_input_file) {
    
    FILE* input_file = fopen(curves_input_file.c_str(), "r");
    if ( !input_file ) vw_throw( IOErr() << "read_curves_file: failed to open file for reading." );
    
    char c_line[2048];
    
    curves.clear();
    while ( !feof(input_file) ) {
      if ( !fgets(c_line, 2048, input_file) )
        break;
      std::string line = c_line;
      boost::trim_left(line); 
      boost::trim_right(line);

      std::vector< std::string > split_vec; // #2: Search for individual values
      boost::split( split_vec, line, boost::is_any_of(" ") );
      Vector<double> curve(split_vec.size());
      for ( unsigned i = 0; i < split_vec.size(); ++i ) {
        curve[i] = atof(split_vec[i].c_str());
      }
      curves.push_back(curve);
    }
    fclose(input_file);
  }


  /*
   * Approximate estimation of the camera response curve as a 
   * polynomial.  This appoarch uses a random sampling of pixel 
   * values from mutliple exposures of the same image.  The 
   * solution can be computed in closed form, but we have actually 
   * found it to be more stable (better conditioned) to compute it 
   * using a nonlinear least squares iterative solver.  
   * 
   * The input to this function, pairs, is an Nx3 matrix where each
   * row contains a pixel value from image 1, the corresponding pixel 
   * value from image 2, and the ratio of exposure between these two
   * images.
   * 
   * So far, we have found that 1 photographic stop (a factor of 2 
   * in shutter speed or 1 f-stop) equates to a difference of sqrt(2) 
   * in the amount of light that falls on the sensor.  Therefore, for 
   * photos that are seperated by 1 f-stop, the ratios in the third 
   * column of pixel_pairs should be powers of sqrt(2).
   */
  void estimate_inverse_camera_curve(Matrix<double> &pixel_pairs,
                                     Vector<double> &poly,
                                     unsigned int polynomial_order) {
  
    /* 
     * Initialize the polynomial coefficients to 1.0
     * 
     * Note: We force the coefficient of the largest power of x to be 
     * 1-sum(other coefficients) to remove the gauge freedom in this 
     * problem, so we do not explicitly optimize over this coefficient 
     * in this optimization.
     */
    Vector<double> theta(polynomial_order);
    fill(theta,1.0);
  
    /* Set up the optimizer and optimizer options */
    CameraCurveCostFn cost_fn(pixel_pairs);
    Vector<double> target(pixel_pairs.rows());
    fill(target,0.0);
    int status;

    Vector<double> result = vw::math::levenberg_marquardt(cost_fn, theta, target, status, 0.1, 1e-16, 1000);
    if (status == vw::math::optimization::eDidNotConverge) 
      vw_throw( LogicErr() << "Levenberg Marquardt did not converge when computing the camera curve." );

    /* Copy the result into poly 
     * 
     * The coefficient associated with x^N is equal to 1.0, which 
     * is a restriction we have placed on the problem by forcing 
     * f(1) = 1.0.  This eliminated the gauge freedom in the camera
     * curve estimation problem.
     * 
     * The polynomial coefficients are from x^0 to x^N.
     */ 
    Vector<double> output_coeff(result.size() + 1);
    for ( unsigned i = 0; i < result.size(); ++i ) {
      output_coeff(i) = result(i);
    }
    output_coeff(polynomial_order) = 1 - sum(result);
    poly = output_coeff;
  }


  /* - - - - - - - - - - - - - - - - - - - - - - - - - -  
   *     CameraCurveCostFunction implementation
   *
   * This is the function f_i(x) inside the nonlinear 
   * least squares problem:
   *
   * error = sum( [ f_i(x) ]^2 )  over all i
   * - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* 
   * Evaluate the polynomial given the set of coefficients 
   * provided in theta as the first N-1 coefficients.  The
   * last coefficient in the polynomial is forced to be equal
   * to 1-sum(theta), which forces the polynomial to meet the
   * restriction that f(1) = 1.  This eliminates the gauge 
   * freedom in the camera response curve estimation problem.
   */
  inline double psi(double x, Vector<double> const &theta) {    
    unsigned int N = theta.size();
  
    // Compute the main term:    (sum(k = 0, N-1, theta_k * x^k)) 
    double L = 0;
    for (unsigned int i = 0; i < theta.size(); i++) {
      L += theta(i) * pow(x, i); 
    }
  
    // Add the constraint term:  (x^N * (1 - sum(k = 0, N-1, theta_k))) 
    double theta_sum = sum(theta);
    L += pow(x,N)*(1 - theta_sum);
    return L;
  }

  // Compute the value of the function f at x
  CameraCurveCostFn::result_type CameraCurveCostFn::operator()( domain_type const& x ) const {
  
    Vector<double> fx(m_pairs.rows());
    for (int32 i = 0; i < m_pairs.rows(); ++i) {
      double x_i = m_pairs(i,0);
      double x_j = m_pairs(i,1);
      double R_ij =  m_pairs(i,2);
      
      // Compute the function F1 
      fx(i) = psi(x_i,x) - R_ij * psi(x_j,x);
    }
    return fx;
  }

  // Compute the Jacobian of f at x 
  Matrix<double> CameraCurveCostFn::jacobian( domain_type const& x ) const {

    Matrix<double> J(m_pairs.rows(), x.size());
    unsigned int N = x.size();

    for (int32 i = 0; i < m_pairs.rows(); i++) {
      double x_i = m_pairs(i,0);
      double x_j = m_pairs(i,1);
      double R_ij = m_pairs(i,2);
    
      // Compute the Jacobian of F: J 
      for (unsigned int nn = 0; nn < N; nn++) 
        J(i,nn) = (pow(x_i,nn) - pow(x_i,N)) - R_ij * (pow(x_j,nn) - pow(x_j,N));

    }
    return J;
  }

  void invert_curve(Vector<double> &curve, Vector<double> &inverse,
                    unsigned int polynomial_order, double lower_bound,
                    double upper_bound) {
    /* 
     * Initialize the polynomial coefficients to 1.0
     * 
     * Note: We force the coefficient of the largest power of x to be 
     * 1-sum(other coefficients) to remove the gauge freedom in this 
     * problem, so we do not explicitly optimize over this coefficient 
     * in this optimization.
     */
    Vector<double> theta(polynomial_order);
    fill(theta,1.0);

    // Generate (scaled radiance, LDR intensity) pairs to fit a curve to 
    const int RESIDUALS = 1000;
    Matrix<double> pairs(RESIDUALS, 2);
    double val = lower_bound;
    double step = (upper_bound - lower_bound) / (RESIDUALS - 1.0);
    for (int i = 0; i < RESIDUALS; i++) {
      pairs(i, 1) = val;
      pairs(i, 0) = psi(val, curve);
      val += step;
    }

    // Set up the optimizer and optimizer options 
    InverseCostFn cost_fn(pairs);
    Vector<double> target(pairs.rows());
    fill(target,0.0);
    int status;

    Vector<double> result = vw::math::levenberg_marquardt(cost_fn, theta, target, status, 0.1, 1e-16, 1000);
    if (status == vw::math::optimization::eDidNotConverge) 
      vw_throw( LogicErr() << "Levenberg Marquardt did not converge when computing the inverse camera curve." );
    
    /* Copy the result into inverse
     * 
     * The coefficient associated with x^N is equal to 1.0, which 
     * is a restriction we have placed on the problem by forcing 
     * f(1) = 1.0.  This eliminated the gauge freedom in the camera
     * curve estimation problem.
     * 
     * The polynomial coefficients are from x^0 to x^N.
     */
    Vector<double> output_coeff(result.size() + 1);
    for ( unsigned i = 0; i < result.size(); ++i ) {
      output_coeff(i) = result(i);
    }
    output_coeff(polynomial_order) = 1 - sum(result);
    inverse = output_coeff;
  }

  // Compute the value of the function f at x 
  InverseCostFn::result_type InverseCostFn::operator()( domain_type const& x ) const {
    Vector<double> fx(m_pairs.rows());
    for ( int32 i = 0; i < m_pairs.rows(); ++i ) {
      double x_i = m_pairs(i,0);
      double x_j = m_pairs(i,1);
    
      // Compute the function F1 
      fx(i) = psi(x_i,x) - x_j;
    }
    return fx;
  }

}} // namespace vw::HDR
