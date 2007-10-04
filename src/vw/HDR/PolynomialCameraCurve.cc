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

/// \file PolynomialCameraCurve.cc
/// 
/// Functions for deducing the camera response curve by comparing
/// relative brightness values from several images of the same scene.
/// The resulting curve is modelled as a polynomial.

#include <vw/HDR/PolynomialCameraCurve.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <boost/algorithm/string.hpp>

namespace vw {
namespace hdr {

  /* - - - - - - - - - - - - - - - - - - - - - - - - - -  
   *     InverseCameraCurveCostFunction implementation
   *
   * This is the function f_i(x) inside the nonlinear 
   * least squares problem:
   *
   * error = sum( [ f_i(x) ]^2 )  over all i
   * - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /// Cost function objects for the nonlinear least squares optimizer
  ///
  /// The VW least squares optimizer will optimize a function that 
  /// is a CRTP subclass of LeastSquaresModelBase.  Here we specialize
  /// for our problem of solving for the camera curve by implementing
  /// the methods operator() and jacobian().
  class InverseCameraCurveCostFn : public vw::math::LeastSquaresModelBase<InverseCameraCurveCostFn> {
    
    Matrix<double> m_pairs;
    
  public:
    typedef Vector<double> domain_type;
    typedef Vector<double> result_type;
    typedef Matrix<double> jacobian_type;
    
    InverseCameraCurveCostFn(Matrix<double> pixel_pairs): m_pairs(pixel_pairs) {}
    result_type operator()( domain_type const& x ) const {
    
      // Remove the gauge freedom from the problem by forcing f(1) = 1.
      domain_type theta = x/sum(x);
      
      double residual = 0;
      for (int32 i = 0; i < m_pairs.rows(); ++i) {
        double x_i = m_pairs(i,0);
        double x_j = m_pairs(i,1);
        double R_ij = m_pairs(i,2);
        
        // Compute the function F1 
        residual += psi(x_i,theta) - R_ij * psi(x_j,theta);
      }
      //      vw_out(0) << "Residual: " << residual << "  " << "\n";
      InverseCameraCurveCostFn::result_type return_val(1);
      return_val(0) = residual;
      return return_val;
    }
  };

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
   */
  vw::Vector<double> estimate_polynomial_camera_curve(vw::Matrix<double> &pixel_pairs,
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
    InverseCameraCurveCostFn cost_fn(pixel_pairs);
    Vector<double> target(1);
    fill(target,0.0);
    int status;

    Vector<double> result = vw::math::levenberg_marquardt(cost_fn, theta, target, status, 0.1, 1e-16, 1000);
    if (status == vw::math::optimization::eDidNotConverge) 
      vw_throw( LogicErr() << "Levenberg Marquardt did not converge when computing the camera curve." );

    /* 
     * Copy the result, but normalize it first so that f(1) = 1
     * (i.e. the sum of the polynomial coefficients sums to 1.  This
     * eliminated the gauge freedom in the camera curve estimation
     * problem.
     */
    return result / sum(result);
  }


  class InverseCostFn : public vw::math::LeastSquaresModelBase<InverseCostFn> {
    Matrix<double> m_pairs;
  public:
    typedef Vector<double> domain_type;
    typedef Vector<double> result_type;
    typedef Matrix<double> jacobian_type;

    InverseCostFn(Matrix<double> pairs) : m_pairs(pairs) {}

    result_type operator()( domain_type const& x ) const {

      // Remove the gauge freedom from the problem by forcing f(1) = 1.
      domain_type theta = x/sum(x);

      double residual = 0;
      for (int32 i = 0; i < m_pairs.rows(); ++i) {
        double x_i = m_pairs(i,0);
        double L_i = m_pairs(i,1);
        // Compute the function F1 
        residual += psi(L_i,theta) - x_i;
      }
      //      vw_out(0) << "Inverse Residual: " << residual << "  " << "\n";
      InverseCostFn::result_type return_val(1);
      return_val(0) = residual;
      return return_val;
    }
  };

  std::vector<vw::Vector<double> > invert_polynomial(std::vector<Vector<double> > &curves,
                                                     double lower_bound, double upper_bound) {

    // Generate curves with the same polynomial order as the source curves.
    std::vector<Vector<double> > inverse_curves(curves.size());
    for (int c = 0; c < curves.size(); ++c) {

      // Generate (scaled radiance, LDR intensity) pairs to fit a curve to 
      const int samples = 3000;
      Matrix<double> pairs(samples, 2);
      double step = (upper_bound - lower_bound) / (samples-1);
      double val = lower_bound;
      for (int i = 0; i < samples; ++i) {
        pairs(i, 0) = val;
        pairs(i, 1) = psi(val, curves[c]);
        val += step;
      }

      Vector<double> theta(curves[0].size());
      fill(theta,1.0);
      
      // Set up the optimizer and optimizer options 
      InverseCostFn cost_fn(pairs);
      Vector<double> target(1);
      fill(target,0.0);
      int status;

      Vector<double> result = vw::math::levenberg_marquardt(cost_fn, theta, target, status, 0.1, 1e-16, 1000);
      if (status == vw::math::optimization::eDidNotConverge) 
        vw_throw( LogicErr() << "Levenberg Marquardt did not converge when inverting the camera curve polynomial." );

      inverse_curves[c] = result / sum(result);
    }
    return inverse_curves;
  }



  //--------------------------- ----------------- ------------------------------
  //--------------------------- File IO Functions ------------------------------
  //--------------------------- ----------------- ------------------------------

  void write_polynomial_curve_coefficients(std::string const& generated_curves_file, 
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

  void read_polynomial_curve_coefficients(std::vector<vw::Vector<double> > &curves, 
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


  /// Read and write the actual camera curves themselves in a
  /// tabulated data format.
  void write_polynomial_curves(std::string const& curves_file, 
                    std::vector<vw::Vector<double> > const &curves) {
    FILE* output_file = fopen(curves_file.c_str(), "w");
    if ( !output_file ) vw_throw( IOErr() << "write_curves: failed to open file for writing." );
    for (double i = 0; i < 1.0; i+=0.01) {
      for ( unsigned j = 0; j < curves.size(); ++j ) {
        fprintf(output_file, "%f ", psi(i,curves[j]));
      }
      fprintf(output_file, "\n");
    }
    fclose(output_file);
  }

}} // namespace vw::HDR

