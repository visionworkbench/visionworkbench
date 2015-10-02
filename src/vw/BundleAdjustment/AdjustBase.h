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


/// \file AdjustBase.h
///
/// This is the generic outline for code that performs bundle
/// adjustment. This provides a template for the interface with the
/// developer. Other bundleadjustment ideas can be built on top of
/// this, like reference impl, sparse impl, and a research impl. This also
/// contains all the code for cost functions.

#ifndef __VW_BUNDLEADJUSTMENT_ADJUST_BASE_H__
#define __VW_BUNDLEADJUSTMENT_ADJUST_BASE_H__

#include <vw/BundleAdjustment/ModelBase.h>
#include <boost/foreach.hpp>

namespace vw {
namespace ba {

  // ROBUST COST FUNCTION
  //----------------------------------------------------------------
  // Robust cost functions.  These cost function can help to reduce
  // the impact of outliers in the bundle adjustment.
  struct PseudoHuberError {
    double m_b;
    double m_b_2;
    PseudoHuberError(double b) : m_b(b) {
      m_b_2 = b*b;
    }

    double operator() (double delta_norm) {
      return 2.0f * m_b_2* (sqrt(1.0f + delta_norm*delta_norm/m_b_2) - 1.0f);
    }

    std::string name_tag () const { return "PsudeoHuberError"; }
    double threshold() const { return m_b; }
  };

  struct HuberError {
    double m_b;
    HuberError(double b) : m_b(b) {}

    double operator() (double delta_norm) {
      if (delta_norm < m_b)
        return delta_norm*delta_norm;
      else
        return 2*m_b*delta_norm - m_b*m_b;
    }

    std::string name_tag () const { return "HuberError"; }
    double threshold() const { return m_b; }
  };

  struct L1Error {
    double operator() (double delta_norm) { return fabs(delta_norm); }

    std::string name_tag () const { return "L1Error"; }
    double threshold() const { return 0.0; }
  };


  struct L2Error {
    double operator() (double delta_norm) { return delta_norm*delta_norm; }

    std::string name_tag () const { return "L2Error"; }
    double threshold() const { return 0.0; }
  };

  // Our Cauchy Error is missing a sigma^2 at the front of the log,
  // this is okay. Now when increasing sigma the weighting on higher error
  // is relaxed, increasing sigma increases the tail thickness of our PDF.
  struct CauchyError {
    double m_sigma;
    double m_sigma_2;
    CauchyError(double sigma) : m_sigma(sigma) {
      m_sigma_2 = m_sigma*m_sigma;
    }

    double operator() (double delta_norm) {
      return log(1+delta_norm*delta_norm/m_sigma_2);
    }

    std::string name_tag () const { return "CauchyError"; }
    double threshold() const { return m_sigma; }
  };


  // CHOLEKSY MATH FUNCTIONS
  //--------------------------------------------------------
  // DEVELOPER NOTE TO SELF:
  // Can we rely on LAPACK for these functions?

  // Returns cholesky factor L, D in lower left hand corner and
  // diagonal modifies in place returns 1 if original matrix was
  // positive definite, 0 otherwise has verbose output if not positive
  // definite
  template<class T>
  inline unsigned mod_cholesky(Matrix<T>& M){
    unsigned n = M.rows();
    for(unsigned i = 0; i < n; i++){
      for(unsigned j = i; j < n; j++){
        T sum = M(i,j);
        for(int k = i-1; k >= 0; k--)
          sum -= M(i,k) * M(j,k);

        if (i == j){
          if (sum <= 0.0){
            vw_out() << " Not positive definite! " << "\n";
            return 0;
          }
          M(i,i) = sqrt(sum);
        }else
          M(j,i) = sum/M(i,i);
      }
    }
    return 1;
  };

  // returns cholesky factor L, D in lower left hand corner and
  // diagonal modifies in place Assumes M is positive definite
  // DEVELOPER NOTE TO SELF:
  // This function is not used in any of the BA code!!!!!!!!!!!!!
  template<class T>
  inline void cholesky(Matrix<T>& M){
    unsigned n = M.rows();
    for (unsigned j = 0; j < n; j++){

      if(j > 0)
        submatrix(M, j, j, n-j, 1) -= submatrix(M, j, 0, n-j, j) * transpose(submatrix(M, j, 0, 1, j));

      submatrix(M, j, j, n-j, 1)/= sqrt(M(j,j));
    }
  };

  // Replaces lower triangle of A with that of L^{-1} Replaces upper
  // triangle of A with zeros Returns 0 if A is not positive definite
  template<class T>
  inline unsigned chol_inverse(Matrix<T>& A){
    unsigned n = A.rows();
    unsigned ret = mod_cholesky(A);
    if (ret == 0)
       return ret;
    T sum;
    for(unsigned i=0; i<n; i++){
      A(i,i) = 1/A(i,i);
      for(unsigned j = i+1; j < n; j++){
        sum = 0.0;
        for(unsigned k = i; k < j; k++)
          sum -= A(j,k) * A(k,i);
        A(j,i) = sum/A(j,j);
        A(i,j) = 0;
      }
    }
    return ret;
  };

  // solves Lx = b
  // modifies in place
  template<class T>
  inline void fsolve(Vector<T>& b, Matrix<T>& L){

    unsigned n = L.rows();

    b(0) /= L(0,0);
    for(unsigned i = 1; i < n; i++){
      Vector<double, 1> temp = submatrix(L, i, 0, 1, i) * subvector(b, 0, i);
      b(i) = (b(i) - temp(0))/L(i,i);
    }

  };

  // solves L^T y = x
  template<class T>
  inline void bsolve(Vector<T>& b, Matrix<T>& U){

    unsigned n = U.rows();

    b(n-1) /= U(n-1, n-1);
    for(unsigned i = b.size()-1; i>=1; i--){
      Vector<double, 1> temp = submatrix(U, i-1, i, 1, b.size()-i) * subvector(b, i, b.size()-i);
      b(i-1) = (b(i-1) - temp(0))/U(i-1,i-1);
    }

  };

  // solves Ax = b for positive definite A returns solution in b, and
  // cholesky decomp. of A in L returns 0 if A is not positive
  // definite
  template<class T>
  inline unsigned solve(Vector<T>&b, Matrix<T>& A){

    unsigned ret = mod_cholesky(A);

    if (ret == 0)
      return ret;

    fsolve(b, A);
    Matrix<T> U = transpose(A);
    bsolve(b, U);

    return ret;
  };

  // BUNDLE ADJUSTMENT BASE
  //--------------------------------------------------------
  // This is a base class for the item which actually performs the
  //bundle adjustment calculations.
  template <class BundleAdjustModelT, class RobustCostT>
  class AdjustBase {

  protected:
    boost::shared_ptr<ControlNetwork> m_control_net;
    BundleAdjustModelT & m_model;
    RobustCostT m_robust_cost_func;

    unsigned m_control; /* 0 = fletcher control, good for L2 cost
                           1 = old style / * by 10s            */
    double m_lambda;
    double m_nu;        // Used by fletcher control
    double g_tol;       // Surprisingly not used (8/22/09)
    double d_tol;       // Surprisingly not used (8/22/09)
    int m_iterations;

    bool m_use_camera_constraint;
    bool m_use_gcp_constraint;

  public:
    // Constructor
    AdjustBase( BundleAdjustModelT &model,
                RobustCostT const& robust_cost_func,
                bool use_camera_constraint=true,
                bool use_gcp_constraint=true ) :
    m_model(model), m_robust_cost_func(robust_cost_func),
      m_use_camera_constraint(use_camera_constraint),
      m_use_gcp_constraint(use_gcp_constraint) {

      m_iterations = 0;
      m_control_net = m_model.control_network();

      m_lambda = 1e-3;
      m_control = 0;
      m_nu = 2;
      g_tol = 1e-10;
      d_tol = 1e-10;

      // Useful declarations
      unsigned num_cam_params = BundleAdjustModelT::camera_params_n;
      unsigned num_pt_params = BundleAdjustModelT::point_params_n;
      unsigned num_cameras = m_model.num_cameras();
      unsigned num_ground_control_points = m_control_net->num_ground_control_points();
      unsigned num_observations = 2*m_model.num_pixel_observations();

      if (m_use_camera_constraint)
        num_observations += num_cameras*num_cam_params;
      if (m_use_gcp_constraint)
        num_observations += num_ground_control_points*num_pt_params;

      // Compute the initial error
      Vector<double> epsilon(num_observations);   // Error vector
      int idx = 0;
      int i = 0;
      for (ControlNetwork::const_iterator cpoint = m_control_net->begin();
           cpoint != m_control_net->end(); ++cpoint ) {
        for (ControlPoint::const_iterator cmeasure = cpoint->begin();
             cmeasure != cpoint->end(); ++cmeasure ) {
          int camera_idx = cmeasure->image_id();

          // Apply robust cost function weighting and populate error
          // vector
          Vector2 unweighted_error;
          try {
            unweighted_error = cmeasure->dominant() -
              m_model.cam_pixel(i, camera_idx,m_model.cam_params(camera_idx),
                      m_model.point_params(i));
          } catch (const camera::PointToPixelErr& e) {
            vw_out(WarningMessage,"ba") << "Unable to calculate starting error for point";
          }
          double mag = norm_2(unweighted_error);
          double weight = sqrt(m_robust_cost_func(mag))/mag;
          subvector(epsilon,2*idx,2) = unweighted_error * weight;

          ++idx;
        }
        i++;
      }

      // Add rows epsilon for a priori position/pose constraints ...
      if (m_use_camera_constraint)
        for (unsigned j=0; j < num_cameras; ++j)
          subvector(epsilon,
                    2*m_model.num_pixel_observations() + j*num_cam_params,
                    num_cam_params) = m_model.cam_target(j)-m_model.cam_params(j);

      // .. and the position of the 3D points to epsilon ...
      if (m_use_gcp_constraint) {
        idx = 2*m_model.num_pixel_observations();
        if ( m_use_camera_constraint )
          idx += num_cameras*num_cam_params;

        for (unsigned i = 0; i < m_model.num_points(); ++i )
          if ( (*m_control_net)[i].type() == ControlPoint::GroundControlPoint ) {
            subvector( epsilon, idx, num_pt_params) =
              m_model.point_target(i)-m_model.point_params(i);
            idx += num_pt_params;
          }
      }
    }

    // Access to inner templates
    typedef BundleAdjustModelT model_type;
    typedef RobustCostT cost_type;

    // Operational Controls
    double lambda() const { return m_lambda; }
    void set_lambda(double lambda) { m_lambda = lambda; }
    unsigned control() const { return m_control; }
    void set_control(unsigned control) { m_control = control; }
    bool camera_constraint() const { return m_use_camera_constraint; }
    bool gcp_constraint() const { return m_use_gcp_constraint; }

    // Additional Information
    int iterations() const { return m_iterations; }
    RobustCostT costfunction() const { return m_robust_cost_func; }
    BundleAdjustModelT& bundle_adjust_model() { return m_model; }

    // This is called repeatedly
    double update( double &abs_tol,
                   double &rel_tol ) {
      vw_throw( vw::NoImplErr() << "Programmer needs to implement AdjustBase::update()" );
      return 0;
    }

  };

}} // namespace vw::ba

#endif//__VW_BUNDLEADJUSTMENT_ADJUST_BASE_H__
