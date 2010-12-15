// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file AdjustRef.h
///
/// Reference implementation of bundle adjustment. Very slow!

#ifndef __VW_BUNDLEADJUSTMENT_ADJUST_REF_H__
#define __VW_BUNDLEADJUSTMENT_ADJUST_REF_H__

#include <vw/BundleAdjustment/AdjustBase.h>
#include <vw/Math/LinearAlgebra.h>

// Boost
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/version.hpp>

// The sparse vectors/matrices are needed
// for the covariance calculation

#if BOOST_VERSION<=103200
// Mapped matrix doesn't exist in 1.32, but Sparse Matrix does
//
// Unfortunately some other tests say this doesn't work
#define boost_sparse_matrix boost::numeric::ublas::sparse_matrix
#define boost_sparse_vector boost::numeric::ublas::sparse_vector
#else
// Sparse Matrix was renamed Mapped Matrix in later editions
#define boost_sparse_matrix boost::numeric::ublas::mapped_matrix
#define boost_sparse_vector boost::numeric::ublas::mapped_vector
#endif

namespace vw {
namespace ba {

  template <class BundleAdjustModelT, class RobustCostT>
  class AdjustRef : public AdjustBase<BundleAdjustModelT,RobustCostT>, private boost::noncopyable {

    // Need to save S for covariance calculations
    math::Matrix<double> m_S;

  public:

    AdjustRef( BundleAdjustModelT & model,
               RobustCostT const& robust_cost_func,
               bool use_camera_constraint=true,
               bool use_gcp_constraint=true ) :
    AdjustBase<BundleAdjustModelT,RobustCostT>( model, robust_cost_func,
                                                use_camera_constraint,
                                                use_gcp_constraint ) {}

    Matrix<double> S() { return m_S; }
    void set_S(const math::Matrix<double>& S) {
      m_S = S;
    }

    // Covariance Calculator
    // __________________________________________________
    // This routine inverts a sparse matrix S, and prints the individual
    // covariance matrices for each camera
    void covCalc(){

      // camera params
      unsigned num_cam_params = BundleAdjustModelT::camera_params_n;
      unsigned num_cameras = this->m_model.num_cameras();

      unsigned inverse_size = num_cam_params * num_cameras;

      typedef Matrix<double, BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n> matrix_camera_camera;

      // final vector of camera covariance matrices
      vw::Vector< matrix_camera_camera > sparse_cov(num_cameras);

      // Get the S matrix from the model
      Matrix<double> S = this->S();
      Matrix<double> Id(inverse_size, inverse_size);
      Id.set_identity();
      Matrix<double> Cov = multi_solve_symmetric(S, Id);

      //pick out covariances of individual cameras
      for ( unsigned i = 0; i < num_cameras; i++ )
         sparse_cov(i) = submatrix(Cov, i*num_cam_params,
                                   i*num_cam_params,
                                   num_cam_params,
                                   num_cam_params);

      std::cout << "Covariance matrices for cameras are:"
                << sparse_cov << "\n\n";

      return;
    }

    // UPDATE IMPLEMENTATION
    //---------------------------------------------------------------
    // This is a simple, non-sparse, unoptimized implementation of LM
    // bundle adjustment.  It is primarily used for validation and
    // debugging.
    //
    // Each entry in the outer vector corresponds to a distinct 3D
    // point.  The inner vector contains a list of image IDs and
    // pixel coordinates where that point was imaged.
    double update(double &abs_tol, double &rel_tol) {
      ++this->m_iterations;

      // Here are some useful variable declarations that make the code
      // below more readable.
      unsigned num_cam_params = BundleAdjustModelT::camera_params_n;
      unsigned num_pt_params = BundleAdjustModelT::point_params_n;
      unsigned num_points = this->m_model.num_points();
      unsigned num_model_parameters = this->m_model.num_cameras()*num_cam_params +
        this->m_model.num_points()*num_pt_params;

      unsigned num_cameras = this->m_model.num_cameras();
      unsigned num_ground_control_points = this->m_control_net->num_ground_control_points();
      unsigned num_observations = 2*this->m_model.num_pixel_observations();
      if (this->m_use_camera_constraint)
        num_observations += num_cameras*num_cam_params;
      if (this->m_use_gcp_constraint)
        num_observations += num_ground_control_points*num_pt_params;

      // The core LM matrices and vectors
      Matrix<double> J(num_observations, num_model_parameters);   // Jacobian Matrix
      Vector<double> error(num_observations);                   // Error vector
      Matrix<double> sigma(num_observations, num_observations);   // Sigma (uncertainty) matrix

      // --- SETUP STEP ----
      // Add rows to J and error for the imaged pixel observations
      int idx = 0;
      for (unsigned i = 0; i < this->m_control_net->size(); ++i) {       // Iterate over control points
        for (unsigned m = 0; m < (*(this->m_control_net))[i].size();
             ++m) {  // Iterate over control measures
          int camera_idx = (*(this->m_control_net))[i][m].image_id();

          Matrix<double> J_a = this->m_model.A_jacobian(i,camera_idx,
                                                        this->m_model.A_parameters(camera_idx),
                                                        this->m_model.B_parameters(i));
          Matrix<double> J_b = this->m_model.B_jacobian(i,camera_idx,
                                                        this->m_model.A_parameters(camera_idx),
                                                        this->m_model.B_parameters(i));

          // Populate the Jacobian Matrix
          submatrix(J, 2*idx, num_cam_params*camera_idx, 2, num_cam_params) = J_a;
          submatrix(J, 2*idx, num_cam_params*num_cameras + i*num_pt_params, 2, num_pt_params) = J_b;

          // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error;
          try {
            unweighted_error = (*(this->m_control_net))[i][m].dominant() -
              this->m_model(i, camera_idx,
                            this->m_model.A_parameters(camera_idx),
                            this->m_model.B_parameters(i));
          } catch ( camera::PixelToRayErr &e ) {}
          double mag = norm_2(unweighted_error);
          double weight = sqrt(this->m_robust_cost_func(mag)) / mag;
          subvector(error,2*idx,2) = unweighted_error * weight;

          // Fill in the entries of the sigma matrix with the uncertainty of the observations.
          Matrix2x2 inverse_cov;
          Vector2 pixel_sigma = (*(this->m_control_net))[i][m].sigma();

          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));
          submatrix(sigma, 2*idx, 2*idx, 2, 2) = inverse_cov;

          ++idx;
        }
      }

      double max = 0.0;
      if (this->m_iterations == 1 && this->m_lambda == 1e-3){
        Matrix<double> hessian = transpose(J) * sigma * J;
        for (unsigned i = 0; i < hessian.rows(); ++i)
          if (fabs(hessian(i,i)) > max)
            max = fabs(hessian(i,i));
        this->m_lambda = 1e-10 * max;
      }

      // Add rows to J and error for a priori camera parameters...
      if ( this->m_use_camera_constraint )
        for (unsigned j=0; j < num_cameras; ++j) {
          Matrix<double> id(num_cam_params, num_cam_params);
          id.set_identity();
          submatrix(J,
                    2*this->m_model.num_pixel_observations() + j*num_cam_params,
                    j*num_cam_params,
                    num_cam_params,
                    num_cam_params) = id;
          Vector<double> unweighted_error = this->m_model.A_target(j) -
            this->m_model.A_parameters(j);
          subvector(error,
                    2*this->m_model.num_pixel_observations() + j*num_cam_params,
                    num_cam_params) = unweighted_error;
          submatrix(sigma,
                    2*this->m_model.num_pixel_observations() + j*num_cam_params,
                    2*this->m_model.num_pixel_observations() + j*num_cam_params,
                    num_cam_params, num_cam_params) = this->m_model.A_inverse_covariance(j);
        }

      // ... and the position of the 3D points to J and error ...
      if ( this->m_use_gcp_constraint ) {
        idx = 0;
        for (unsigned i=0; i < this->m_model.num_points(); ++i) {
          if ((*this->m_control_net)[i].type() == ControlPoint::GroundControlPoint) {
            Matrix<double> id(num_pt_params,num_pt_params);
            id.set_identity();
            submatrix(J,2*this->m_model.num_pixel_observations() +
                      num_cameras*num_cam_params + idx*num_pt_params,
                      num_cameras*num_cam_params + idx*num_pt_params,
                      num_pt_params,
                      num_pt_params) = id;
            Vector<double> unweighted_error = this->m_model.B_target(i)-this->m_model.B_parameters(i);
            subvector(error,
                      2*this->m_model.num_pixel_observations() +
                      num_cameras*num_cam_params + idx*num_pt_params,
                      num_pt_params) = unweighted_error;
            submatrix(sigma,
                      2*this->m_model.num_pixel_observations() +
                      num_cameras*num_cam_params + idx*num_pt_params,
                      2*this->m_model.num_pixel_observations() +
                      num_cameras*num_cam_params + idx*num_pt_params,
                      num_pt_params, num_pt_params) = this->m_model.B_inverse_covariance(i);
            ++idx;
          }
        }
      }


      // --- SOLVE UPDATE STEP -----------------------------------
      Matrix<double> JTS = transpose(J) * sigma;

      // Build up the right side of the normal equation...
      Vector<double> epsilon = -1.0 * (JTS * error);

      // ... and the left side.  (Remembering to rescale the diagonal
      // entries of the approximated hessian by lambda)
      Matrix<double> hessian = JTS * J;

      // initialize m_lambda on first iteration, ignore if user has
      // changed it.
      unsigned num_cam_entries = num_cam_params * num_cameras;
      unsigned num_pt_entries = num_pt_params * num_points;

      //WARNING: debugging only (uncomment two lines below)
      for ( unsigned i=0; i < hessian.rows(); ++i )
         hessian(i,i) +=  this->m_lambda;

      //Cholesky decomposition. Returns Cholesky matrix in lower left hand corner.
      Vector<double> delta = epsilon;

      // Here we want to make sure that if we apply Schur methods as
      // on p. 604, we can get the same answer as in the general delta.
      Matrix<double> U = submatrix( hessian, 0, 0,
                                    num_cam_entries, num_cam_entries);
      Matrix<double> W = submatrix( hessian, 0, num_cam_entries,
                                    num_cam_entries, num_pt_entries);
      Matrix<double> Vinv = submatrix(hessian, num_cam_entries, num_cam_entries,
                                      num_pt_entries, num_pt_entries);
      chol_inverse(Vinv);
      Matrix<double> Y = W * transpose(Vinv) * Vinv;
      Vector<double> e = subvector(delta, 0, num_cam_entries) -
        W * transpose(Vinv) * Vinv * subvector(delta, num_cam_entries, num_pt_entries);
      Matrix<double> S = U - Y * transpose(W);

      // set S
      this->set_S(S);
      solve(e, S); // using cholesky
      solve(delta, hessian);

      double nsq_x = 0;
      for ( unsigned j=0; j<this->m_model.num_cameras(); ++j )
        nsq_x += norm_2( this->m_model.A_parameters(j) );
      for ( unsigned i=0; i<this->m_model.num_points(); ++i )
        nsq_x += norm_2( this->m_model.B_parameters(i) );

      // --- EVALUATE POTENTIAL UPDATE STEP ---
      Vector<double> new_error(num_observations);                  // Error vector
      idx = 0;
      for (unsigned i = 0; i < this->m_control_net->size(); ++i) {
        for (unsigned m = 0; m < (*this->m_control_net)[i].size(); ++m) {
          int camera_idx = (*this->m_control_net)[i][m].image_id();

          Vector<double> cam_delta = subvector(delta, num_cam_params*camera_idx, num_cam_params);
          Vector<double> pt_delta = subvector(delta, num_cam_params*num_cameras + num_pt_params*i,
                                              num_pt_params);

          // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error;
          try {
            unweighted_error = (*this->m_control_net)[i][m].dominant() -
              this->m_model(i, camera_idx,
                            this->m_model.A_parameters(camera_idx)-cam_delta,
                            this->m_model.B_parameters(i)-pt_delta);
          } catch ( camera::PixelToRayErr &e ) {}
          double mag = norm_2(unweighted_error);
          double weight = sqrt(this->m_robust_cost_func(mag)) / mag;
          subvector(new_error,2*idx,2) = unweighted_error * weight;

          ++idx;
        }
      }

      // Add rows to J and error for a priori position/pose constraints...
      if (this->m_use_camera_constraint)
        for (unsigned j=0; j < num_cameras; ++j) {
          Vector<double> cam_delta = subvector(delta, num_cam_params*j, num_cam_params);
          subvector(new_error,
                    2*this->m_model.num_pixel_observations() + j*num_cam_params,
                    num_cam_params) = this->m_model.A_target(j)-
            (this->m_model.A_parameters(j)-cam_delta);
        }

      // ... and the position of the 3D points to J and error ...
      if (this->m_use_gcp_constraint) {
        idx = 0;
        for (unsigned i=0; i < this->m_model.num_points(); ++i) {
          if ((*this->m_control_net)[i].type() == ControlPoint::GroundControlPoint) {
            Vector<double> pt_delta = subvector(delta, num_cam_params*num_cameras + num_pt_params*i,
                                                num_pt_params);
            subvector(new_error,
                      2*this->m_model.num_pixel_observations() + num_cameras*num_cam_params +
                      idx*num_pt_params,
                      num_pt_params) = this->m_model.B_target(i)-
              (this->m_model.B_parameters(i)-pt_delta);
            ++idx;
          }
        }
      }

      //Fletcher modification
      double Splus = .5*transpose(new_error)*sigma*new_error; //new objective
      double SS = .5*transpose(error)*sigma*error;            //old objective
      double dS = .5 * transpose(delta) *(this->m_lambda * delta + epsilon);

      double R = (SS - Splus)/dS; // Compute ratio

      rel_tol = math::max(abs(epsilon) );
      abs_tol = Splus;

      if (R > 0) {

        for (unsigned j=0; j<this->m_model.num_cameras(); ++j)
          this->m_model.set_A_parameters(j, this->m_model.A_parameters(j) -
                                         subvector(delta, num_cam_params*j, num_cam_params));
        for (unsigned i=0; i<this->m_model.num_points(); ++i)
          this->m_model.set_B_parameters(i, this->m_model.B_parameters(i) -
                                         subvector(delta, num_cam_params*num_cameras + num_pt_params*i,
                                                   num_pt_params));

        if (this->m_control==0){
          double temp = 1 - pow((2*R - 1),3);
          if (temp < 1.0/3.0)
            temp = 1.0/3.0;

          this->m_lambda *= temp;
          this->m_nu = 2;
        } else if (this->m_control == 1)
          this->m_lambda /= 10;

        return SS-Splus;
      }

      // Didn't make progress ...
      if (this->m_control == 0){
        this->m_lambda *= this->m_nu;
        this->m_nu*=2;
      } else if (this->m_control == 1)
        this->m_lambda *= 10;

      return 0;
    }

  };

}}

#endif//__VW_BUNDLEADJUSTMENT_ADJUST_REF_H__
