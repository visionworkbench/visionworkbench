// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file AdjustRobustRef.h
///
/// Reference implementation of robust bundle adjustment. Very slow!

#ifndef __VW_BUNDLEADJUSTMENT_ADJUST_ROBUST_REF_H__
#define __VW_BUNDLEADJUSTMENT_ADJUST_ROBUST_REF_H__

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
  class AdjustRobustRef : public AdjustBase<BundleAdjustModelT,RobustCostT>, private boost::noncopyable {

    // Need to save S for covariance calculations
    math::Matrix<double> m_S;

  public:

    AdjustRobustRef( BundleAdjustModelT & model,
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
      for(unsigned i = 0; i < num_cameras; i++){
         sparse_cov(i) = submatrix(Cov, i*num_cam_params, i*num_cam_params, num_cam_params, num_cam_params);
      }

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
      unsigned num_model_parameters = this->m_model.num_cameras()*num_cam_params + this->m_model.num_points()*num_pt_params;

      unsigned num_cameras = this->m_model.num_cameras();
      unsigned num_ground_control_points = this->m_control_net->num_ground_control_points();

      // Need to keep pixel obs separate from 'initials' for robust code
      // unsigned num_pixel_observations = 2*this->m_model.num_pixel_observations();
      // unsigned num_pixels = this->m_model.num_pixel_observations(); // delete ?

      unsigned num_observations = 2*this->m_model.num_pixel_observations();
      if (this->m_use_camera_constraint)
        num_observations += num_cameras*num_cam_params;
      if (this->m_use_gcp_constraint)
        num_observations += num_ground_control_points*num_pt_params;

      // Jacobian Matrix for Robust algorithm
      Matrix<double> J(num_observations, num_model_parameters);
      // Modified error vector for Robust algorithm
      Vector<double> error(num_observations);
      // Degrees of freedom for data (can be modified later)
      double t_df = 4;
      // dimension of pixels
      double t_dim_pixel = 2;
      // dimension of camera params
      double t_dim_cam = 6;
      // dimension of world point
      double t_dim_pt = 3;
      // Sigma (uncertainty) matrix
      Matrix<double> sigma(num_observations, num_observations);
      // robust objective
      double robust_objective = 0.0;

      // --- SETUP STEP ----
      // Add rows to J and error for the imaged pixel observations

      int idx = 0;
      // Iterate over control points
      for (unsigned i = 0; i < this->m_control_net->size(); ++i) {
        // Iterate over control measures
        for (unsigned m = 0; m < (*(this->m_control_net))[i].size(); ++m) {
          int camera_idx = (*(this->m_control_net))[i][m].image_id();

          Matrix<double> J_a = this->m_model.A_jacobian(i,camera_idx,
                                                        this->m_model.A_parameters(camera_idx),
                                                        this-> m_model.B_parameters(i));
          Matrix<double> J_b = this->m_model.B_jacobian(i,camera_idx,
                                                        this->m_model.A_parameters(camera_idx),
                                                        this->m_model.B_parameters(i));

          // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error;
          try {
            unweighted_error = (*(this->m_control_net))[i][m].dominant() -
              this->m_model(i, camera_idx,
                            this->m_model.A_parameters(camera_idx),
                            this->m_model.B_parameters(i));
          } catch ( camera::PixelToRayErr &e ) {}
          // Fill in the entries of the sigma matrix with the uncertainty of the observations.
          Matrix2x2 inverse_cov;
          Vector2 pixel_sigma = (*(this->m_control_net))[i][m].sigma();
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));
          submatrix(sigma, 2*idx, 2*idx, 2, 2) = inverse_cov;

          // Populate the S_weights, mu_weights vectors
          double S_weight = transpose(unweighted_error) * inverse_cov * unweighted_error;
          double mu_weight = (t_df + t_dim_pixel)/(t_df + S_weight);

          robust_objective += 0.5*(t_df + t_dim_pixel)*log(1 + S_weight/t_df);


          // Populate the robust error vector
          subvector(error,2*idx,2) = unweighted_error * sqrt(mu_weight);

          // Populate the robust Jacobian Matrix
          submatrix(J, 2*idx, num_cam_params*camera_idx, 2, num_cam_params) = sqrt(mu_weight)*J_a;
          submatrix(J, 2*idx, num_cam_params*num_cameras + i*num_pt_params, 2, num_pt_params) = sqrt(mu_weight)*J_b;


          ++idx;
        }
      }

      // initialize m_lambda on first iteration, ignore if user has
      // changed it.
      double max = 0.0;
      if ( this->m_iterations == 1 && this->m_lambda == 1e-3 ) {
        Matrix<double> hessian = transpose(J)*sigma*J;
        for ( unsigned i = 0; i < hessian.rows(); ++i )
          if ( fabs(hessian(i,i)) > max )
            max = fabs(hessian(i,i));
        this->m_lambda = 1e-10 * max;
      }

      // Add rows to J and error for a priori camera parameters...
      if (this->m_use_camera_constraint)
        for ( unsigned j=0; j < num_cameras; ++j ) {
          Matrix<double> id( num_cam_params, num_cam_params );
          id.set_identity();

          Vector<double> unweighted_error = this->m_model.A_target(j)-this->m_model.A_parameters(j);

          double S_weight = transpose(unweighted_error) * (this->m_model.A_inverse_covariance(j)) * unweighted_error;
          double mu_weight = (t_df + t_dim_cam)/(t_df + S_weight);

          robust_objective += 0.5*(t_df + t_dim_cam)*log(1 + S_weight/t_df);

          // Here the J is modified exactly as J was
          submatrix(J,
                    2*this->m_model.num_pixel_observations() + j*num_cam_params,
                    j*num_cam_params,
                    num_cam_params,
                    num_cam_params) = id*sqrt(mu_weight);

          // Here error is modified exactly as error was
          subvector(error,
                    2*this->m_model.num_pixel_observations() + j*num_cam_params,
                    num_cam_params) = unweighted_error*sqrt(mu_weight);

          submatrix(sigma,
                    2*this->m_model.num_pixel_observations() + j*num_cam_params,
                    2*this->m_model.num_pixel_observations() + j*num_cam_params,
                    num_cam_params, num_cam_params) = this->m_model.A_inverse_covariance(j);

          // add the initial constraints on the regular quadratic scale
        }

      // ... and the position of the 3D points to J and error ...
      if (this->m_use_gcp_constraint) {
        idx = 0;
        for (unsigned i=0; i < this->m_model.num_points(); ++i) {
          if ((*(this->m_control_net))[i].type() == ControlPoint::GroundControlPoint) {
            Matrix<double> id(num_pt_params,num_pt_params);
            id.set_identity();

            Vector<double> unweighted_error = this->m_model.B_target(i)-this->m_model.B_parameters(i);
            // Here the J is modified exactly as J was
            double S_weight = transpose(unweighted_error) * this->m_model.B_inverse_covariance(i) * unweighted_error;
            double mu_weight = (t_df + t_dim_pt)/(t_df + S_weight);

            robust_objective += 0.5*(t_df + t_dim_pt)*log(1 + S_weight/t_df);

            submatrix(J, 2*this->m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                      num_cameras*num_cam_params + idx*num_pt_params,
                      num_pt_params,
                      num_pt_params) = id*sqrt(mu_weight);

            // Here error is modified exactly as error was
            subvector(error,
                      2*this->m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                      num_pt_params) = unweighted_error*sqrt(mu_weight);

            submatrix(sigma,
                      2*this->m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                      2*this->m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
                      num_pt_params, num_pt_params) = this->m_model.B_inverse_covariance(i);

             ++idx;
          }
        }
      }

      // --- SOLVE UPDATE STEP ----------------------------------------

      // Build up the right side of the normal equation for robust algorithm
      Vector<double> epsilon = -1.0 * transpose(J) * sigma * error;

      // also build up left side for robust
      Matrix<double> hessian = transpose(J) * sigma * J;

      for ( unsigned i=0; i < hessian.rows(); ++i )
        hessian(i,i) +=  this->m_lambda;

      // Cholesky decomposition. Returns Cholesky matrix in lower left
      // hand corner.
      Vector<double> delta = epsilon;

      // Here we want to make sure that if we apply Schur methods as on p. 604, we can get the same answer as in the general delta.
      unsigned num_cam_entries = num_cam_params * num_cameras;
      unsigned num_pt_entries = num_pt_params * num_points;

      Matrix<double> U = submatrix( hessian, 0, 0,
                                    num_cam_entries, num_cam_entries);
      Matrix<double> W = submatrix( hessian, 0, num_cam_entries,
                                    num_cam_entries, num_pt_entries);
      Matrix<double> Vinv = submatrix( hessian, num_cam_entries, num_cam_entries,
                                       num_pt_entries, num_pt_entries);
      chol_inverse(Vinv);
      Matrix<double> Y = W * transpose(Vinv) * Vinv;
      Vector<double> e = subvector(delta, 0, num_cam_entries) -
        W * transpose(Vinv) * Vinv * subvector(delta, num_cam_entries, num_pt_entries);

      Matrix<double> S = U - Y * transpose(W);

      // Set S
      this->set_S(S);
      solve(e, S); // using cholesky
      solve(delta, hessian);

      // Solve for update
      double new_objective = 0.0;

      idx = 0;
      // Iterate over control points
      for (unsigned i = 0; i < this->m_control_net->size(); ++i) {
        // Iterate over control measures
        for (unsigned m = 0; m < (*(this->m_control_net))[i].size(); ++m) {
          int camera_idx = (*(this->m_control_net))[i][m].image_id();

          Vector<double> cam_delta = subvector(delta, num_cam_params*camera_idx, num_cam_params);

          Vector<double> pt_delta = subvector(delta, num_cam_params*num_cameras + num_pt_params*i, num_pt_params);

          // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error;
          try {
            unweighted_error = (*(this->m_control_net))[i][m].dominant() -
              this->m_model(i, camera_idx,
                            this->m_model.A_parameters(camera_idx)-cam_delta,
                            this-> m_model.B_parameters(i)-pt_delta);
          } catch ( camera::PixelToRayErr &e ) {}
          Matrix2x2 inverse_cov = submatrix(sigma, 2*idx, 2*idx, 2, 2);

          // Populate the S_weights, mu_weights vectors
          double S_weight = transpose(unweighted_error) * inverse_cov * unweighted_error;

          new_objective += 0.5*(t_df + t_dim_pixel)*log(1 + S_weight/t_df);

          ++idx;
        }
      }

      // Add rows to J and error for a priori position/pose constraints...
      if (this->m_use_camera_constraint)
        for (unsigned j=0; j < num_cameras; ++j) {
          Vector<double> cam_delta = subvector(delta, num_cam_params*j, num_cam_params);

          Vector<double> unweighted_error = this->m_model.A_target(j)- (this->m_model.A_parameters(j) - cam_delta);

          double S_weight = transpose(unweighted_error) * this->m_model.A_inverse_covariance(j) * unweighted_error;

           new_objective += 0.5*(t_df + t_dim_cam)*log(1 + S_weight/t_df);
        }

      // ... and the position of the 3D points to J and error ...
      if (this->m_use_gcp_constraint) {
        idx = 0;
        for (unsigned i=0; i < this->m_model.num_points(); ++i) {
          if ((*(this->m_control_net))[i].type() == ControlPoint::GroundControlPoint) {
            Vector<double> pt_delta = subvector(delta, num_cam_params*num_cameras + num_pt_params*i, num_pt_params);

            Vector<double> unweighted_error = this->m_model.B_target(i)-(this->m_model.B_parameters(i) - pt_delta);

            double S_weight = transpose(unweighted_error)*this->m_model.B_inverse_covariance(i)*unweighted_error;

            new_objective += 0.5*(t_df + t_dim_pt)*log(1 + S_weight/t_df);

            ++idx;
          }
        }
      }

      //Fletcher modification for robust case
      double dS = .5 * transpose(delta)*(this->m_lambda*delta + epsilon);
      double R = (robust_objective - new_objective)/dS;

      rel_tol = math::max(abs(epsilon));
      abs_tol = new_objective;

      if ( R > 0 ) {

        for (unsigned j=0; j<this->m_model.num_cameras(); ++j)
          this->m_model.set_A_parameters(j, this->m_model.A_parameters(j) - subvector(delta, num_cam_params*j, num_cam_params));
        for (unsigned i=0; i<this->m_model.num_points(); ++i)
          this->m_model.set_B_parameters(i, this->m_model.B_parameters(i) - subvector(delta, num_cam_params*num_cameras + num_pt_params*i, num_pt_params));


        if (this->m_control==0){
          double temp = 1 - pow((2*R - 1),3);
          if (temp < 1.0/3.0)
            temp = 1.0/3.0;

          this->m_lambda *= temp;
          this->m_nu = 2;
        } else if (this->m_control == 1)
          this->m_lambda /= 10;

        return robust_objective-new_objective;
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

#endif//__VW_BUNDLEADJUSTMENT_ADJUST_ROBUST_REF_H__
