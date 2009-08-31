// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

/// \file BundleAdjustmentRobustRef.h
/// 
/// Reference implementation of robust bundle adjustment. Very slow!

#ifndef __VW_CAMERA_BUNDLE_ADJUSTMENT_ROBUST_REF_H__
#define __VW_CAMERA_BUNDLE_ADJUSTMENT_ROBUST_REF_H__

#include <vw/Camera/BundleAdjustmentBase.h>

namespace vw {
namespace camera {
  
  template <class BundleAdjustModelT, class RobustCostT>
  class BundleAdjustmentRobustRef : public BundleAdjustmentBase<BundleAdjustModelT,RobustCostT> {

  public:

    BundleAdjustmentRobustRef( BundleAdjustModelT & model,
			 RobustCostT const& robust_cost_func,
			 bool use_camera_constraint=true,
			 bool use_gcp_constraint=true ) :
    BundleAdjustmentBase<BundleAdjustModelT,RobustCostT>( model, robust_cost_func,
							  use_camera_constraint,
							  use_gcp_constraint ) {}
    
   
    // UPDATE IMPLEMENTATION
    //---------------------------------------------------------------
    // This is a simple, non-sparse, unoptimized implementation of LM
    // bundle adjustment.  It is primarily used for validation and
    // debugging.
    //
    // Each entry in the outer vector corresponds to a distinct 3D
    // point.  The inner vector contains a list of image IDs and
    // pixel coordinates where that point was imaged.
    virtual double update(double &abs_tol, double &rel_tol) { 
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
      unsigned num_pixel_observations = 2*this->m_model.num_pixel_observations();
      unsigned num_pixels = this->m_model.num_pixel_observations();

      unsigned num_observations = 2*this->m_model.num_pixel_observations();
      if (this->m_use_camera_constraint)
	num_observations += num_cameras*num_cam_params;
      if (this->m_use_gcp_constraint)
	num_observations += num_ground_control_points*num_pt_params;

      // The core LM matrices and vectors
      Matrix<double> J_robust(num_observations, num_model_parameters); // Jacobian Matrix for Robust algorithm

      Vector<double> epsilon_robust(num_observations);            // Modified error vector for Robust algorithm 
      
      double t_df = 4;                                            // Degrees of freedom for data (can be modified later)
      double t_dim = 2;                                           // dimension of pixels


      Matrix<double> sigma(num_observations, num_observations);   // Sigma (uncertainty) matrix

      double robust_objective = 0.0;                               // robust objective



      // --- SETUP STEP ----
      // Add rows to J and epsilon for the imaged pixel observations
     

      int idx = 0;
      for (unsigned i = 0; i < this->m_control_net->size(); ++i) {       // Iterate over control points
        for (unsigned m = 0; m < (*(this->m_control_net))[i].size(); ++m) {  // Iterate over control measures
          int camera_idx = (*(this->m_control_net))[i][m].image_id();

          Matrix<double> J_a = this->m_model.A_jacobian(i,camera_idx, 
                                                  this->m_model.A_parameters(camera_idx), 
                                                 this-> m_model.B_parameters(i));
          Matrix<double> J_b = this->m_model.B_jacobian(i,camera_idx, 
                                                  this->m_model.A_parameters(camera_idx), 
                                                  this->m_model.B_parameters(i));

	  // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error = (*(this->m_control_net))[i][m].dominant() - this->m_model(i, camera_idx, 
                                                                              this->m_model.A_parameters(camera_idx), 
                                                                              this->m_model.B_parameters(i));
	  // Fill in the entries of the sigma matrix with the uncertainty of the observations.
          Matrix2x2 inverse_cov;
          Vector2 pixel_sigma = (*(this->m_control_net))[i][m].sigma();
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));
          submatrix(sigma, 2*idx, 2*idx, 2, 2) = inverse_cov;

	  // Populate the S_weights, mu_weights vectors
	  double S_weight = transpose(unweighted_error) * inverse_cov * unweighted_error;
	  double mu_weight = (t_df + t_dim)/(t_df + S_weight);

	  robust_objective += 0.5*(t_df + t_dim)*log(1 + S_weight/t_df);

	  
	  // Populate the robust epsilon vector
	  subvector(epsilon_robust,2*idx,2) = unweighted_error * sqrt(mu_weight);
	  
	  // Populate the robust Jacobian Matrix 
          submatrix(J_robust, 2*idx, num_cam_params*camera_idx, 2, num_cam_params) = sqrt(mu_weight)*J_a;
          submatrix(J_robust, 2*idx, num_cam_params*num_cameras + i*num_pt_params, 2, num_pt_params) = sqrt(mu_weight)*J_b;


          ++idx;
        }
      }
      
      // Add rows to J and epsilon for a priori camera parameters...
      if (this->m_use_camera_constraint)
	for (unsigned j=0; j < num_cameras; ++j) {
	  Matrix<double> id(num_cam_params, num_cam_params);
	  id.set_identity();
	
	  Vector<double> unweighted_error = this->m_model.A_initial(j)-this->m_model.A_parameters(j);
	 
	  double S_weight = transpose(unweighted_error) * (this->m_model.A_inverse_covariance(j)) * unweighted_error;
	  double mu_weight = (t_df + t_dim)/(t_df + S_weight);

	  robust_objective += 0.5*(t_df + t_dim)*log(1 + S_weight/t_df);

	  // Here the J_robust is modified exactly as J was
	  submatrix(J_robust,
		    2*this->m_model.num_pixel_observations() + j*num_cam_params,
		    j*num_cam_params,
		    num_cam_params,
		    num_cam_params) = id*sqrt(mu_weight);


	  // Here epsilon_robust is modified exactly as epsilon was
	  subvector(epsilon_robust,
		    2*this->m_model.num_pixel_observations() + j*num_cam_params,
		    num_cam_params) = unweighted_error*sqrt(mu_weight);

	  submatrix(sigma,
		    2*this->m_model.num_pixel_observations() + j*num_cam_params,
		    2*this->m_model.num_pixel_observations() + j*num_cam_params,
		    num_cam_params, num_cam_params) = this->m_model.A_inverse_covariance(j);

	  // add the initial constraints on the regular quadratic scale
	}

      // ... and the position of the 3D points to J and epsilon ...
      if (this->m_use_gcp_constraint) {
	idx = 0;
	for (unsigned i=0; i < this->m_model.num_points(); ++i) {
	  if ((*(this->m_control_net))[i].type() == ControlPoint::GroundControlPoint) {
	    Matrix<double> id(num_pt_params,num_pt_params);
	    id.set_identity();


	    Vector<double> unweighted_error = this->m_model.B_initial(i)-this->m_model.B_parameters(i);
	    // Here the J_robust is modified exactly as J was

	    double S_weight = transpose(unweighted_error) * this->m_model.B_inverse_covariance(i) * unweighted_error;
	    double mu_weight = (t_df + t_dim)/(t_df + S_weight);

	    robust_objective += 0.5*(t_df + t_dim)*log(1 + S_weight/t_df);
	    

	    submatrix(J_robust, 2*this->m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
		      num_cameras*num_cam_params + idx*num_pt_params,
		      num_pt_params,
		      num_pt_params) = id*sqrt(mu_weight);
	    

	    
	    
	    // Here epsilon_robust is modified exactly as epsilon was
	    subvector(epsilon_robust,
		      2*this->m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
		      num_pt_params) = unweighted_error*sqrt(mu_weight);

	    submatrix(sigma,
		      2*this->m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
		      2*this->m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
		      num_pt_params, num_pt_params) = this->m_model.B_inverse_covariance(i);
	    
	    // Add 3D points on the regular scale
	  

	    ++idx;
	  }
	}
      }


    
      // --- UPDATE STEP ----


// Build up the right side of the normal equation for robust algorithm

      Matrix<double> JTS = transpose(J_robust) * sigma;

      Vector<double> del_J_robust = -1.0 * (JTS * epsilon_robust);

    

      // also build up left side for robust
      Matrix<double> hessian_robust = JTS * J_robust;

    
      // initialize m_lambda on first iteration, ignore if user has
      // changed it.
      double max = 0.0;
      if (this->m_iterations == 1 && this->m_lambda == 1e-3){
	for (unsigned i = 0; i < hessian_robust.rows(); ++i){
	  if (fabs(hessian_robust(i,i)) > max)
	    max = fabs(hessian_robust(i,i));
	}
	this->m_lambda = 1e-10 * max;
      }

      for ( unsigned i=0; i < hessian_robust.rows(); ++i )
        hessian_robust(i,i) +=  this->m_lambda;

      //Cholesky decomposition. Returns Cholesky matrix in lower left hand corner.
      Vector<double> delta_robust = del_J_robust;


      // Here we want to make sure that if we apply Schur methods as on p. 604, we can get the same answer as in the general delta.  
      unsigned num_cam_entries = num_cam_params * num_cameras;
      unsigned num_pt_entries = num_pt_params * num_points;

      Matrix<double> U_robust = submatrix(hessian_robust, 0, 0, num_cam_entries, num_cam_entries);
      Matrix<double> W_robust = submatrix(hessian_robust, 0, num_cam_entries,  num_cam_entries, num_pt_entries);
      Matrix<double> Vinv_robust = submatrix(hessian_robust, num_cam_entries, num_cam_entries, num_pt_entries, num_pt_entries); 
      chol_inverse(Vinv_robust);

      Matrix<double> Y_robust = W_robust * transpose(Vinv_robust) * Vinv_robust;

      Vector<double> e_robust = subvector(delta_robust, 0, num_cam_entries) - 
	W_robust * transpose(Vinv_robust) * Vinv_robust * subvector(delta_robust, num_cam_entries, num_pt_entries);

      Matrix<double> S_robust = U_robust - Y_robust * transpose(W_robust);

    
      solve(e_robust, S_robust); // using cholesky

      solve(delta_robust, hessian_robust);

      
      // Solve for update

      // --- EVALUATE UPDATE STEP ---
      Vector<double> new_epsilon_robust(num_observations);           // Robust Error vector
      
      double new_robust_objective = 0.0;

      idx = 0;
      for (unsigned i = 0; i < this->m_control_net->size(); ++i) {         // Iterate over control points
        for (unsigned m = 0; m < (*(this->m_control_net))[i].size(); ++m) {  // Iterate over control measures
          int camera_idx = (*(this->m_control_net))[i][m].image_id();

      	  Vector<double> cam_delta_robust = subvector(delta_robust, num_cam_params*camera_idx, num_cam_params);

      	  Vector<double> pt_delta_robust = subvector(delta_robust, num_cam_params*num_cameras + num_pt_params*i, num_pt_params);


          // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error = (*(this->m_control_net))[i][m].dominant() - this->m_model(i, camera_idx,
                                                                              this->m_model.A_parameters(camera_idx)-cam_delta_robust, 
                                                                             this-> m_model.B_parameters(i)-pt_delta_robust);
	  Matrix2x2 inverse_cov = submatrix(sigma, 2*idx, 2*idx, 2, 2);

	  // Populate the S_weights, mu_weights vectors
	  double S_weight = transpose(unweighted_error) * inverse_cov * unweighted_error;
	  double mu_weight = (t_df + t_dim)/(t_df + S_weight);


	  // Populate the robust epsilon vector
	  subvector(new_epsilon_robust,2*idx,2) = unweighted_error * sqrt(mu_weight);
	  
	  new_robust_objective += 0.5*(t_df + t_dim)*log(1 + S_weight/t_df);
        
          ++idx;
        }
      }

      // Add rows to J and epsilon for a priori position/pose constraints...
      if (this->m_use_camera_constraint)
	for (unsigned j=0; j < num_cameras; ++j) {
	  Vector<double> cam_delta_robust = subvector(delta_robust, num_cam_params*j, num_cam_params);

	  Vector<double> unweighted_error = this->m_model.A_initial(j)- this->m_model.A_parameters(j);
	  
	  double S_weight = transpose(unweighted_error) * this->m_model.A_inverse_covariance(j) * unweighted_error;
	  double mu_weight = (t_df + t_dim)/(t_df + S_weight);

	  subvector(new_epsilon_robust,
		    2*this->m_model.num_pixel_observations() + j*num_cam_params,
		    num_cam_params) = unweighted_error*sqrt(mu_weight);
	  
	  new_robust_objective += 0.5*(t_df + t_dim)*log(1 + S_weight/t_df);
	  
	}

      // ... and the position of the 3D points to J and epsilon ...
      if (this->m_use_gcp_constraint) {
	idx = 0;
	for (unsigned i=0; i < this->m_model.num_points(); ++i) {
	  if ((*(this->m_control_net))[i].type() == ControlPoint::GroundControlPoint) {
	    Vector<double> pt_delta_robust = subvector(delta_robust, num_cam_params*num_cameras + num_pt_params*i, num_pt_params);

	    Vector<double> unweighted_error = this->m_model.B_initial(i)-this->m_model.B_parameters(i);

	    double S_weight = transpose(unweighted_error)*this->m_model.B_inverse_covariance(i)*unweighted_error;
	    double mu_weight = (t_df + t_dim)/(t_df + S_weight);
	  
	    subvector(new_epsilon_robust,
		      2*this->m_model.num_pixel_observations() + num_cameras*num_cam_params + idx*num_pt_params,
		      num_pt_params) = unweighted_error*sqrt(mu_weight);

	    new_robust_objective += 0.5*(t_df + t_dim)*log(1 + S_weight/t_df);
	  
	    ++idx;
	  }
	}
      }



      //Fletcher modification for robust case
      double dS_robust = .5 * transpose(delta_robust)*(this->m_lambda*delta_robust + del_J_robust);

      double R_robust = (robust_objective - new_robust_objective)/dS_robust;

      unsigned ret = 0;

      std::cout << " Old robust objective: " << robust_objective << "\n";
      std::cout << " New robust objective: " << new_robust_objective << "\n";
      std::cout << " Lambda " << this->m_lambda << "\n";
         
      if (R_robust > 0){
	ret = 1;	
	for (unsigned j=0; j<this->m_model.num_cameras(); ++j) 
          this->m_model.set_A_parameters(j, this->m_model.A_parameters(j) - subvector(delta_robust, num_cam_params*j, num_cam_params));
        for (unsigned i=0; i<this->m_model.num_points(); ++i)
          this->m_model.set_B_parameters(i, this->m_model.B_parameters(i) - subvector(delta_robust, num_cam_params*num_cameras + num_pt_params*i, num_pt_params));
	
	// Summarize the stats from this step in the iteration
        double overall_norm = sqrt(.5 * transpose(new_epsilon_robust) * sigma * new_epsilon_robust);
	
        //double overall_delta = sqrt(.5 * transpose(epsilon) * sigma * epsilon) - sqrt(.5 * transpose(new_epsilon) * sigma * new_epsilon) ; 
	double overall_delta = robust_objective - new_robust_objective; 


	abs_tol = overall_norm;
        rel_tol = overall_delta;	

	if (this->m_control==0){
	  double temp = 1 - pow((2*R_robust - 1),3);		
	  if (temp < 1.0/3.0)
	    temp = 1.0/3.0;
	  
	  this->m_lambda *= temp;
	  this->m_nu = 2;
	} else if (this->m_control == 1)
	  this->m_lambda /= 10;
	
	return overall_delta;
      
      } else { // R_robust <= 0
	
	if (this->m_control == 0){
	  this->m_lambda *= this->m_nu; 
	  this->m_nu*=2; 
	} else if (this->m_control == 1)
	  this->m_lambda *= 10;

	// double overall_delta = sqrt(.5 * transpose(epsilon) * sigma * epsilon) - sqrt(.5 * transpose(new_epsilon) * sigma * new_epsilon);
 
	return ScalarTypeLimits<double>::highest();
      }
    }





  };

}}

#endif//__VW_CAMERA_BUNDLE_ADJUSTMENT_ROBUST_REF_H__
