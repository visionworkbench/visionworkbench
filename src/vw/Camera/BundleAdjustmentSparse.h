// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

/// \file BundleAdjustmentSparse.h
/// 
/// Sparse implementation of bundle adjustment. Fast yo!

#ifndef __VW_CAMERA_BUNDLE_ADJUSTMENT_SPARSE_H__
#define __VW_CAMERA_BUNDLE_ADJUSTMENT_SPARSE_H__

// Vision Workbench
#include <vw/Camera/BundleAdjustmentBase.h>
#include <vw/Math/SparseSkylineMatrix.h>

// Boost 
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/make_shared.hpp>
#include <boost/version.hpp>
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
namespace camera {

  template <class BundleAdjustModelT, class RobustCostT>
  class BundleAdjustmentSparse : public BundleAdjustmentBase<BundleAdjustModelT, RobustCostT> {

    // Need to save S for covariance calculations
    boost::shared_ptr<math::SparseSkylineMatrix<double> > m_S; 

  public:
    math::SparseSkylineMatrix<double> S() { return *m_S; }
    void set_S(math::SparseSkylineMatrix<double> S) { 
      m_S = boost::make_shared<math::SparseSkylineMatrix<double> >(S); 
    }

    BundleAdjustmentSparse( BundleAdjustModelT & model,
                            RobustCostT const& robust_cost_func,
                            bool use_camera_constraint=true,
                            bool use_gcp_constraint=true) :
    BundleAdjustmentBase<BundleAdjustModelT,RobustCostT>( model, robust_cost_func,
                                                          use_camera_constraint,
                                                          use_gcp_constraint ) {}

   




     // Covariance Calculator
    // ___________________________________________________________
    // This routine inverts a sparse matrix S, and prints the individual 
    // covariance matrices for each camera

      virtual void covCalc(){

      // camera params
      unsigned num_cam_params = BundleAdjustModelT::camera_params_n;
      unsigned num_cameras = this->m_model.num_cameras();

      unsigned inverse_size = num_cam_params * num_cameras;

      typedef Matrix<double, BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n> matrix_camera_camera;
    
      // final vector of camera covariance matrices
      boost_sparse_vector< matrix_camera_camera > sparse_cov(num_cameras);
      
      
      
      // Get the S matrix from the model
      math::SparseSkylineMatrix<double> S = this->S();
      
      
       Matrix<double> Id(inverse_size, inverse_size);
       Id.set_identity();
       
       Matrix<double> Cov = multi_sparse_solve(S, Id);

   
       //pick out covariances of individual cameras
       for(int i = 0; i < num_cameras; i++){
	 sparse_cov(i) = submatrix(Cov, i*num_cam_params, i*num_cam_params, num_cam_params, num_cam_params);
       }

       std::cout << "Covariance matrices for cameras are:" << sparse_cov << "\n\n";
      
       return;
       }





    // UPDATE IMPLEMENTATION
    //-------------------------------------------------------------
    // This is the sparse levenberg marquardt update step.  Returns
    // the average improvement in the cost function.
    virtual double update(double &abs_tol, double &rel_tol) {
      ++this->m_iterations;
      
      VW_DEBUG_ASSERT(this->m_control_net->size() == this->m_model.num_points(), LogicErr() << "BundleAdjustment::update() : Number of bundles does not match the number of points in the bundle adjustment model.");

      // Jacobian Matrices and error values
      typedef Matrix<double, 2, BundleAdjustModelT::camera_params_n> matrix_2_camera;
      typedef Matrix<double, 2, BundleAdjustModelT::point_params_n> matrix_2_point;
      boost_sparse_matrix< matrix_2_camera  > A(this->m_model.num_points(), this->m_model.num_cameras());
      boost_sparse_matrix< matrix_2_point > B(this->m_model.num_points(), this->m_model.num_cameras());
      boost_sparse_matrix<Vector2> epsilon(this->m_model.num_points(), this->m_model.num_cameras());
      boost_sparse_matrix<Vector2> new_epsilon(this->m_model.num_points(), this->m_model.num_cameras());
      
      // Data structures necessary for Fletcher modification
      // boost::numeric::ublas::mapped_matrix<Vector2> Jp(m_model.num_points(), m_model.num_cameras());
           
      // Intermediate Matrices and vectors
      typedef Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::camera_params_n> matrix_camera_camera;
      typedef Matrix<double,BundleAdjustModelT::point_params_n,BundleAdjustModelT::point_params_n> matrix_point_point;
      typedef Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::point_params_n> matrix_camera_point;
      boost_sparse_vector< matrix_camera_camera > U(this->m_model.num_cameras());
      boost_sparse_vector< matrix_point_point > V(this->m_model.num_points());  
      boost_sparse_matrix< matrix_camera_point > W(this->m_model.num_cameras(), this->m_model.num_points());     
     
      // Copies of Intermediate Marices
      typedef Vector<double,BundleAdjustModelT::camera_params_n> vector_camera;
      typedef Vector<double,BundleAdjustModelT::point_params_n> vector_point;
      boost_sparse_vector< vector_camera > epsilon_a(this->m_model.num_cameras());
      boost_sparse_vector< vector_point > epsilon_b(this->m_model.num_points());
      boost_sparse_matrix< matrix_camera_point > Y(this->m_model.num_cameras(), this->m_model.num_points());

      unsigned num_cam_params = BundleAdjustModelT::camera_params_n;
      unsigned num_pt_params = BundleAdjustModelT::point_params_n;

      unsigned delta_length = U.size() * num_cam_params + V.size() * num_pt_params;
      
      Vector<double> g(delta_length);  
      unsigned current_g_length = 0;
      
      Vector<double> delta(delta_length);
      unsigned current_delta_length = 0;

      // Fletcher LM parameteres
      double dS = 0; //Predicted improvement for Fletcher modification
      
      // Populate the Jacobian, which is broken into two sparse
      // matrices A & B, as well as the error matrix and the W 
      // matrix.
      vw_out(DebugMessage, "bundle_adjustment") << "Image Error: " << std::endl;
      unsigned i = 0;
      double error_total = 0; // assume this is r^T\Sigma^{-1}r


      for (typename ControlNetwork::const_iterator iter = this->m_control_net->begin(); 
	   iter != this->m_control_net->end(); ++iter) {
	for (typename ControlPoint::const_iterator measure_iter = (*iter).begin(); 
	     measure_iter != (*iter).end(); ++measure_iter) {


     
      //for (unsigned i = 0; i < this->m_control_net->size(); ++i) {       // Iterate over control points
      //for (unsigned m = 0; m < (*(this->m_control_net))[i].size(); 
      //     ++m) {  // Iterate over control measures




          unsigned j = (*measure_iter).image_id();
	  //unsigned j = (*(this->m_control_net))[i][m].image_id();

	   VW_DEBUG_ASSERT(j >=0 && j < this->m_model.num_cameras(), 
	                ArgumentErr() << "BundleAdjustment::update() : image index out of bounds.");
          
          // Store jacobian values
	   if ( i == 106144 ) std::cout << "A Jacobian!\n";
          A(i,j) = this->m_model.A_jacobian(i,j,this->m_model.A_parameters(j),
                                            this->m_model.B_parameters(i));
	   if ( i == 106144 ) std::cout << "B Jacobian!\n";
          B(i,j) = this->m_model.B_jacobian(i,j,this->m_model.A_parameters(j),
                                            this->m_model.B_parameters(i));

          // Apply robust cost function weighting
	  Vector2 unweighted_error = measure_iter->dominant() - 
	    this->m_model(i,j,this->m_model.A_parameters(j),this->m_model.B_parameters(i));

	  //Vector2 unweighted_error = (*(this->m_control_net))[i][m].dominant() - 
	  //this->m_model(i, j, 
	  //              this->m_model.A_parameters(j), 
	  //              this->m_model.B_parameters(i));

          double mag = norm_2(unweighted_error);
          double weight = sqrt(this->m_robust_cost_func(mag)) / mag;

          epsilon(i,j) = unweighted_error * weight;
	  
	  //  std::cout << "i: " << i << "j: " << j << "epsilon(i,j): " << unweighted_error * weight << "\n";
                            
          Matrix2x2 inverse_cov;

	  Vector2 pixel_sigma = measure_iter->sigma();
	  //Vector2 pixel_sigma = (*(this->m_control_net))[i][m].sigma();
        

          Vector2 epsilon_inst = epsilon(i,j);
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));
          error_total += .5 * transpose(epsilon_inst) * 
            inverse_cov * epsilon_inst;
          
          // Store intermediate values
          U(j) += transpose(static_cast< matrix_2_camera >(A(i,j))) * 
            inverse_cov * static_cast< matrix_2_camera >(A(i,j));
          V(i) += transpose(static_cast< matrix_2_point >(B(i,j))) * 
            inverse_cov * static_cast< matrix_2_point >(B(i,j));
          W(j,i) = transpose(static_cast< matrix_2_camera >(A(i,j))) * 
            inverse_cov * static_cast< matrix_2_point >(B(i,j));

          epsilon_a(j) += transpose(static_cast< matrix_2_camera >(A(i,j))) * 
            inverse_cov * epsilon_inst;
          epsilon_b(i) += transpose(static_cast< matrix_2_point >(B(i,j))) * 
            inverse_cov * epsilon_inst;

        }
	 ++i;
      }

 // set initial lambda, and ignore if the user has touched it
      if (this->m_iterations == 1 && this->m_lambda == 1e-3){
        double max = 0.0;
        for (unsigned i = 0; i < U.size(); ++i)  
          for (unsigned j = 0; j < BundleAdjustModelT::camera_params_n; ++j){
            if (fabs(static_cast<matrix_camera_camera>(U(i))(j,j)) > max)
              max = fabs(static_cast<matrix_camera_camera>(U(i))(j,j));
          }
        for (unsigned i = 0; i < V.size(); ++i) 
          for (unsigned j = 0; j < BundleAdjustModelT::point_params_n; ++j) {
            if ( fabs(static_cast<matrix_point_point>(V(i))(j,j)) > max)
              max = fabs(static_cast<matrix_point_point>(V(i))(j,j));
          }
        this->m_lambda = max * 1e-10;

	//	std::cout << "Computed lambda is: " << max*1e-10 << "\n";
	//std::cout << "other order: " << 1e-10*max << "\n";
      }
      // std::cout << "\n (new) U Matrix array prior to camera positions is: " << U << "\n";
    

      // Add in the camera position and pose constraint terms and covariances.
      if (this->m_use_camera_constraint) {
        for (unsigned j = 0; j < U.size(); ++j) {
          
          matrix_camera_camera inverse_cov;
          inverse_cov = this->m_model.A_inverse_covariance(j);
          
          matrix_camera_camera C;
	  C.set_identity();

          U(j) += transpose(C) * inverse_cov * C;
          
          vector_camera eps_a = this->m_model.A_initial(j)-this->m_model.A_parameters(j);
          
          error_total += .5  * transpose(eps_a) * inverse_cov * eps_a;
          
          epsilon_a(j) += transpose(C) * inverse_cov * eps_a;
          
        }
      }
      
    
      //  std::cout << "\n U Matrix array after camera positions is: " << U << "\n"; 


      // Add in the 3D point position constraint terms and
      // covariances. We only add constraints for Ground Control
      // Points (GCPs), not for 3D tie points.
      if (this->m_use_gcp_constraint) {
        for (unsigned i = 0; i < V.size(); ++i) {
          if ((*this->m_control_net)[i].type() == ControlPoint::GroundControlPoint) {
            
            matrix_point_point inverse_cov;
            inverse_cov = this->m_model.B_inverse_covariance(i);
            
            matrix_point_point D;
            D.set_identity();
            V(i) +=  transpose(D) * inverse_cov * D;
            
            vector_point eps_b = this->m_model.B_initial(i)-this->m_model.B_parameters(i);
            
            error_total += .5 * transpose(eps_b) * inverse_cov * eps_b;
            
            epsilon_b(i) += transpose(D) * inverse_cov * eps_b;
            
          }
        }
      }
      
      // flatten both epsilon_b and epsilon_a into a vector
      for (unsigned j = 0; j < U.size(); j++){
        subvector(g, current_g_length, num_cam_params) = static_cast<vector_camera>(epsilon_a(j));
        current_g_length += num_cam_params;
      }

      for (unsigned i = 0; i < V.size(); i++){
        subvector(g, current_g_length, num_pt_params) = static_cast<vector_point>(epsilon_b(i));
        current_g_length += num_pt_params;
      }

      //e at this point should be -g_a
      
      //std::cout << "\n (new) U Matrix array before lambda adjustment " << U << "\n";
     
      // "Augment" the diagonal entries of the U and V matrices with
      // the parameter lambda.
      {
        matrix_camera_camera u_lambda;
        u_lambda.set_identity();
        u_lambda *= this->m_lambda;

        for ( i = 0; i < U.size(); ++i ){
	  // static_cast<matrix_camera_camera>(U(i)) += u_lambda; 
	  U(i) += u_lambda;  
	}
      }
      //  std::cout << "\n (new) U Matrix array after lambda adjustment " << U << "\n";

      {
        matrix_point_point v_lambda;
        v_lambda.set_identity();
        v_lambda *= this->m_lambda;
        for ( i = 0; i < V.size(); ++i ){
          //static_cast<matrix_point_point>(V(i)) += v_lambda;
	  V(i) += v_lambda;
	}
      }
      
      // Create the 'e' vector in S * delta_a = e.  The first step is
      // to "flatten" our block structure to a vector that contains
      // scalar entries.
      Vector<double> e(this->m_model.num_cameras() * BundleAdjustModelT::camera_params_n);
      for (unsigned j = 0; j < epsilon_a.size(); ++j) {
        subvector(e, j*BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n) = 
          static_cast<vector_camera>(epsilon_a(j));
      }
               
      //Second Pass.  Compute Y and finish constructing e.
      i = 0; 
      for (typename ControlNetwork::const_iterator iter = this->m_control_net->begin(); 
           iter != this->m_control_net->end(); ++iter) { 
        for (typename ControlPoint::const_iterator measure_iter = (*iter).begin(); 
             measure_iter != (*iter).end(); ++measure_iter) { 
          unsigned j = measure_iter->image_id(); 
         
          // Compute the blocks of Y 
          Matrix<double> V_temp = static_cast<matrix_point_point>(V(i)); 
          chol_inverse(V_temp); 
          Y(j,i) = static_cast<matrix_camera_point>(W(j,i)) * transpose(V_temp) * V_temp; 

          // "Flatten the block structure to compute 'e'. 
          vector_camera temp = static_cast<matrix_camera_point>(Y(j,i))*static_cast<vector_point>(epsilon_b(i)); 
          subvector(e, j*num_cam_params, num_cam_params) -= temp;
        } 
        ++i; 
      } 
      
      // The S matrix is a m x m block matrix with blocks that are
      // camera_params_n x camera_params_n in size.  It has a sparse
      // skyline structure, which makes it more efficient to solve
      // through L*D*L^T decomposition and forward/back substitution
      // below.
      math::SparseSkylineMatrix<double> S(this->m_model.num_cameras()*num_cam_params, 
                                          this->m_model.num_cameras()*num_cam_params);
      
      i = 0;
      for (typename ControlNetwork::const_iterator iter = this->m_control_net->begin(); 
           iter != this->m_control_net->end(); ++iter) {
        for (typename ControlPoint::const_iterator j_measure_iter = (*iter).begin(); 
             j_measure_iter != (*iter).end(); ++j_measure_iter) {
          unsigned j = j_measure_iter->image_id();
          
          for (typename ControlPoint::const_iterator k_measure_iter = (*iter).begin(); 
               k_measure_iter != (*iter).end(); ++k_measure_iter) {
            unsigned k = k_measure_iter->image_id();
            
            // Compute the block entry...
            matrix_camera_camera temp = -static_cast< matrix_camera_point >(Y(j,i)) * 
              transpose( static_cast<matrix_camera_point>(W(k,i)) );
            // ... and "flatten" this matrix into the scalar entries of S
            for (unsigned aa = 0; aa < num_cam_params; ++aa) {
              for (unsigned bb = 0; bb < num_cam_params; ++bb) {
                // FIXME: This if clause is required at the moment to
                // ensure that we do not use the += on the symmetric
                // entries of the SparseSkylineMatrix.  These
                // symmetric entries are shallow, hence this code
                // would add the value twice if we're not careful
                // here.
                if (k*num_cam_params + bb <=
                    j*num_cam_params + aa) {
                  S(j*num_cam_params + aa,
                    k*num_cam_params + bb) += temp(aa,bb);

                  //  S_old(j*BundleAdjustModelT::camera_params_n + aa,
                  //  k*BundleAdjustModelT::camera_params_n + bb) += temp_old(aa,bb);
                  
                }
              }
            }
          }
        }
        ++i;
      }
      
      // Augment the diagonal entries S(i,i) with U(i)
      for (unsigned i = 0; i < this->m_model.num_cameras(); ++i) {
        // ... and "flatten" this matrix into the scalar entries of S
        for (unsigned aa = 0; aa < num_cam_params; ++aa) {
          for (unsigned bb = 0; bb < num_cam_params; ++bb) {
            // FIXME: This if clause is required at the moment to
            // ensure that we do not use the += on the symmetric
            // entries of the SparseSkylineMatrix.  These
            // symmetric entries are shallow, hence this code
            // would add the value twice if we're not careful
            // here.
            if (i*num_cam_params + bb <= 
                i*num_cam_params + aa) {
              S(i*num_cam_params + aa,
                i*num_cam_params + bb) += static_cast<matrix_camera_camera>(U(i))(aa,bb);
            }
          }
        }
      } 

      
      this->set_S(S);
      // Compute the LDL^T decomposition and solve using sparse methods.
      Vector<double> delta_a = sparse_solve(S, e);

      //std::cout << "delta_a: " << delta_a << "\n\n";

      // Save S; used for covariance calculations
     
      subvector(delta, current_delta_length, e.size()) = delta_a;
      current_delta_length += e.size();

      boost_sparse_vector<vector_point > delta_b(this->m_model.num_points());

      i = 0;
      for (typename ControlNetwork::const_iterator iter = this->m_control_net->begin(); 
           iter != this->m_control_net->end(); ++iter) {

        vector_point temp;

        for (typename ControlPoint::const_iterator j_measure_iter = (*iter).begin(); 
             j_measure_iter != (*iter).end(); ++j_measure_iter) {
          unsigned j = j_measure_iter->image_id();

          temp += transpose( static_cast<matrix_camera_point>(W(j,i)) ) * 
            subvector( delta_a, j*num_cam_params, num_cam_params );
        }

        Vector<double> delta_temp = static_cast<vector_point>(epsilon_b(i)) - temp;
        
        Matrix<double> hessian = static_cast<matrix_point_point>(V(i));
        
        solve(delta_temp, hessian);
        delta_b(i) = delta_temp;

        subvector(delta, current_delta_length, num_pt_params) = delta_temp;
        current_delta_length += num_pt_params;

        ++i;
      }
      
      dS = .5 * transpose(delta) *(this->m_lambda * delta + g);

      // -------------------------------
      // Compute the update error vector and predicted change
      // -------------------------------
      i = 0;
      double new_error_total = 0;
      for (typename ControlNetwork::const_iterator iter = this->m_control_net->begin(); 
           iter != this->m_control_net->end(); ++iter) {
        for (typename ControlPoint::const_iterator measure_iter = (*iter).begin(); 
             measure_iter != (*iter).end(); ++measure_iter) {

          unsigned j = measure_iter->image_id();
          
          // Compute error vector
          vector_camera new_a = this->m_model.A_parameters(j) + 
            subvector(delta_a, num_cam_params*j, num_cam_params);
          Vector<double> del_a = subvector(delta_a, num_cam_params*j, num_cam_params);

          vector_point new_b = this->m_model.B_parameters(i) + 
            static_cast<vector_point>(delta_b(i));

         
          // Apply robust cost function weighting
          Vector2 unweighted_error = measure_iter->dominant() - 
            this->m_model(i,j,new_a,new_b);
          double mag = norm_2(unweighted_error);
          double weight = sqrt(this->m_robust_cost_func(mag)) / mag;
          new_epsilon(i,j) = weight * unweighted_error;
          
          Matrix2x2 inverse_cov;
          Vector2 pixel_sigma = measure_iter->sigma();
          Vector2 epsilon_inst = new_epsilon(i,j);
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));
         
          new_error_total += .5 * transpose(epsilon_inst) * 
            inverse_cov * epsilon_inst;
        }
        ++i;
      }

      // Camera Constraints
      if (this->m_use_camera_constraint)
        for (unsigned j = 0; j < U.size(); ++j) {
          
          vector_camera new_a = this->m_model.A_parameters(j) + 
            subvector(delta_a, num_cam_params*j, num_cam_params);
          vector_camera eps_a = this->m_model.A_initial(j)-new_a;
          
          matrix_camera_camera inverse_cov;
          inverse_cov = this->m_model.A_inverse_covariance(j);
          
          new_error_total += .5 * transpose(eps_a) * inverse_cov * eps_a;
        }
      
      // GCP Error
      if (this->m_use_gcp_constraint)
        for (unsigned i = 0; i < V.size(); ++i) {
          if ((*this->m_control_net)[i].type() == ControlPoint::GroundControlPoint) {
            
            vector_point new_b = this->m_model.B_parameters(i) + 
              static_cast<vector_point>(delta_b(i));
            vector_point eps_b = this->m_model.B_initial(i)-new_b;
            
            matrix_point_point inverse_cov;
            inverse_cov = this->m_model.B_inverse_covariance(i);
            
            new_error_total += .5 * transpose(eps_b) * inverse_cov * eps_b;
          }
        }
      
      //Fletcher modification
      double Splus = new_error_total;     //Compute new objective
      double SS = error_total;            //Compute old objective
      double R = (SS - Splus)/dS;         // Compute ratio
      
      
  /*     std::cout << "New Objective: " << new_error_total << "\n";  */
/*       std::cout << "Old Objective: " << error_total << "\n";  */
/*       std::cout << "Lambda: " << this->m_lambda << "\n";  */
     


      if (R>0){

        for (unsigned j=0; j<this->m_model.num_cameras(); ++j)
          this->m_model.set_A_parameters(j, this->m_model.A_parameters(j) + 
                                         subvector(delta_a, num_cam_params*j,num_cam_params));
        for (unsigned i=0; i<this->m_model.num_points(); ++i)
          this->m_model.set_B_parameters(i, this->m_model.B_parameters(i) + 
                                         static_cast<vector_point>(delta_b(i)));
        
        // Summarize the stats from this step in the iteration
        //double overall_norm = sqrt(new_error_total);
        //double overall_delta = sqrt(error_total) - sqrt(new_error_total);

        abs_tol = vw::math::max(g) + vw::math::max(-g);
        rel_tol = transpose(delta)*delta;
        
        if(this->m_control == 0){
          double temp = 1 - pow((2*R - 1),3);
          if (temp < 1.0/3.0)
            temp = 1.0/3.0;
          
          this->m_lambda *= temp;
          this->m_nu = 2;

        } else if (this->m_control == 1){
          this->m_lambda /= 10;
        }

        return rel_tol;

      } else { // here we didn't make progress

        abs_tol = vw::math::max(g) + vw::math::max(-g);
        rel_tol = transpose(delta)*delta;

        if (this->m_control == 0){
          this->m_lambda *= this->m_nu;
          this->m_nu*=2;
        } else if (this->m_control == 1){
          this->m_lambda *= 10;
        }
        // double overall_delta = sqrt(error_total) - sqrt(new_error_total);

        /* // Taken over mostly by Bundle Adjust Report
        std::cout <<"\n" << "Sparse LM Iteration " << m_iterations << "  delta  "<< overall_delta << "  lambda: " << m_lambda << "   Ratio:  " << R <<"\n";
        */

        return ScalarTypeLimits<double>::highest();
      }

      return 0;
    }

  };

}} // namespace vw::camera

#endif//__VW_CAMERA_BUNDLE_ADJUSTMENT_SPARSE_H__
