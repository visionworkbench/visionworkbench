// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file AdjustRobustSparse.h
///
/// Robust Sparse implementation of bundle adjustment. Fast yo!

#ifndef __VW_BUNDLEADJUSTMENT_ADJUST_ROBUST_SPARSE_H__
#define __VW_BUNDLEADJUSTMENT_ADJUST_ROBUST_SPARSE_H__

// Vision Workbench
#include <vw/BundleAdjustment/AdjustBase.h>
#include <vw/Math/MatrixSparseSkyline.h>

// Boost
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
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
namespace ba {

  template <class BundleAdjustModelT, class RobustCostT>
  class AdjustRobustSparse : public AdjustBase<BundleAdjustModelT, RobustCostT> {

    // Common Types
    typedef Matrix<double, 2, BundleAdjustModelT::camera_params_n> matrix_2_camera;
    typedef Matrix<double, 2, BundleAdjustModelT::point_params_n> matrix_2_point;
    typedef Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::camera_params_n> matrix_camera_camera;
    typedef Matrix<double,BundleAdjustModelT::point_params_n,BundleAdjustModelT::point_params_n> matrix_point_point;
    typedef Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::point_params_n> matrix_camera_point;
    typedef Vector<double,BundleAdjustModelT::camera_params_n> vector_camera;
    typedef Vector<double,BundleAdjustModelT::point_params_n> vector_point;

    math::MatrixSparseSkyline<double> m_S;
    std::vector<size_t> m_ideal_ordering;
    Vector<size_t> m_ideal_skyline;
    bool m_found_ideal_ordering;
    CameraRelationNetwork<JFeature> m_crn;
    typedef CameraNode<JFeature>::iterator crn_iter;

    // Reused structures
    std::vector< matrix_camera_camera > U;
    std::vector< matrix_point_point > V, V_inverse;
    std::vector< vector_camera > epsilon_a;
    std::vector< vector_point > epsilon_b;

  public:

    AdjustRobustSparse( BundleAdjustModelT & model,
                        RobustCostT const& robust_cost_func,
                        bool use_camera_constraint=true,
                        bool use_gcp_constraint=true) :
    AdjustBase<BundleAdjustModelT,RobustCostT>( model, robust_cost_func,
                                                use_camera_constraint,
                                                use_gcp_constraint ),
      U( this->m_model.num_cameras() ), V( this->m_model.num_points() ),
      V_inverse( this->m_model.num_points() ),
      epsilon_a( this->m_model.num_cameras() ), epsilon_b( this->m_model.num_points() ) {
      m_crn.read_controlnetwork( *(this->m_control_net).get() );
      m_found_ideal_ordering = false;
    }

    math::MatrixSparseSkyline<double> S() const { return m_S; }

    // Covariance Calculator
    // __________________________________________________
    // This routine inverts a sparse matrix S, and prints the individual
    // covariance matrices for each camera
    void covCalc() {
      // camera params
      size_t num_cam_params = BundleAdjustModelT::camera_params_n;
      size_t num_cameras = this->m_model.num_cameras();

      size_t inverse_size = num_cam_params * num_cameras;

      typedef Matrix<double, BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n> matrix_camera_camera;

      // final vector of camera covariance matrices
      vw::Vector< matrix_camera_camera > sparse_cov(num_cameras);

      // Get the S matrix from the model
      math::MatrixSparseSkyline<double> S = this->S();
      Matrix<double> Id(inverse_size, inverse_size);
      Id.set_identity();
      Matrix<double> Cov = multi_sparse_solve(S, Id);

      //pick out covariances of indvidual cameras
      for ( size_t i = 0; i < num_cameras; i++ )
        sparse_cov(i) = submatrix(Cov, i*num_cam_params,
                                  i*num_cam_params,
                                  num_cam_params,
                                  num_cam_params);

      std::cout << "Covariance matrices for cameras are:"
                << sparse_cov << "\n\n";
    }

    // UPDATE IMPLEMENTATION
    //-------------------------------------------------------------
    // This is the sparse levenberg marquardt update step.  Returns
    // the average improvement in the cost function.
    double update(double &abs_tol, double &rel_tol) {
      ++this->m_iterations;
      boost::scoped_ptr<Timer> time;

      VW_DEBUG_ASSERT(this->m_control_net->size() == this->m_model.num_points(), LogicErr() << "BundleAdjustment::update() : Number of bundles does not match the number of points in the bundle adjustment model.");

      // Reseting the values for U, V, and epsilon.
      for ( uint32 j = 0; j < this->m_model.num_cameras(); j++ ) {
        U[j] = matrix_camera_camera();
        epsilon_a[j] = vector_camera();
      }
      for ( uint32 i = 0; i < this->m_model.num_points(); i++ ) {
        V[i] = matrix_point_point();
        epsilon_b[i] = vector_point();
      }

      size_t num_cam_params = BundleAdjustModelT::camera_params_n;
      size_t num_pt_params = BundleAdjustModelT::point_params_n;

      // Degrees of freedom for data (can be modified later)
      double t_df = 4;
      // dimension of pixels
      double t_dim_pixel = 2;
      // dimension of cameras
      double t_dim_cam   = 6;
      // dimension of world points
      double t_dim_pt    = 3;

      // Populate the Jacobian, which is broken into two sparse
      // matrices A & B, as well as the error matrix and the W
      // matrix.
      time.reset(new Timer("Solve for Image Error, Jacobian, U, V, and W:", DebugMessage, "ba"));
      double robust_objective = 0.0;
      for ( uint32 j = 0; j < m_crn.size(); j++ ) {
        for ( crn_iter fiter = m_crn[j].begin();
              fiter != m_crn[j].end(); fiter++ ) {
          uint32 i = (**fiter).m_point_id;

          matrix_2_camera A = this->m_model.A_jacobian(i,j,this->m_model.A_parameters(j),
                                                       this->m_model.B_parameters(i));
          matrix_2_point B = this->m_model.B_jacobian(i,j,this->m_model.A_parameters(j),
                                                      this->m_model.B_parameters(i));

          // Apply robust cost function weighting
          Vector2 unweighted_error;
          try {
            unweighted_error = (**fiter).m_location -
              this->m_model(i,j,this->m_model.A_parameters(j),
                            this->m_model.B_parameters(i));
          } catch ( camera::PixelToRayErr &e ) {}

          Vector2 pixel_sigma = (**fiter).m_scale;
          Matrix2x2 inverse_cov;
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));

          double S_weight = transpose(unweighted_error) * inverse_cov * unweighted_error;
          double mu_weight = (t_df + t_dim_pixel)/(t_df + S_weight);

          // do NOT want epsilon_inst scaled
          robust_objective += 0.5*(t_df + t_dim_pixel)*log(1 + S_weight/t_df);

          // Storing intermediate values
          U[j] += mu_weight * transpose(A) * inverse_cov * A;
          V[i] += mu_weight * transpose(B) * inverse_cov * B;
          (**fiter).m_w = mu_weight * transpose(A) * inverse_cov * B;
          epsilon_a[j] += mu_weight * transpose(A) * inverse_cov * unweighted_error;
          epsilon_b[i] += mu_weight * transpose(B) * inverse_cov * unweighted_error;
        }
      }
      time.reset();

      // Add in the camera position and pose constraint terms and covariances.
      time.reset(new Timer("Solving for Camera and GCP error:",DebugMessage,"ba"));
      if ( this->m_use_camera_constraint )
        for (size_t j = 0; j < U.size(); ++j) {
          matrix_camera_camera inverse_cov;
          inverse_cov = this->m_model.A_inverse_covariance(j);
          vector_camera eps_a = this->m_model.A_target(j)-this->m_model.A_parameters(j);
          double S_weight = transpose(eps_a) * inverse_cov * eps_a;
          double mu_weight = (t_df + t_dim_cam)/(t_df + S_weight);
          U[j] += mu_weight * inverse_cov;
          robust_objective += 0.5*(t_df + t_dim_cam)*log(1 + S_weight/t_df);
          epsilon_a[j] += mu_weight * inverse_cov * eps_a;
        }

      // Add in the 3D point position constraint terms and
      // covariances. We only add constraints for Ground Control
      // Points (GCPs), not for 3D tie points.
      if ( this->m_use_gcp_constraint )
        for (size_t i = 0; i < V.size(); ++i) {
          if ((*this->m_control_net)[i].type() == ControlPoint::GroundControlPoint) {
            matrix_point_point inverse_cov;
            inverse_cov = this->m_model.B_inverse_covariance(i);
            vector_point eps_b = this->m_model.B_target(i)-this->m_model.B_parameters(i);
            double S_weight = transpose(eps_b) * inverse_cov * eps_b;
            double mu_weight = (t_df + t_dim_pt)/(t_df + S_weight);
            V[i] +=  mu_weight * inverse_cov;
            robust_objective += 0.5*(t_df + t_dim_pt)*log(1 + S_weight/t_df);
            epsilon_b[i] += mu_weight * inverse_cov * eps_b;
          }
        }
      time.reset();

      // set initial lambda, and ignore if the user has touched it
      if (this->m_iterations == 1 && this->m_lambda == 1e-3){
        time.reset(new Timer("Solving for Lambda:",DebugMessage,"ba"));
        double max = 0.0;
        for (size_t i = 0; i < U.size(); ++i)
          for (size_t j = 0; j < BundleAdjustModelT::camera_params_n; ++j){
            if (fabs(U[i](j,j)) > max)
              max = fabs(U[i](j,j));
          }
        for (size_t i = 0; i < V.size(); ++i)
          for (size_t j = 0; j < BundleAdjustModelT::point_params_n; ++j) {
            if ( fabs(V[i](j,j)) > max)
              max = fabs(V[i](j,j));
          }
        this->m_lambda = max * 1e-10;
        time.reset();
      }

      time.reset(new Timer("Augmenting with Lambda",DebugMessage,"ba"));
      //e at this point should be -g_a

      // "Augment" the diagonal entries of the U and V matrices with
      // the parameter lambda.
      {
        matrix_camera_camera u_lambda;
        u_lambda.set_identity();
        u_lambda *= this->m_lambda;
        for ( uint32 i = 0; i < U.size(); ++i )
          U[i] += u_lambda;
      }
      {
        matrix_point_point v_lambda;
        v_lambda.set_identity();
        v_lambda *= this->m_lambda;
        for ( uint32 i = 0; i < V.size(); ++i )
          V[i] += v_lambda;
      }
      time.reset();

      // Create the 'e' vector in S * delta_a = e.  The first step is
      // to "flatten" our block structure to a vector that contains
      // scalar entries.
      time.reset(new Timer("Create special e vector", DebugMessage, "ba"));
      Vector<double> e(this->m_model.num_cameras() * BundleAdjustModelT::camera_params_n);
      for ( size_t j = 0; j < epsilon_a.size(); ++j ) {
        subvector(e, j*BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n) =
          epsilon_a[j];
      }

      // Compute V inverse
      for ( uint32 i = 0; i < this->m_model.num_points(); i++ ) {
        Matrix<double> V_temp = V[i];
        chol_inverse( V_temp );
        V_inverse[i] = transpose(V_temp)*V_temp;
      }

      // Compute Y and finish constructing e.
      for ( uint32 j = 0; j < m_crn.size(); j++ ) {
        for ( crn_iter fiter = m_crn[j].begin();
              fiter != m_crn[j].end(); fiter++ ) {
          uint32 i = (**fiter).m_point_id;
          // Compute the blocks of Y
          (**fiter).m_y = (**fiter).m_w * V_inverse[i];
          // Flatten the block structure to compute 'e'
          subvector(e, j*num_cam_params, num_cam_params) -= (**fiter).m_y
            * epsilon_b[i];
        }
      }
      time.reset();

      // --- BUILD SPARSE, SOLVE A'S UPDATE STEP -------------------------
      time.reset(new Timer("Build Sparse", DebugMessage, "ba"));

      // The S matrix is a m x m block matrix with blocks that are
      // camera_params_n x camera_params_n in size.  It has a sparse
      // skyline structure, which makes it more efficient to solve
      // through L*D*L^T decomposition and forward/back substitution
      // below.
      math::MatrixSparseSkyline<double> S(this->m_model.num_cameras()*num_cam_params,
                                          this->m_model.num_cameras()*num_cam_params);
      for ( uint32 j = 0; j < m_crn.size(); j++ ) {
        { // Filling in diagonal
          matrix_camera_camera S_jj;

          // Iterate across all features seen by the camera
          for ( crn_iter fiter = m_crn[j].begin();
                fiter != m_crn[j].end(); fiter++ ) {
            S_jj -= (**fiter).m_y*transpose((**fiter).m_w);
          }

          // Augmenting Diagonal
          S_jj += U[j];

          // Loading into sparse matrix
          uint32 offset = j * num_cam_params;
          for ( uint32 aa = 0; aa < num_cam_params; aa++ ) {
            for ( uint32 bb = aa; bb < num_cam_params; bb++ ) {
              S( offset+bb, offset+aa ) = S_jj(aa,bb);  // Transposing
            }
          }
        }

        // Filling in off diagonal
        for ( uint32 k = j+1; k < m_crn.size(); k++ ) {
          typedef boost::shared_ptr<JFeature> f_ptr;
          typedef std::multimap< uint32, f_ptr >::iterator mm_iterator;
          std::pair< mm_iterator, mm_iterator > feature_range;
          feature_range = m_crn[j].map.equal_range( k );

          // Iterating through all features in camera j that have
          // connections to camera k.
          matrix_camera_camera S_jk;
          bool found = false;
          for ( mm_iterator f_j_iter = feature_range.first;
                f_j_iter != feature_range.second; f_j_iter++ ) {
            f_ptr f_k = (*f_j_iter).second->m_map[k];
            found = true;
            S_jk -= (*f_j_iter).second->m_y *
              transpose( (*f_k).m_w );
          }

          // Loading into sparse matrix
          // - if it seems we are loading in oddly, it's because the sparse
          //   matrix is row major.
          if ( found ) {
            submatrix( S, k*num_cam_params, j*num_cam_params,
                       num_cam_params, num_cam_params ) = transpose(S_jk);
          }
        }
      }

      m_S = S; // S is modified in sparse solve. Keeping a copy;
      time.reset();

      // Computing ideal ordering of sparse matrix
      if ( !m_found_ideal_ordering ) {
        time.reset(new Timer("Solving Cuthill-Mckee", DebugMessage, "ba"));
        m_ideal_ordering = cuthill_mckee_ordering(S,num_cam_params);
        math::MatrixReorganize<math::MatrixSparseSkyline<double> > mod_S( S, m_ideal_ordering );
        m_ideal_skyline = solve_for_skyline( mod_S );

        m_found_ideal_ordering = true;
        time.reset();
      }

      time.reset(new Timer("Solve Delta A", DebugMessage, "ba"));

      // Compute the LDL^T decomposition and solve using sparse methods.
      math::MatrixReorganize<math::MatrixSparseSkyline<double> > modified_S( S, m_ideal_ordering );
      Vector<double> delta_a = sparse_solve( modified_S,
                                             reorganize(e, m_ideal_ordering),
                                             m_ideal_skyline );
      delta_a = reorganize( delta_a, modified_S.inverse() );
      time.reset();

      // --- SOLVE B'S UPDATE STEP ---------------------------------

      // Back Solving for Delta B
      time.reset(new Timer("Solve Delta B", DebugMessage, "ba"));
      Vector<double> delta_b( this->m_model.num_points() * num_pt_params );

      {
        // delta_b = inverse(V)*( epsilon_b - sum_across_cam( WijT * delta_aj ) )

        // Building right half, sum( WijT * delta_aj )
        std::vector< vector_point > right_delta_b( this->m_model.num_points() );
        for ( uint32 j = 0; j < m_crn.size(); j++ ) {
          for ( crn_iter fiter = m_crn[j].begin();
                fiter != m_crn[j].end(); fiter++ ) {
            right_delta_b[ (**fiter).m_point_id ] += transpose( (**fiter).m_w ) *
              subvector( delta_a, j*num_cam_params, num_cam_params );
          }
        }

        // Solving for delta b
        for ( uint32 i = 0; i < this->m_model.num_points(); i++ ) {
          Vector<double> delta_temp = epsilon_b[i] - right_delta_b[i];
          Matrix<double> hessian = V[i];
          solve( delta_temp, hessian );
          subvector( delta_b, i*num_pt_params, num_pt_params ) = delta_temp;
        }
      }
      time.reset();

      //Predicted improvement for Fletcher modification
      double dS = 0;
      for ( uint32 j = 0; j < this->m_model.num_cameras(); j++ )
        dS += transpose(subvector(delta_a,j*num_cam_params,num_cam_params))
          * ( this->m_lambda * subvector(delta_a,j*num_cam_params,num_cam_params) +
              epsilon_a[j] );
      for ( uint32 i = 0; i < this->m_model.num_points(); i++ )
        dS += transpose(subvector(delta_b,i*num_pt_params,num_pt_params))
          * ( this->m_lambda * subvector(delta_b,i*num_pt_params,num_pt_params) +
              epsilon_b[i] );
      dS *= 0.5;

      // -------------------------------
      // Compute the update error vector and predicted change
      // -------------------------------
      time.reset(new Timer("Solve for Updated Error", DebugMessage, "ba"));
      double new_robust_objective = 0;
      for ( uint32 j = 0; j < m_crn.size(); j++ ) {
        for ( crn_iter fiter = m_crn[j].begin();
              fiter != m_crn[j].end(); fiter++ ) {
          uint32 i = (**fiter).m_point_id;

          // Computer error vector
          vector_camera new_a = this->m_model.A_parameters(j) +
            subvector( delta_a, num_cam_params*j, num_cam_params );
          vector_point new_b = this->m_model.B_parameters(i) +
            subvector( delta_b, num_pt_params*i, num_pt_params );

          // Apply robust cost function weighting
          Vector2 unweighted_error;
          try {
            unweighted_error = (**fiter).m_location -
              this->m_model(i,j,new_a,new_b);
          } catch ( camera::PixelToRayErr &e ) {}

          Matrix2x2 inverse_cov;
          Vector2 pixel_sigma = (**fiter).m_scale;
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));

          // Populate the S_weights, mu_weights vectors
          double S_weight = transpose(unweighted_error) * inverse_cov * unweighted_error;

          new_robust_objective += 0.5*(t_df + t_dim_pixel)*log(1 + S_weight/t_df);
        }
      }

      // Camera Constraints
      if ( this->m_use_camera_constraint )
        for ( size_t j = 0; j < U.size(); ++j ) {
          // note the signs here: should be +
          vector_camera new_a = this->m_model.A_parameters(j) +
            subvector(delta_a, num_cam_params*j, num_cam_params);
          vector_camera eps_a = this->m_model.A_target(j)-new_a;
          matrix_camera_camera inverse_cov;
          inverse_cov = this->m_model.A_inverse_covariance(j);
          double S_weight = transpose(eps_a) * inverse_cov * eps_a;
          new_robust_objective += 0.5*(t_df + t_dim_cam)*log(1 + S_weight/t_df);
        }

      // GCP Error
      if ( this->m_use_gcp_constraint )
        for ( size_t i = 0; i < V.size(); ++i )
          if ((*this->m_control_net)[i].type() == ControlPoint::GroundControlPoint) {
            // note the signs here: should be +
            vector_point new_b = this->m_model.B_parameters(i) +
              subvector( delta_b, num_pt_params*i, num_pt_params );
            vector_point eps_b = this->m_model.B_target(i)-new_b;
            matrix_point_point inverse_cov;
            inverse_cov = this->m_model.B_inverse_covariance(i);
            double S_weight = transpose(eps_b)*inverse_cov*eps_b;
            new_robust_objective += 0.5*(t_df + t_dim_pt)*log(1 + S_weight/t_df);
          }
      time.reset();

      //Fletcher modification
      double Splus = new_robust_objective;     //Compute new objective
      double SS = robust_objective;            //Compute old objective
      double R = (SS - Splus)/dS;         // Compute ratio

      rel_tol = -1e30;
      BOOST_FOREACH( vector_camera const& a, epsilon_a )
        rel_tol = std::max( rel_tol, vw::math::max( abs( a ) ) );
      BOOST_FOREACH( vector_point const& b, epsilon_b )
        rel_tol = std::max( rel_tol, vw::math::max( abs( b ) ) );
      abs_tol = Splus;

      if ( R > 0 ) {

        time.reset(new Timer("Setting Parameters",DebugMessage,"ba"));
        for (size_t j=0; j<this->m_model.num_cameras(); ++j)
          this->m_model.set_A_parameters(j, this->m_model.A_parameters(j) +
                                         subvector(delta_a, num_cam_params*j,num_cam_params));
        for (size_t i=0; i<this->m_model.num_points(); ++i)
          this->m_model.set_B_parameters(i, this->m_model.B_parameters(i) +
                                         subvector(delta_b, num_pt_params*i,num_pt_params));
        time.reset();

        if(this->m_control == 0){
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

}} // namespace vw::ba

#endif//__VW_BUNDLEADJUSTMENT_ADJUST_ROBUST_SPARSE_H__
