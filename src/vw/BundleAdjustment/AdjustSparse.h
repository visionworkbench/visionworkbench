// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file AdjustSparse.h
///
/// Sparse implementation of bundle adjustment. Fast yo!

#ifndef __VW_BUNDLEADJUSTMENT_ADJUST_SPARSE_H__
#define __VW_BUNDLEADJUSTMENT_ADJUST_SPARSE_H__

// Vision Workbench
#include <vw/Math/MatrixSparseSkyline.h>
#include <vw/Core/Debugging.h>
#include <vw/BundleAdjustment/AdjustBase.h>
#include <vw/BundleAdjustment/CameraRelation.h>

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
  class AdjustSparse : public AdjustBase<BundleAdjustModelT, RobustCostT> {

    // Common Types
    typedef Matrix<double, 2, BundleAdjustModelT::camera_params_n> matrix_2_camera;
    typedef Matrix<double, 2, BundleAdjustModelT::point_params_n> matrix_2_point;
    typedef Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::camera_params_n> matrix_camera_camera;
    typedef Matrix<double,BundleAdjustModelT::point_params_n,BundleAdjustModelT::point_params_n> matrix_point_point;
    typedef Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::point_params_n> matrix_camera_point;
    typedef Vector<double,BundleAdjustModelT::camera_params_n> vector_camera;
    typedef Vector<double,BundleAdjustModelT::point_params_n> vector_point;

    math::MatrixSparseSkyline<double> m_S;
    std::vector<unsigned> m_ideal_ordering;
    Vector<unsigned> m_ideal_skyline;
    bool m_found_ideal_ordering;
    CameraRelationNetwork<JFeature> m_crn;
    typedef CameraNode<JFeature>::iterator crn_iter;

    // Reused structures
    std::vector< matrix_camera_camera > U;
    std::vector< matrix_point_point > V, V_inverse;
    std::vector< vector_camera > epsilon_a;
    std::vector< vector_point > epsilon_b;

  public:

    AdjustSparse( BundleAdjustModelT & model,
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
    // ___________________________________________________________
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
      math::MatrixSparseSkyline<double> S = this->S();  // Make copy as solve is destructive
      Matrix<double> Id(inverse_size, inverse_size);
      Id.set_identity();
      Matrix<double> Cov = multi_sparse_solve(S, Id);

      //pick out covariances of individual cameras
      for ( unsigned i = 0; i < num_cameras; i++ )
        sparse_cov(i) = submatrix(Cov, i*num_cam_params, i*num_cam_params, num_cam_params, num_cam_params);

      std::cout << "Covariance matrices for cameras are:"
                << sparse_cov << "\n\n";
    }

    // UPDATE IMPLEMENTATION
    //-------------------------------------------------------------
    // This is the sparse levenberg marquardt update step.  Returns
    // the average improvement in the cost function.
    double update(double &abs_tol, double &rel_tol) {
      ++this->m_iterations;
      Timer* time;

      VW_DEBUG_ASSERT(this->m_control_net->size() == this->m_model.num_points(), LogicErr() << "BundleAdjustment::update() : Number of bundles does not match the number of points in the bundle adjustment model.");

      // Reseting the values for U and V
      for ( uint32 j = 0; j < this->m_model.num_cameras(); j++ ) {
        U[j] = matrix_camera_camera();
        epsilon_a[j] = vector_camera();
      }
      for ( uint32 i = 0; i < this->m_model.num_points(); i++ ) {
        V[i] = matrix_point_point();
        epsilon_b[i] = vector_point();
      }

      unsigned num_cam_params = BundleAdjustModelT::camera_params_n;
      unsigned num_pt_params = BundleAdjustModelT::point_params_n;

      // Fletcher LM parameteres
      double dS = 0; //Predicted improvement for Fletcher modification

      // Populate the Jacobian, which is broken into two sparse
      // matrices A & B, as well as the error matrix and the W
      // matrix.
      time = new Timer("Solve for Image Error, Jacobian, U, V, and W:", DebugMessage, "ba");

      double error_total = 0; // assume this is r^T\Sigma^{-1}r
      for ( uint32 j = 0; j < m_crn.size(); j++ ) {
        for ( crn_iter fiter = m_crn[j].begin();
              fiter != m_crn[j].end(); fiter++ ) {

          matrix_2_camera A = this->m_model.A_jacobian( (**fiter).m_point_id, j,
                                                        this->m_model.A_parameters(j),
                                                        this->m_model.B_parameters((**fiter).m_point_id) );
          matrix_2_point B = this->m_model.B_jacobian( (**fiter).m_point_id, j,
                                                       this->m_model.A_parameters(j),
                                                       this->m_model.B_parameters((**fiter).m_point_id) );
          // Apply robust cost function weighting
          Vector2 error = (**fiter).m_location -
            this->m_model((**fiter).m_point_id,j,this->m_model.A_parameters(j),
                          this->m_model.B_parameters((**fiter).m_point_id));

          double mag = norm_2(error);
          double weight = sqrt(this->m_robust_cost_func(mag)) / mag;
          error *= weight;

          Matrix2x2 inverse_cov;
          Vector2 pixel_sigma = (**fiter).m_scale;
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));
          error_total += .5 * transpose(error) *
            inverse_cov * error;

          // Storing intermediate values
          U[j] += transpose(A) * inverse_cov * A;
          V[(**fiter).m_point_id] += transpose(B) * inverse_cov * B;
          epsilon_a[j] += transpose(A) * inverse_cov * error;
          epsilon_b[(**fiter).m_point_id] += transpose(B) * inverse_cov * error;
          (**fiter).m_w = transpose(A) * inverse_cov * B;
        }
      }
      delete time;

      // set initial lambda, and ignore if the user has touched it
      if ( this->m_iterations == 1 && this->m_lambda == 1e-3 ) {
        time = new Timer("Solving for Lambda:", DebugMessage, "ba");
        double max = 0.0;
        for (unsigned i = 0; i < U.size(); ++i)
          for (unsigned j = 0; j < BundleAdjustModelT::camera_params_n; ++j){
            if (fabs(U[i](j,j)) > max)
              max = fabs(U[i](j,j));
          }
        for (unsigned i = 0; i < V.size(); ++i)
          for (unsigned j = 0; j < BundleAdjustModelT::point_params_n; ++j) {
            if ( fabs(V[i](j,j)) > max)
              max = fabs(V[i](j,j));
          }
        this->m_lambda = max * 1e-10;
        delete time;
      }

      // Add in the camera position and pose constraint terms and covariances.
      time = new Timer("Solving for Camera and GCP error:",DebugMessage,"ba");
      if ( this->m_use_camera_constraint )
        for ( unsigned j = 0; j < U.size(); ++j ) {
          matrix_camera_camera inverse_cov;
          inverse_cov = this->m_model.A_inverse_covariance(j);
          U[j] += inverse_cov;
          vector_camera eps_a = this->m_model.A_initial(j)-this->m_model.A_parameters(j);
          error_total += .5  * transpose(eps_a) * inverse_cov * eps_a;
          epsilon_a[j] += inverse_cov * eps_a;
        }

      // Add in the 3D point position constraint terms and
      // covariances. We only add constraints for Ground Control
      // Points (GCPs), not for 3D tie points.
      if (this->m_use_gcp_constraint)
        for ( unsigned i = 0; i < V.size(); ++i )
          if ((*this->m_control_net)[i].type() == ControlPoint::GroundControlPoint) {
            matrix_point_point inverse_cov;
            inverse_cov = this->m_model.B_inverse_covariance(i);
            V[i] += inverse_cov;
            vector_point eps_b = this->m_model.B_initial(i)-this->m_model.B_parameters(i);
            error_total += .5 * transpose(eps_b) * inverse_cov * eps_b;
            epsilon_b[i] += inverse_cov * eps_b;
          }
      delete time;

      // flatten both epsilon_b and epsilon_a into a vector
      time = new Timer("Augmenting with lambda",DebugMessage,"ba");
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
      delete time;

      // Create the 'e' vector in S * delta_a = e.  The first step is
      // to "flatten" our block structure to a vector that contains
      // scalar entries.
      time = new Timer("Create special e vector", DebugMessage, "ba");
      Vector<double> e(this->m_model.num_cameras() * BundleAdjustModelT::camera_params_n);
      for (unsigned j = 0; j < epsilon_a.size(); ++j) {
        subvector(e, j*BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n) =
          epsilon_a[j];
      }

      // Compute V inverse
      for ( uint32 i = 0; i < this->m_model.num_points(); i++ ) {
        Matrix<double> V_temp = V[i];
        chol_inverse( V_temp );
        V_inverse[i] = V_temp;
      }

      // Compute Y and finish constructing e.
      for ( uint32 j = 0; j < m_crn.size(); j++ ) {
        for ( crn_iter fiter = m_crn[j].begin();
              fiter != m_crn[j].end(); fiter++ ) {
          // Compute the blocks of Y
          (**fiter).m_y = (**fiter).m_w * transpose(V_inverse[(**fiter).m_point_id])
            * V_inverse[ (**fiter).m_point_id ];
          // Flatten the block structure to compute 'e'
          subvector(e, j*num_cam_params, num_cam_params) -= (**fiter).m_y
            * epsilon_b[ (**fiter).m_point_id ];
        }
      }

      delete time;

      // --- BUILD SPARSE, SOLVE A'S UPDATE STEP -------------------------
      time = new Timer("Build Sparse", DebugMessage, "ba");

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

        { // Filling in off diagonal
          for ( uint32 k = j+1; k < m_crn.size(); k++ ) {
            matrix_camera_camera S_jk;

            // Iterate across all features seens by this camera and k
            for ( crn_iter f_j_iter = m_crn[j].begin();
                  f_j_iter != m_crn[j].end(); f_j_iter++ ) {
              for ( crn_iter f_k_iter = (**f_j_iter).m_connections.begin();
                    f_k_iter != (**f_j_iter).m_connections.end(); f_k_iter++ ) {
                if ( (**f_k_iter).m_camera_id == k ) {
                  S_jk -= (**f_j_iter).m_y
                    * transpose((**f_k_iter).m_w);
                  continue;
                }
              }
            } // done iteration through features

            // Loading into sparse matrix
            // - if it seems we are loading in oddly, it's because the sparse
            //   matrix is row major.
            if ( S_jk != matrix_camera_camera() ) {
              submatrix( S, k*num_cam_params, j*num_cam_params,
                         num_cam_params, num_cam_params ) = transpose(S_jk);
            }
          }
        }
      }

      m_S = S; // S is modified in sparse solve. Keeping a copy.
      delete time;

      // Computing ideal ordering
      if (!m_found_ideal_ordering) {
        time = new Timer("Solving Cuthill-Mckee", DebugMessage, "ba");
        m_ideal_ordering = cuthill_mckee_ordering(S,num_cam_params);
        math::MatrixReorganize<math::MatrixSparseSkyline<double> > mod_S( S, m_ideal_ordering );
        m_ideal_skyline = solve_for_skyline(mod_S);

        m_found_ideal_ordering = true;
        delete time;
      }

      time = new Timer("Solve Delta A", DebugMessage, "ba");

      // Compute the LDL^T decomposition and solve using sparse methods.
      math::MatrixReorganize<math::MatrixSparseSkyline<double> > modified_S( S, m_ideal_ordering );
      Vector<double> delta_a = sparse_solve( modified_S,
                                             reorganize(e, m_ideal_ordering),
                                             m_ideal_skyline );
      delta_a = reorganize(delta_a, modified_S.inverse());
      delete time;

      // --- SOLVE B'S UPDATE STEP ---------------------------------

      // Back Solving for Delta B
      time = new Timer("Solve Delta B", DebugMessage, "ba");
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
      delete time;

      dS = 0;
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
      time = new Timer("Solve for Updated Error", DebugMessage, "ba");
      double new_error_total = 0;
      for ( uint32 j = 0; j < m_crn.size(); j++ ) {
        for ( crn_iter fiter = m_crn[j].begin();
              fiter != m_crn[j].end(); fiter++ ) {
          // Compute error vector
          vector_camera new_a = this->m_model.A_parameters(j) +
            subvector( delta_a, num_cam_params*j, num_cam_params );
          vector_point new_b = this->m_model.B_parameters((**fiter).m_point_id) +
            subvector( delta_b, num_pt_params*(**fiter).m_point_id, num_pt_params );

          // Apply robust cost function weighting
          Vector2 error = (**fiter).m_location -
            this->m_model((**fiter).m_point_id,j,new_a,new_b);
          double mag = norm_2( error );
          double weight = sqrt( this->m_robust_cost_func(mag)) / mag;
          error *= weight;

          Matrix2x2 inverse_cov;
          Vector2 pixel_sigma = (**fiter).m_scale;
          inverse_cov(0,0) = 1/(pixel_sigma(0)*pixel_sigma(0));
          inverse_cov(1,1) = 1/(pixel_sigma(1)*pixel_sigma(1));

          new_error_total += .5 * transpose(error) *
            inverse_cov * error;
        }
      }

      // Camera Constraints
      if ( this->m_use_camera_constraint )
        for (unsigned j = 0; j < U.size(); ++j) {

          vector_camera new_a = this->m_model.A_parameters(j) +
            subvector(delta_a, num_cam_params*j, num_cam_params);
          vector_camera eps_a = this->m_model.A_initial(j)-new_a;

          matrix_camera_camera inverse_cov;
          inverse_cov = this->m_model.A_inverse_covariance(j);
          new_error_total += .5 * transpose(eps_a) * inverse_cov * eps_a;
        }

      // GCP Error
      if ( this->m_use_gcp_constraint )
        for ( unsigned i = 0; i < V.size(); ++i )
          if ( (*this->m_control_net)[i].type() ==
               ControlPoint::GroundControlPoint) {

            vector_point new_b = this->m_model.B_parameters(i) +
              subvector( delta_b, num_pt_params*i, num_pt_params );
            vector_point eps_b = this->m_model.B_initial(i)-new_b;
            matrix_point_point inverse_cov;
            inverse_cov = this->m_model.B_inverse_covariance(i);
            new_error_total += .5 * transpose(eps_b) * inverse_cov * eps_b;
          }
      delete time;

      //Fletcher modification
      double Splus = new_error_total;     //Compute new objective
      double SS = error_total;            //Compute old objective
      double R = (SS - Splus)/dS;         // Compute ratio

      {
        double gvec_max = -1e30, gvec_nmax = -1e30;
        for ( uint32 j = 0; j < this->m_model.num_cameras(); j++ ) {
          gvec_max = std::max( gvec_max, vw::math::max(epsilon_a[j]) );
          gvec_nmax = std::max( gvec_nmax, vw::math::max(-epsilon_a[j]) );
        }
        for ( uint32 i = 0; i < this->m_model.num_points(); i++ ) {
          gvec_max = std::max( gvec_max, vw::math::max(epsilon_b[i]) );
          gvec_nmax = std::max( gvec_nmax, vw::math::max(-epsilon_b[i]) );
        }
        abs_tol = gvec_max + gvec_nmax;
      }
      rel_tol = transpose(delta_a)*delta_a + transpose(delta_b)*delta_b;

      if ( R > 0 ) {

        time = new Timer("Setting Parameters",DebugMessage,"ba");
        for (unsigned j = 0; j < this->m_model.num_cameras(); ++j)
          this->m_model.set_A_parameters(j, this->m_model.A_parameters(j) +
                                         subvector(delta_a, num_cam_params*j,num_cam_params));
        for (unsigned i = 0; i < this->m_model.num_points(); ++i)
          this->m_model.set_B_parameters(i, this->m_model.B_parameters(i) +
                                         subvector(delta_b, num_pt_params*i,num_pt_params));
        delete time;

        if ( this->m_control == 0 ) {
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
      if ( this->m_control == 0 ) {
        this->m_lambda *= this->m_nu;
        this->m_nu*=2;
      } else if ( this->m_control == 1 )
        this->m_lambda *= 10;

      return 0;
    }

  };

}} // namespace vw::ba

#endif//__VW_BUNDLEADJUSTMENT_ADJUST_SPARSE_H__
