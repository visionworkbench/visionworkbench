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

/// \file BundleAdjust.h
/// 
/// Optimization classes for carrying out bundle adjustment of many
/// camera images.

#ifndef __VW_CAMERA_BUNDLE_ADJUST_H__
#define __VW_CAMERA_BUNDLE_ADJUST_H__

// Vision Workbench
#include <vw/Camera/ControlNetwork.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Core/Debugging.h>

// Boost 
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace vw {
namespace camera {
  
  // CRTP Base class for Bundle Adjustment functors.
  // 
  // The child class must implement this method:
  //
  //   Vector2 operator() ( Vector<double, CameraParamsN> const& a_j, Vector<double, PointParamsN> const& b_i ) const;
  //
  template <class ImplT, unsigned CameraParamsN, unsigned PointParamsN>
  struct BundleAdjustmentModelBase {

    static const unsigned camera_params_n = CameraParamsN;
    static const unsigned point_params_n = PointParamsN;

    /// \cond INTERNAL
    // Methods to access the derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
    /// \endcond

    virtual ~BundleAdjustmentModelBase() {}
   

    // Approximate the jacobian for small variations in the a_j
    // parameters (camera parameters). 
    inline Matrix<double, 2, CameraParamsN> A_jacobian ( unsigned i, unsigned j,
                                                         Vector<double, CameraParamsN> const& a_j, 
                                                         Vector<double, PointParamsN> const& b_i ) {
      // Get nominal function value
      Vector2 h0 = impl()(i,j,a_j,b_i);

      // Jacobian is #outputs x #params
      Matrix<double, 2, CameraParamsN> J;

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( unsigned i=0; i < CameraParamsN; ++i ){
        Vector<double, CameraParamsN> a_j_prime = a_j;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(a_j(i)*1e-7);
        a_j_prime(i) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        Vector2 hi = impl()(i,j,a_j_prime, b_i);
        select_col(J,i) = (hi-h0)/epsilon;
      }
      return J;
    }

    // Approximate the jacobian for small variations in the b_i
    // parameters (3d point locations). 
    inline Matrix<double, 2, PointParamsN> B_jacobian ( unsigned i, unsigned j,
                                                        Vector<double, CameraParamsN> const& a_j, 
                                                        Vector<double, PointParamsN> const& b_i ) {
      // Get nominal function value
      Vector2 h0 = impl()(i,j,a_j, b_i);

      // Jacobian is #outputs x #params
      Matrix<double, 2, PointParamsN> J;

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( unsigned i=0; i < PointParamsN; ++i ){
        Vector<double, PointParamsN> b_i_prime = b_i;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(b_i(i)*1e-7);
        b_i_prime(i) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        Vector2 hi = impl()(i,j,a_j,b_i_prime);
        select_col(J,i) = (hi-h0)/epsilon;
      }
      return J;
    }
  };

  class CameraBundleAdjustmentModel : public BundleAdjustmentModelBase<CameraBundleAdjustmentModel, 6, 3> {
    std::vector<boost::shared_ptr<AdjustedCameraModel> > m_cameras;

  public:
    typedef AdjustedCameraModel camera_type;
    
    CameraBundleAdjustmentModel(std::vector<boost::shared_ptr<CameraModel> > const& cameras) : 
      m_cameras(cameras.size()) {

      for (unsigned i = 0; i < cameras.size(); ++i)
        m_cameras[i] = boost::shared_ptr<AdjustedCameraModel>(new AdjustedCameraModel(cameras[i]));
    }

    void update(std::vector<Vector<double, camera_params_n> > a) {
      for (unsigned j=0; j < a.size(); ++j) {
        Vector3 q1 = subvector(a[j], 0, 3);
        m_cameras[j]->set_axis_angle_rotation(q1);
        m_cameras[j]->set_translation(subvector(a[j], 3,3));
      }
    }
   
    std::vector<boost::shared_ptr<CameraModel> > adjusted_cameras() const {
      std::vector<boost::shared_ptr<CameraModel> > result(m_cameras.size());
      for (unsigned i = 0; i < m_cameras.size(); ++i) {
        result[i] = m_cameras[i];
      }
      return result;
    }

    void write_adjustments(std::vector<std::string> const& filenames) {
      for (unsigned int i=0; i < filenames.size(); ++i) {
        m_cameras[i]->write(filenames[i]+".adjust");
      }
    }

    // Given the 'a' vector (camera model parameters) for the j'th
    // image, and the 'b' vector (3D point location) for the i'th
    // point, return the location of b_i on imager j in pixel
    // coordinates.
    Vector2 operator() ( unsigned i, unsigned j, Vector<double,6> const& a_j, Vector<double,3> const& b_i ) {

      Vector3 q1 = subvector(a_j, 0, 3);
      m_cameras[j]->set_axis_angle_rotation(q1);
      m_cameras[j]->set_translation(subvector(a_j, 3,3));
      
      return m_cameras[j]->point_to_pixel(b_i);
    }    
    
    Vector<double,6> initial_parameters(unsigned j) const { return Vector<double,6>(); }
  };

  //------------------------------------------------------------------
  //                 Sparse Skyline Matrix 
  //
  // The hessian in the sparse LM algorithm takes the form of a Sparse
  // skyline matrix.  Under this sparseness condition, efficient
  // solutions to a linear system can be efficiently computed.  The
  // code below defines a class for a sparse skyline matrix as well as
  // methods for it's decomposition and for linear system solving.
  // -----------------------------------------------------------------

  /// An extremely simple sparse matrix class that wraps around a
  /// boost compessed_matrix<>, but keeps track of the first non-zero
  /// element in each row (i.e. the "skyline") for more efficient
  /// processing in some algorithms.
  template <class ElemT>
  class SparseSkylineMatrix {  
    typedef boost::numeric::ublas::compressed_matrix<ElemT> sparse_matrix_type;
    // Haven't decided yet which boost sparse matrix type is
    // fastest... this one might be good, too.
    //
    // typedef boost::numeric::ublas::generalized_vector_of_vector<ElemT, row_major, vector<coordinate_vector<ElemT> > > sparse_matrix_type;
    
    sparse_matrix_type m_matrix;
    std::vector<uint32> m_skyline;

  public:
    SparseSkylineMatrix(unsigned cols, unsigned rows) : 
      m_matrix(cols,rows), m_skyline(cols) {
      VW_ASSERT(cols == rows, ArgumentErr() << "SparseSkylineMatrix must be square and symmetric.\n");
      for (unsigned i = 0; i < cols; ++i)
        m_skyline[i] = i;
    }
  

    // Returns the "skyline" of the sparse, symmetric matrix.  That is,
    // each index i in the returned vector contains the index of the
    // first valid entry in row i (or equivelently, column i) of the
    // skyline matrix.
    const std::vector<uint32>& skyline() const { return m_skyline; }
  
    uint32 rows() const { return m_matrix.size1(); }
    uint32 cols() const { return m_matrix.size2(); }

    const ElemT* find_element (uint32 i, uint32 j) const { 
      if (j > i) 
        return m_matrix.find_element(j,i); 
      else
        return m_matrix.find_element(i,j); 
    }

    // Some boost sparse matrix types define this method...
    void push_back (uint32 i, uint32 j, ElemT const& val) {
      return m_matrix.push_back(i,j,val);
    }

    typename sparse_matrix_type::const_reference operator () (uint32 i, uint32 j) const {
      VW_DEBUG_ASSERT(i < 0 || i >= this->rows() || j < 0 || j >= this->cols(),
                      ArgumentErr() << "SparseSkylineMatrix: index " << i << " " << j << " out of bounds.");

      // Force symmetry by reflecting all points to the lower left
      // triangle.
      if (j > i) 
        return m_matrix(j,i);
      else
        return m_matrix(i,j);
    }
    
    typename sparse_matrix_type::reference operator () (uint32 i, uint32 j) {
      VW_DEBUG_ASSERT(i < 0 || i >= this->rows() || j < 0 || j >= this->cols(),
                      ArgumentErr() << "SparseSkylineMatrix: index " << i << " " << j << " out of bounds.");
      
      // Force symmetry by reflecting all points to the lower left
      // triangle.
      if (j > i) {
        uint32 temp = j; 
        j = i; 
        i = temp;
      }
      
      if (j < m_skyline[i])
        m_skyline[i] = j;
      return m_matrix(i,j);
    }
    
    // Handy for debugging...
    void print_sparse_structure() {
      std::cout << "SPARSE STRUCTURE: \n";
      for (unsigned i = 0; i < this->rows(); ++i) {
        for (unsigned j = 0; j < this->cols(); ++j) {
          const ElemT* e = this->find_element(i,j);
          if (e) 
            std::cout << "* ";
          else 
            std::cout << ". ";
        }
        std::cout << "\n";
      }
    }    

    // Handy for debugging...
    void print_matrix() {
      std::cout << "SPARSE SKYLINE MATRIX: \n";
      for (unsigned i = 0; i < this->rows(); ++i) {
        for (unsigned j = 0; j < this->cols(); ++j) {
          const ElemT* e = this->find_element(i,j);
          if (e) 
            std::cout << this->operator()(i,j) << "  ";
          else 
            std::cout << "0         ";
        }
        std::cout << "\n";
      }
    }    

  };

  //--------------------------------------------------------------------
  //        L*D*L^T Decompostion for Symmetric, Spares, Skyline Matrices
  //--------------------------------------------------------------------
  
  // Perform L*D*L^T decomposition on a sparse, skyline symmetric
  // semi-definite matrix.  WARNING: The results are stored in place,
  // so this operation destroys the previous contents of A.
  //
  // Once this operation is complete, the diagonal entries of A
  // contain the values from D, and the lower left block diagonal of A
  // contains L.  (The diagonal entries of L are always 1, so those
  // are assumed here...)
  template <class ElemT>
  void ldl_decomposition(SparseSkylineMatrix<ElemT>& A) {
    VW_ASSERT(A.cols() == A.rows(), ArgumentErr() << "ldl_decomposition: argument must be square and symmetric.\n");
    
    const std::vector<uint32>& skyline = A.skyline();
    
    for (unsigned j = 0; j < A.cols(); ++j) {

      // Compute v(1:j)
      std::vector<double> v(j+1);
      v[j] = A(j,j);
      for (unsigned i = skyline[j]; i < j; ++i) {
        v[i] = A(j,i)*A(i,i);
        v[j] -= A(j,i)*v[i];
      }
      
      // Store d(j) and compute L(j+1:n,j)
      A(j,j) = v[j];
      for (unsigned i = j+1; i < A.cols(); ++i) {
        double row_sum = 0;
        for (unsigned jj = skyline[i]; jj < j; ++jj) 
          row_sum += A(i,jj)*v[jj];
        if (j >= skyline[i])
          A(i,j) = ( A(i,j)-row_sum ) / v[j];
      }
    }
  }

  //--------------------------------------------------------------
  //            Solve Spare Skyline Linear System: Ax=b
  //--------------------------------------------------------------

  /// Perform L*D*L^T decomposition on a sparse skyline symmetric
  /// semi-definite matrix A (in place) and then solves an equation of the
  /// form Ax=b using forward and backward substitution.
  /// 
  /// WARNING: Modifies the contents of the matrix A.
  template <class ElemT, class VectorT>
  vw::Vector<double> sparse_solve(SparseSkylineMatrix<ElemT>& A, VectorT const& b) {
    VW_ASSERT(A.cols() == A.rows(), ArgumentErr() << "sparse_solve: matrix must be square and symmetric.\n");

    // Compute the L*D*L^T decomposition of A
    ldl_decomposition(A);
    
    const std::vector<uint32>& skyline = A.skyline();
    std::vector<uint32> inverse_skyline(skyline.size());
    
    // Construct the inverse skyline matrix, which is used to optimize the final
    // back substitution step below.
    for (unsigned j = 0; j < inverse_skyline.size(); ++j) {
      inverse_skyline[j] = 0;
      for (int i = skyline.size()-1; i>=0; --i) {
        if (j < skyline[i]) 
          ++(inverse_skyline[j]);
        else
          break; // Break out of the inner loop
      }
    }
   
    // Forward Substitution Step ( L*x'=b )
    vw::Vector<double> x_prime(A.cols());
    for (unsigned i = 0; i < x_prime.size(); ++i) {
      double sum = 0;
      for (unsigned j = skyline[i]; j < i; ++j) 
        sum += A(i,j)*x_prime(j);
      x_prime(i) = b(i)-sum;
    }
    
    // Divide by D ( D*x''=x' )
    vw::Vector<double> x_doubleprime(A.cols());
    for (unsigned i = 0; i < x_doubleprime.size(); ++i) {
      x_doubleprime(i) = x_prime(i)/A(i,i);
    }

    // Back Substitution step ( L^T*x=x'' )
    vw::Vector<double> x(A.cols());
    for (int32 i = x.size()-1; i >= 0; --i) {
      double sum = 0;
      for (unsigned j = i+1; j < A.cols()-inverse_skyline[i]; ++j) 
        sum += A(j,i)*x(j);
      x(i) = x_doubleprime(i) - sum;
    }
    return x;
  }


  //--------------------------------------------------------------
  //                     BundleAdjustment
  //--------------------------------------------------------------

  // Use this type to describe a "bundle" for a given 3D point.
  // Each element 'i' of the vector contains 1) the index 'j' of a
  // camera where that point was imaged (corresponding to a camera
  // in camera_models), and 2) the pixel location where point 'i'
  // was imaged by imager 'j'.
  struct Bundle {
    typedef std::vector<std::pair<uint32, Vector2> > pixel_list_type;
    Vector3 position;
    pixel_list_type imaged_pixel_locations;
  };
  
  template <class BundleAdjustModelT>
  class BundleAdjustment {   
    
  private:
    ControlNetwork m_control_net;
    //    std::vector<Bundle> m_bundles;
    BundleAdjustModelT &m_model;
    double m_lambda;
    int m_num_observations;
    int m_iterations;
    
    // These are the paramaters to be estimated.  Their values are
    // updated as the optimization progresses.
    std::vector<Vector<double, BundleAdjustModelT::camera_params_n> > a;
    std::vector<Vector<double, BundleAdjustModelT::point_params_n> > b;

  public:
    
    BundleAdjustment(BundleAdjustModelT &model, ControlNetwork const& cnet) : 
      m_control_net(cnet), m_model(model), a(model.size()), b(cnet.size()) {

      m_iterations = 0;

      // Set up the a and b vectors.
      for (unsigned i = 0; i < cnet.size(); ++i) {
        b[i] = cnet[i].position();
      }

      for (unsigned j = 0; j < model.size(); ++j) {
        a[j] = m_model.initial_parameters(j);
      }

      // Set an initial value for lambda.  This value varies as the
      // algorithm converges.
      m_lambda = 1e-25;

      // Compute the number of observations from the bundle.
      m_num_observations = 0;
      for (unsigned i = 0; i < cnet.size(); ++i)
        m_num_observations += cnet[i].size();
    }
//     BundleAdjustment(BundleAdjustModelT &model, std::vector<Bundle> const& bundles) : 
//       m_bundles(bundles), m_model(model), a(model.size()), b(bundles.size()) {

//       m_iterations = 0;

//       // Set up the a and b vectors.
//       for (unsigned i = 0; i < bundles.size(); ++i) {
//         b[i] = m_bundles[i].position;
//       }

//       for (unsigned j = 0; j < model.size(); ++j) {
//         a[j] = m_model.initial_parameters(j);
//       }

//       // Set an initial value for lambda.  This value varies as the
//       // algorithm converges.
//       m_lambda = 1e-25;

//       // Compute the number of observations from the bundle.
//       m_num_observations = 0;
//       for (unsigned i = 0; i < bundles.size(); ++i)
//         m_num_observations += bundles[i].imaged_pixel_locations.size();
//     }
    
    void set_lambda(double lambda) { m_lambda = lambda; }

    // Robust cost functions.  These cost function can help to reduce
    // the impact of outliers in the bundle adjustment.
    double pseudo_huber_error(double delta_norm, double b) {
      return 2.0f * pow(b,2) * (sqrt(1.0f + pow(delta_norm/b,2)) - 1.0f);
    }

    double huber_error(double delta_norm, double b) {
      if (delta_norm < b)
        return delta_norm*delta_norm;
      else
        return 2*b*delta_norm - b*b;
    }

    double squared_error(double delta_norm) {
      return delta_norm*delta_norm;
    }
    
    double L1_error(double delta_norm) {
      return abs(delta_norm);
    }

    // This is a simple, non-sparse, unoptimized implementation of LM
    // bundle adjustment.  It is primarily used for validation and
    // debugging.
    //
    // Each entry in the outer vector corresponds to a distinct 3D
    // point.  The inner vector contains a list of image IDs and
    // pixel coordinates where that point was imaged.
    double update_reference_impl(double &abs_tol, double &rel_tol) { 
      ++m_iterations;
      int num_observations = 2*m_num_observations+6*a.size()+3*b.size();
      int num_parameters = 6*a.size() + 3*b.size();
      
      Matrix<double> J(num_observations, num_parameters);
      Vector<double> epsilon(num_observations);
      
      // --- SETUP STEP ----

      // Add rows to J and epsilon for the imaged pixel observations
      int idx = 0;
      int num_cameras = m_model.size();
      for (unsigned i = 0; i < m_control_net.size(); ++i) {
        Vector<double, BundleAdjustModelT::point_params_n> point_params = b[i];
        for (unsigned j = 0; j < m_control_net[i].size(); ++j) {
          int camera_idx = m_control_net[i][j].image_id();
          Vector<double, BundleAdjustModelT::camera_params_n> camera_params = a[camera_idx];
          
          Matrix<double> J_a = m_model.A_jacobian(i,camera_idx, camera_params, point_params);
          Matrix<double> J_b = m_model.B_jacobian(i,camera_idx, camera_params, point_params);

          // Populate the Jacobian Matrix
          submatrix(J, 2*idx, camera_idx*J_a.cols(), 2, J_a.cols()) = J_a;
          submatrix(J, 2*idx, num_cameras*J_a.cols() + i*J_b.cols(), 2, J_b.cols()) = J_b;

          // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error = m_control_net[i][j].position() - m_model(i, camera_idx, camera_params, point_params);
          double weight = sqrt(huber_error(norm_2(unweighted_error),1.0)) / norm_2(unweighted_error);
          subvector(epsilon,2*idx,2) = unweighted_error * weight;

          ++idx;
        }
      }
      
      // Add rows to J for a priori position/pose constraints
      Matrix<double> id(6*a.size(),6*a.size());
      id.set_identity();
      for (unsigned i=0; i < a.size(); ++i) {
        for (unsigned j=0; j < 3; ++j) {
          unsigned pos = 6*i+j;
          id(pos,pos) *= 0.1;   // position
        }
        for (unsigned j=3; j < 6; ++j) {
          unsigned pos = 6*i+j;
          id(pos,pos) *= 0.9;   // pose
        }
      }
      submatrix(J,2*m_num_observations,0,6*a.size(),6*a.size()) = id;

      // And add error measurement to epsilon
      for (unsigned i=0; i < a.size(); ++i) {
        subvector(epsilon,2*m_num_observations + BundleAdjustModelT::camera_params_n*i,6) = -a[i];
      }


      // --- UPDATE STEP ----
      
      // Build up the right side of the normal equation...
      Vector<double> del_J = -1.0 * (transpose(J) * epsilon);

      // ... and the left side.  (Remembering to rescale the diagonal
      // entries of the approximated hessian by lambda)
      Matrix<double> hessian = (transpose(J) * J);
      for ( unsigned i=0; i < hessian.rows(); ++i )
          hessian(i,i) += hessian(i,i)*m_lambda + m_lambda;

      // Solve for update
      Vector<double> delta = least_squares(hessian, del_J);


      // --- EVALUATE UPDATE STEP ---

      // Recompute error vector
      idx = 0;
      Vector<double> new_epsilon(2*m_num_observations + 6*a.size() + 6*b.size());
      for (unsigned i = 0; i < m_control_net.size(); ++i) {
        Vector<double, BundleAdjustModelT::point_params_n> point_params = b[i];
        for (unsigned j = 0; j < m_control_net[i].size(); ++j) {
          int camera_idx = m_control_net[i][j].image_id();
          Vector<double, BundleAdjustModelT::camera_params_n> camera_params = a[camera_idx];
          
          Vector<double> cam_delta = subvector(delta, BundleAdjustModelT::camera_params_n*camera_idx, BundleAdjustModelT::camera_params_n);
          Vector<double> point_delta = subvector(delta, BundleAdjustModelT::camera_params_n*num_cameras + BundleAdjustModelT::point_params_n*i, BundleAdjustModelT::point_params_n);

          // Apply robust cost function weighting
          Vector2 unweighted_error = m_control_net[i][j].position() - m_model(i, camera_idx, camera_params - cam_delta, point_params - point_delta);
          double weight = sqrt(huber_error(norm_2(unweighted_error),1.0)) / norm_2(unweighted_error);

          // Populate the new error vector
          subvector(new_epsilon,2*idx,2) = unweighted_error * weight;

          ++idx;
        }
      }

      // Populate new error from camera params
      for (unsigned i=0; i < a.size(); ++i) {
        Vector<double> cam_delta = subvector(delta, BundleAdjustModelT::camera_params_n*i, BundleAdjustModelT::camera_params_n);
        subvector(new_epsilon,2*m_num_observations + BundleAdjustModelT::camera_params_n*i, 6) = -(a[i] - cam_delta);
      }

      if (norm_2(new_epsilon) < norm_2(epsilon))  {
        for (unsigned j=0; j<a.size(); ++j) 
          a[j] -= subvector(delta, BundleAdjustModelT::camera_params_n*j, BundleAdjustModelT::camera_params_n);
        for (unsigned i=0; i<b.size(); ++i)
          b[i] -= subvector(delta, BundleAdjustModelT::camera_params_n*num_cameras + BundleAdjustModelT::point_params_n*i, BundleAdjustModelT::point_params_n);
        m_lambda /= 10;

        // Update the camera model params
        m_model.update(a);

        // Summarize the stats from this step in the iteration
        double improvement = (norm_2(epsilon) - norm_2(new_epsilon))/m_num_observations;
        double avg_err = norm_2(new_epsilon)/m_num_observations;
        std::cout << "Reference LM Iteration " << m_iterations << " -- Weighted Error: " << avg_err << "  Delta: " << improvement << "  Lambda: " << m_lambda << "             \n" << std::flush;

        abs_tol = avg_err;
        rel_tol = improvement;
        
        return improvement;
      } else {
        m_lambda *= 10;
        return ScalarTypeLimits<double>::highest();
      }
    }

    // Each entry in the outer vector corresponds to a distinct 3D
    // point.  The inner vector contains a list of image IDs and
    // pixel coordinates where that point was imaged.
    double update_reference_impl2(double &abs_tol, double &rel_tol) { 
      ++m_iterations;
      int num_observations = 2*m_num_observations+6*a.size()+3*b.size();
      int num_parameters = 6*a.size() + 3*b.size();
      
      Matrix<double> J(num_observations, num_parameters);
      Vector<double> epsilon(num_observations);
      
      // --- SETUP STEP ----

      // Add rows to J and epsilon for the imaged pixel observations
      int idx = 0;
      int num_cameras = m_model.size();
      for (unsigned i = 0; i < m_control_net.size(); ++i) {
        Vector<double, BundleAdjustModelT::point_params_n> point_params = b[i];
        for (unsigned j = 0; j < m_control_net[i].size(); ++j) {
          int camera_idx = m_control_net[i][j].image_id();
          Vector<double, BundleAdjustModelT::camera_params_n> camera_params = a[camera_idx];
          
          Matrix<double> J_a = m_model.A_jacobian(i,camera_idx, camera_params, point_params);
          Matrix<double> J_b = m_model.B_jacobian(i,camera_idx, camera_params, point_params);

          // Populate the Jacobian Matrix
          submatrix(J, 2*idx, camera_idx*J_a.cols(), 2, J_a.cols()) = J_a;
          submatrix(J, 2*idx, num_cameras*J_a.cols() + i*J_b.cols(), 2, J_b.cols()) = J_b;

          // Apply robust cost function weighting and populate the error vector
          Vector2 unweighted_error = m_control_net[i][j].position() - m_model(i, camera_idx, camera_params, point_params);
          double weight = sqrt(huber_error(norm_2(unweighted_error),1.0)) / norm_2(unweighted_error);
          subvector(epsilon,2*idx,2) = unweighted_error * weight;

          ++idx;
        }
      }
      
      // Add rows to J for a priori position/pose constraints
      Matrix<double> id(6*a.size(),6*a.size());
      id.set_identity();
      submatrix(J,2*m_num_observations,0,6*a.size(),6*a.size()) = id;

      // And add error measurement to epsilon
      for (unsigned i=0; i < a.size(); ++i) {
        subvector(epsilon,2*m_num_observations + BundleAdjustModelT::camera_params_n*i,6) = -a[i];
      }

      // --- UPDATE STEP ----
      Matrix<double> A = submatrix(J, 0, 0,          2*m_num_observations, 6*a.size());
      Matrix<double> B = submatrix(J, 0, 6*a.size(), 2*m_num_observations, 3*b.size());
      Matrix<double> C = submatrix(J,2*m_num_observations, 0, 6*a.size(), 6*a.size());
      Matrix<double> D = submatrix(J,2*m_num_observations+6*a.size(), 6*a.size(), 3*b.size(), 3*b.size());
      Vector<double> epsilon_o = subvector(epsilon, 0, 2*m_num_observations);
      Vector<double> epsilon_a = subvector(epsilon, 2*m_num_observations, 6*a.size());
      Vector<double> epsilon_b = subvector(epsilon, 2*m_num_observations+6*a.size(), 3*b.size());

      // Add rows to J for a priori position/pose constraints
      Matrix<double> inverse_cov(6*a.size(),6*a.size());
      inverse_cov.set_identity();
      for (unsigned i=0; i < a.size(); ++i) {
        for (unsigned j=0; j < 3; ++j) {
          unsigned pos = 6*i+j;
          inverse_cov(pos,pos) = 1/0.50;   // position
        }
        for (unsigned j=3; j < 6; ++j) {
          unsigned pos = 6*i+j;
          inverse_cov(pos,pos) = 1/1.0;   // pose
        }
      }

      Matrix<double> U = transpose(A)*A + transpose(C) * inverse_cov * C;
      Matrix<double> V = transpose(B)*B;// + transpose(D) * D;
      Matrix<double> W = transpose(A)*B;
      Vector<double> epa = transpose(A)*epsilon_o + transpose(C) * inverse_cov * epsilon_a;
      Vector<double> epb = transpose(B)*epsilon_o;// + transpose(D)*epsilon_b;
      // Augment the diagonals
      for ( unsigned i=0; i < U.rows(); ++i )
        U(i,i) += U(i,i)*m_lambda + m_lambda;
      for ( unsigned i=0; i < V.rows(); ++i )
        V(i,i) += V(i,i)*m_lambda + m_lambda;

      Matrix<double> Y = W*inverse(V);

      // Solve for update
      Vector<double> delta_a = inverse(U-Y*transpose(W)) * (epa - Y*epb);
      Vector<double> delta_b = inverse(V)*(epb-transpose(W)*delta_a);
      Vector<double> delta(delta_a.size() + delta_b.size());
      subvector(delta, 0, delta_a.size()) = delta_a;
      subvector(delta, delta_a.size(), delta_b.size()) = delta_b;


      // --- EVALUATE UPDATE STEP ---

      // Recompute error vector
      idx = 0;
      Vector<double> new_epsilon(2*m_num_observations + 6*a.size() + 6*b.size());
      for (unsigned i = 0; i < m_control_net.size(); ++i) {
        Vector<double, BundleAdjustModelT::point_params_n> point_params = b[i];
        for (unsigned j = 0; j < m_control_net[i].size(); ++j) {
          int camera_idx = m_control_net[i][j].image_id();
          Vector<double, BundleAdjustModelT::camera_params_n> camera_params = a[camera_idx];
          
          Vector<double> cam_delta = subvector(delta, BundleAdjustModelT::camera_params_n*camera_idx, BundleAdjustModelT::camera_params_n);
          Vector<double> point_delta = subvector(delta, BundleAdjustModelT::camera_params_n*num_cameras + BundleAdjustModelT::point_params_n*i, BundleAdjustModelT::point_params_n);

          // Apply robust cost function weighting
          Vector2 unweighted_error = m_control_net[i][j].position() - m_model(i, camera_idx, camera_params + cam_delta, point_params + point_delta);
          double weight = sqrt(huber_error(norm_2(unweighted_error),1.0)) / norm_2(unweighted_error);

          // Populate the new error vector
          subvector(new_epsilon,2*idx,2) = unweighted_error * weight;

          ++idx;
        }
      }

      // Populate new error from camera params
      for (unsigned i=0; i < a.size(); ++i) {
        Vector<double> cam_delta = subvector(delta, BundleAdjustModelT::camera_params_n*i, BundleAdjustModelT::camera_params_n);
        subvector(new_epsilon,2*m_num_observations + BundleAdjustModelT::camera_params_n*i, 6) = -(a[i] + cam_delta);
      }

      if (norm_2(new_epsilon) < norm_2(epsilon))  {
        for (unsigned j=0; j<a.size(); ++j) 
          a[j] += subvector(delta, BundleAdjustModelT::camera_params_n*j, BundleAdjustModelT::camera_params_n);
        for (unsigned i=0; i<b.size(); ++i)
          b[i] += subvector(delta, BundleAdjustModelT::camera_params_n*num_cameras + BundleAdjustModelT::point_params_n*i, BundleAdjustModelT::point_params_n);
        m_lambda /= 10;

        // Update the camera model params
        m_model.update(a);

        // Summarize the stats from this step in the iteration
        double improvement = (norm_2(epsilon) - norm_2(new_epsilon))/m_num_observations;
        double avg_err = norm_2(new_epsilon)/m_num_observations;
        std::cout << "Reference LM Iteration " << m_iterations << " -- Weighted Error: " << avg_err << "  Delta: " << improvement << "  Lambda: " << m_lambda << "             \n" << std::flush;

        abs_tol = avg_err;
        rel_tol = improvement;
        
        return improvement;
      } else {
        m_lambda *= 10;
        return ScalarTypeLimits<double>::highest();
      }
    }

    //------------------------

    // This is the sparse levenberg marquardt update step.  Returns
    // the average improvement in the cost function.
    double update(double &abs_tol, double &rel_tol) { 
      ++m_iterations;

      VW_DEBUG_ASSERT(m_control_net.size() == b.size(), LogicErr() << "BundleAdjustment::update() : Number of bundles does not match the size of the b vector.");

      // Jacobian Matrices and error values
      boost::numeric::ublas::mapped_matrix<Matrix<double, 2, BundleAdjustModelT::camera_params_n> > A(b.size(), a.size());
      boost::numeric::ublas::mapped_matrix<Matrix<double, 2, BundleAdjustModelT::point_params_n> > B(b.size(), a.size());
      boost::numeric::ublas::mapped_matrix<Vector2> epsilon(b.size(), a.size());

      // Intermediate Matrices and vectors
      boost::numeric::ublas::mapped_vector< Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::camera_params_n> > U(a.size()); 
      boost::numeric::ublas::mapped_vector< Matrix<double,BundleAdjustModelT::point_params_n,BundleAdjustModelT::point_params_n> > V(b.size());  
      boost::numeric::ublas::mapped_matrix<Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::point_params_n> > W(b.size(), a.size());

      boost::numeric::ublas::mapped_vector< Vector<double,BundleAdjustModelT::camera_params_n> > epsilon_a(a.size());
      boost::numeric::ublas::mapped_vector< Vector<double,BundleAdjustModelT::point_params_n > > epsilon_b(b.size());
      boost::numeric::ublas::mapped_matrix<Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::point_params_n> > Y(b.size(), a.size());
    
      // Populate the Jacobian, which is broken into two sparse
      // matrices A & B, as well as the error matrix and the W 
      // matrix.
      unsigned i = 0;
      unsigned n = 0;
      double error_total = 0;
      for (typename ControlNetwork::const_iterator iter = m_control_net.begin(); iter != m_control_net.end(); ++iter) {
        for (typename ControlPoint::const_iterator point_iter = (*iter).begin(); point_iter != (*iter).end(); ++point_iter) {
          unsigned j = (*point_iter).image_id();
          VW_DEBUG_ASSERT(j >=0 && j < a.size(), ArgumentErr() << "BundleAdjustment::update() : image index out of bounds.");

          // Store jacobian values
          A(i,j) = m_model.A_jacobian(i,j,a[j],b[i]);
          B(i,j) = m_model.B_jacobian(i,j,a[j],b[i]);

          // Compute error vector
          Vector2 unweighted_error = (*point_iter).position() - m_model(i,j,a[j],b[i]);

          // Apply robust cost function weighting
          double weight = sqrt(huber_error(norm_2(unweighted_error),1.0)) / norm_2(unweighted_error);
          epsilon(i,j).ref() = unweighted_error * weight;
          error_total += pow(epsilon(i,j).ref()(0),2);
          error_total += pow(epsilon(i,j).ref()(1),2);

          // Store intermediate values
          U(j) += transpose(A(i,j).ref()) * /* inverse(SIGMA)* */ A(i,j).ref();
          V(i) += transpose(B(i,j).ref()) * /* inverse(SIGMA)* */ B(i,j).ref();
          W(i,j) = transpose(A(i,j).ref()) * /* inverse(SIGMA)* */ B(i,j).ref();
          epsilon_a(j) += transpose(A(i,j).ref()) * /* inverse(SIGMA)* */ epsilon(i,j).ref();
          epsilon_b(i) += transpose(B(i,j).ref()) * /* inverse(SIGMA)* */ epsilon(i,j).ref();
          ++n;
        }
        ++i;
      }
      error_total = sqrt(error_total);

      
      // Add in the camera position and pose constraint terms.
      for (unsigned j = 0; j < U.size(); ++j) {        
        Matrix<double,BundleAdjustModelT::camera_params_n,BundleAdjustModelT::camera_params_n> inverse_cov;
        for (unsigned jj=0; jj < 3; ++jj) {
          inverse_cov(jj,jj) = 1/0.4;   // position is know to a fraction of a meter
        }
        for (unsigned jj=3; jj < 6; ++jj) {
          inverse_cov(jj,jj) = 1/0.4;   // pose is known to a fraction of a degree
        }

        Matrix<double, BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n> C;
        C.set_identity();
        U(j) += transpose(C) * inverse_cov * C;

        Vector<double, BundleAdjustModelT::camera_params_n> eps_a = -a[j];
        epsilon_a(j) += transpose(C) * inverse_cov * eps_a;
      }


      // Add in the 3D point position constraint terms
      for (unsigned i = 0; i < V.size(); ++i) {        
        Matrix<double, BundleAdjustModelT::point_params_n, BundleAdjustModelT::point_params_n> D;
        V(i) += transpose(D) * D;
      
        Vector<double, BundleAdjustModelT::point_params_n> eps_b = -b[i];
        epsilon_b(i) += transpose(D) * eps_b;
      }

      //      std::cout << "Average error on the image plane for " << n << " points: " << (error_total/n) << " pixels   (total error: " << error_total << ")\n";
 
      // "Augment" the diagonal entries of the U and V matrices with
      // the parameter lambda.
      for (i = 0; i < U.size(); ++i)  
        for (unsigned j = 0; j < BundleAdjustModelT::camera_params_n; ++j)
          U(i).ref()(j,j) += U(i).ref()(j,j)*m_lambda + m_lambda;

      for (i = 0; i < V.size(); ++i) 
        for (unsigned j = 0; j < BundleAdjustModelT::point_params_n; ++j) 
          V(i).ref()(j,j) += V(i).ref()(j,j)*m_lambda + m_lambda;

      // Create the 'e' vector in S * delta_a = e.  The first step is
      // to "flatten" our block structure to a vector that contains
      // scalar entries.
      Vector<double> e(a.size() * BundleAdjustModelT::camera_params_n);
      for (unsigned j = 0; j < epsilon_a.size(); ++j)
        for (unsigned aa = 0; aa < BundleAdjustModelT::camera_params_n; ++aa) {
          e(j*BundleAdjustModelT::camera_params_n + aa) = epsilon_a(j).ref()(aa);
        }

      // Second Pass.  Compute Y and finish constructing e.
      i = 0;
      for (typename ControlNetwork::const_iterator iter = m_control_net.begin(); iter != m_control_net.end(); ++iter) {
        for (typename ControlPoint::const_iterator point_iter = (*iter).begin(); point_iter != (*iter).end(); ++point_iter) {
          unsigned j = (*point_iter).image_id();

          // Compute the blocks of Y
          Y(i,j) = W(i,j).ref() * inverse( V(i).ref() ); 
          
          // "Flatten the block structure to compute 'e'.
          Vector<double, BundleAdjustModelT::camera_params_n> temp = Y(i,j).ref()*epsilon_b(i).ref();
          for (unsigned aa = 0; aa < temp.size(); ++aa)
            e(j*BundleAdjustModelT::camera_params_n + aa) -= temp(aa);
        }
        ++i;
      }
      
      // The S matrix is a m x m block matrix with blocks that are
      // camera_params_n x camera_params_n in size.  It has a sparse
      // skyline structure, which makes it more efficient to solve
      // through L*D*L^T decomposition and forward/back substitution
      // below.
      SparseSkylineMatrix<double> S(a.size()*BundleAdjustModelT::camera_params_n, 
                                    a.size()*BundleAdjustModelT::camera_params_n);

      i = 0;
      for (typename ControlNetwork::const_iterator iter = m_control_net.begin(); iter != m_control_net.end(); ++iter) {
        for (typename ControlPoint::const_iterator j_point_iter = (*iter).begin(); j_point_iter != (*iter).end(); ++j_point_iter) {
          unsigned j = (*j_point_iter).image_id();
          
          //          for (typename ControlPoint::const_iterator k_point_iter = j_point_iter; k_point_iter != (*iter).end(); ++k_point_iter) {
          for (typename ControlPoint::const_iterator k_point_iter = (*iter).begin(); k_point_iter != (*iter).end(); ++k_point_iter) {
            unsigned k = (*k_point_iter).image_id();
            
            // Compute the block entry...
            Matrix<double, BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n> temp = -Y(i,j).ref() * transpose( W(i,k).ref() );
            
            // ... and "flatten" this matrix into the scalar entries of S
            for (unsigned aa = 0; aa < BundleAdjustModelT::camera_params_n; ++aa) {
              for (unsigned bb = 0; bb < BundleAdjustModelT::camera_params_n; ++bb) {
                // FIXME: This if clause is required at the moment to
                // ensure that we do not use the += on the symmetric
                // entries of the SparseSkylineMatrix.  These
                // symmetric entries are shallow, hence this code
                // would add the value twice if we're not careful
                // here.
                if (k*BundleAdjustModelT::camera_params_n + bb <=
                    j*BundleAdjustModelT::camera_params_n + aa) {
                  S(j*BundleAdjustModelT::camera_params_n + aa,
                    k*BundleAdjustModelT::camera_params_n + bb) += temp(aa,bb);
                }
              }
            }
          }
        }
        ++i;
      }
      
      // Augment the diagonal entries S(i,i) with U(i)
      for (unsigned i = 0; i < a.size(); ++i) {
        // ... and "flatten" this matrix into the scalar entries of S
        for (unsigned aa = 0; aa < BundleAdjustModelT::camera_params_n; ++aa) {
          for (unsigned bb = 0; bb < BundleAdjustModelT::camera_params_n; ++bb) {
            // FIXME: This if clause is required at the moment to
            // ensure that we do not use the += on the symmetric
            // entries of the SparseSkylineMatrix.  These
            // symmetric entries are shallow, hence this code
            // would add the value twice if we're not careful
            // here.
            if (i*BundleAdjustModelT::camera_params_n + bb <=
                i*BundleAdjustModelT::camera_params_n + aa) {
              S(i*BundleAdjustModelT::camera_params_n + aa,
                i*BundleAdjustModelT::camera_params_n + bb) += U(i).ref()(aa,bb);
            }
          }
        }
      } 

      // Compute the LDL^T decomposition and solve using sparse methods.
      Vector<double> delta_a = sparse_solve(S, e);
      boost::numeric::ublas::mapped_vector<Vector<double, BundleAdjustModelT::point_params_n> > delta_b(b.size());

      i = 0;
      for (typename ControlNetwork::const_iterator iter = m_control_net.begin(); iter != m_control_net.end(); ++iter) {
        Vector<double, BundleAdjustModelT::point_params_n> temp;
        for (typename ControlPoint::const_iterator j_point_iter = (*iter).begin(); j_point_iter != (*iter).end(); ++j_point_iter) {
          unsigned j = (*j_point_iter).image_id();
          Vector<double, BundleAdjustModelT::camera_params_n> delta_a_j = subvector(delta_a, j*BundleAdjustModelT::camera_params_n, BundleAdjustModelT::camera_params_n);
          temp += transpose( W(i,j).ref() ) * delta_a_j;
        }
        
        delta_b(i) = inverse( V(i).ref() ) * (epsilon_b(i).ref() - temp);
        //        std::cout << "Delta B: " << delta_b(i).ref() << "\n";
        ++i;
      }

      // Compute the update error vector
      i = 0;
      double new_error_total = 0;
      for (typename ControlNetwork::const_iterator iter = m_control_net.begin(); iter != m_control_net.end(); ++iter) {
        for (typename ControlPoint::const_iterator point_iter = (*iter).begin(); point_iter != (*iter).end(); ++point_iter) {
          unsigned j = (*point_iter).image_id();
          
          // Compute error vector
          Vector<double, BundleAdjustModelT::camera_params_n> new_a = a[j] + subvector(delta_a, BundleAdjustModelT::camera_params_n*j, BundleAdjustModelT::camera_params_n);
          Vector<double, BundleAdjustModelT::point_params_n> new_b = b[i] + delta_b(i).ref();
          Vector2 unweighted_error = (*point_iter).position() - m_model(i,j,new_a,new_b);

          // Apply robust cost function weighting
          double weight = sqrt(huber_error(norm_2(unweighted_error),1.0)) / norm_2(unweighted_error);
          Vector2 weighted_error = weight * unweighted_error;
          new_error_total += pow(weighted_error(0),2);
          new_error_total += pow(weighted_error(1),2);
        }
        ++i;
      }
      new_error_total = sqrt(new_error_total);
      
      if (new_error_total < error_total)  {
        for (unsigned j=0; j<a.size(); ++j) 
          a[j] += subvector(delta_a, BundleAdjustModelT::camera_params_n*j, BundleAdjustModelT::camera_params_n);
        for (unsigned i=0; i<b.size(); ++i)
          b[i] += delta_b(i).ref();
        m_lambda /= 10;

        // Update the camera model params
        m_model.update(a);

        // Summarize the stats from this step in the iteration
        double improvement = (error_total - new_error_total)/m_num_observations;
        double avg_err = new_error_total/m_num_observations;
        std::cout << "Sparse LM Iteration " << m_iterations << " -- Error: " << avg_err << "  Delta: " << improvement << "  Lambda: " << m_lambda << "             \n" << std::flush;

        abs_tol = avg_err;
        rel_tol = improvement;
        
        return improvement;
      } else {
        m_lambda *= 10;
        return ScalarTypeLimits<double>::highest();
      }

    }

    BundleAdjustModelT& bundle_adjust_model() { return m_model; }
  };
  
}} // namespace vw::camera

#endif // __VW_CAMERA_BUNDLE_ADJUST_H__
