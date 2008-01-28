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
  template <class ImplT, int CameraParamsN, int PointParamsN>
  struct BundleAdjustmentModelBase {

    /// \cond INTERNAL
    // Methods to access the derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
    /// \endcond
   

    inline Vector2 error( unsigned i, unsigned j, 
                          Vector2 const& x_ij, 
                          Vector<double, CameraParamsN> const& a_j, 
                          Vector<double, PointParamsN> const& b_i ) {
      return x_ij - impl()(i,j,a_j,b_i);
    }

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
        Vector<double> b_i_prime = b_i;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(b_i(i)*1e-7);
        b_i_prime(i) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        Vector2 hi = impl()(i,j,b_i_prime, b_i);
        select_col(J,i) = (hi-h0)/epsilon;
      }
      return J;
    }
  };

  class CameraBundleAdjustmentModel : BundleAdjustmentModelBase<CameraBundleAdjustmentModel, 7, 3> {
    std::vector<boost::shared_ptr<AdjustedCameraModel> > m_cameras;

  public:
    typedef AdjustedCameraModel camera_type;
    
    CameraBundleAdjustmentModel(std::vector<boost::shared_ptr<CameraModel> > cameras) : m_cameras(cameras.size()) {
      
      // populate...

    }
   
    // Given the 'a' vector (camera model parameters) for the j'th
    // image, and the 'b' vector (3D point location) for the i'th
    // point, return the location of b_i on imager j in pixel
    // coordinates.
    Vector2 operator() ( unsigned i, unsigned j, Vector<double,7> const& a_j, Vector<double,3> const& b_i ) {

      Vector4 q1 = normalize(subvector(a_j, 0, 4));
      m_cameras[j]->set_rotation(Quaternion<double>(q1[0], q1[1], q1[2], q1[3]));
      m_cameras[j]->set_translation(subvector(a_j, 4,3));
      
      return m_cameras[j]->point_to_pixel(b_i);
    }    
    
    Vector<double,7> camera_params_as_vector(unsigned j) const {
      Vector<double,7> result;
      Quaternion<double> rot = m_cameras[j]->rotation();
      Vector3 trans = m_cameras[j]->translation();
      result(0) = rot[0];
      result(1) = rot[1];
      result(2) = rot[2];
      result(3) = rot[3];
      result(4) = trans(0);
      result(5) = trans(1);
      result(6) = trans(2);
      return result;
    }

    std::vector<boost::shared_ptr<AdjustedCameraModel> > camera_models() {
      return m_cameras;
    }
  };

  //--------------------------------------------------------------
  //                 Sparse Skyline Matrix 
  //--------------------------------------------------------------

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

    const ElemT* find_element (uint32 i, uint32 j) const { return m_matrix.find_element(i,j); }

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
  };

  //--------------------------------------------------------------
  //        L*D*L^T Decompostion for Symmetric Matrices
  //--------------------------------------------------------------
  

  // Perform L*D*L^T decomposition on a symmetric semi-definite
  // matrix.  WARNING: The results are stored in place, so this
  // operation destroys the previous contents of A.
  //
  // Once this operation is complete, the diagonal entries of A
  // contain the values from D, and the lower left block diagonal of A
  // contains L.  (The diagonal entries of L are always 1, so those
  // are assumed here...)
  template <class MatrixT>
  void ldl_decomposition(MatrixT& A) {
    VW_ASSERT(A.cols() == A.rows(), ArgumentErr() << "ldl_decomposition: argument must be square and symmetric.\n");
    for (unsigned j = 0; j < A.cols(); ++j) {
      
    // Compute v(1:j)
    std::vector<double> v(j+1);
    v[j] = A(j,j);
    for (unsigned i = 0; i < j; ++i) {
      v[i] = A(j,i)*A(i,i);
      v[j] -= A(j,i)*v[i];
    }
    
    // Store d(j) and compute L(j+1:n,j)
    A(j,j) = v[j];
    for (unsigned i = j+1; i < A.cols(); ++i) {
      double row_sum = 0;
      for (unsigned jj = 0; jj < j; ++jj) 
        row_sum += A(i,jj)*v[jj];
      A(i,j) = ( A(i,j)-row_sum ) / v[j];
    }
  }
}

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

  template <class BundleAdjustModelT>
  class BundleAdjustment {   

  public:

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

    typedef std::vector<std::pair<uint32, Vector2> > bundle_type;

  private:
    std::vector<Bundle> m_bundles;
    BundleAdjustModelT m_model;

    // These are the paramaters to be estimated.  Their values are
    // updated as the optimization progresses.
    std::vector<Vector<double, 7> > a;
    std::vector<Vector<double, 3> > b;
  public:

    BundleAdjustment(std::vector<boost::shared_ptr<CameraModel> > const& camera_models,
                     std::vector<Bundle> const& bundles) : 
      m_bundles(bundles), m_model(camera_models), a(camera_models.size()), b(bundles.size()) {

      // Set up the a and b vectors.
      for (unsigned i = 0; i < bundles.size(); ++i) {
        b[i] = m_bundles[i].position;
      }

      for (unsigned j = 0; j < camera_models.size(); ++j) {
        a[j] = m_model.camera_params_as_vector(j);
      }
    }

    // Returns a list of "adjusted" camera models that resulted from
    // the bundle adjustment.
    std::vector<boost::shared_ptr<typename BundleAdjustModelT::camera_type> > camera_models() {
      return m_model.camera_models();
    }    

//     // Each entry in the outer vector corresponds to a distinct 3D
//     // point.  The inner vector contains a list of image IDs and
//     // pixel coordinates where that point was imaged.
    
//     void update(std::vector<Vector<double> > const& a, std::vector<Vector<double> > const& b, bundles_type const& points, int num_images, double lambda) const {

//       mapped_matrix<Matrix<double> > A(points.size(), num_images);
//       mapped_matrix<Matrix<double> > B(points.size(), num_images);
//       mapped_matrix<Vector2> epsilon(points.size(), num_images);

//       mapped_vector< Matrix<double> > U(A.size1());  // size1() == rows
//       mapped_vector< Matrix<double> > V(B.size1());  
//       mapped_matrix<Matrix<double> > W(points.size(), num_images);
//       mapped_vector< Vector2 > epsilon_a(A.size1());
//       mapped_vector< Vector2 > epsilon_b(B.size1());
//       mapped_matrix<Matrix<double> > Y(points.size(), num_images);
//       //      mapped_matrix<Matrix<double> > S(points.size(), num_image);

    
//       // Populate the Jacobian, which is broken into two sparse
//       // matrices A & B, as well as the error matrix and the W 
//       // matrix.
//       int i = 0;
//       for (bundles_type::const_iterator iter = points.begin(); iter != points.end(); ++iter) {
//         for (bundle_type::const_iterator point_iter = (*iter).begin(); point_iter != (*iter).end(); ++point_iter) {
//           int j = (*point_iter).first;

//           // Store jacobian values
//           A(i,j) = A_jacobian(a[j],b[i]);
//           B(i,j) = B_jacobian(a[j],b[i]);

//           // Compute error vector
//           epsilon(i,j) = this->error((*point_iter).second,a[j],b[i]);

//           // Store intermediate values
//           U(j) += transpose(A(i,j).ref()) * /* inverse(SIGMA)* */ A(i,j).ref();
//           V(i) += transpose(B(i,j).ref()) * /* inverse(SIGMA)* */ B(i,j).ref();
//           W(i,j) = transpose(A(i,j).ref()) * /* inverse(SIGMA)* */ B(i,j).ref();
//           epsilon_a(j) += transpose(A(i,j).ref()) * /* inverse(SIGMA)* */ epsilon(i,j).ref();
//           epsilon_b(i) += transpose(B(i,j).ref()) * /* inverse(SIGMA)* */ epsilon(i,j).ref();

//           /// MUST ONLY AUGMENT DIAGONAL ENTRIES!!!
//           //          Y(i,j) = W(i,j).ref() * inverse( V(i).ref() * (1.0f+lambda) ); 

//           ++i;
//         }
//       }

//       // Second pass
//       //
//       //  NEED TO INCORPORATE THE CORRECT USE OF K!!!
//       //

// //       i = 0;
// //       for (bundles_type::iterator iter = points.begin(); iter != points.end(); ++iter) {
// //         for (bundle_type::iterator point_iter = *iter.begin(); point_iter != *iter.end(); ++point_iter) {
// //           int j = (*point_iter).first;
// //           if (i==j) 
//           /// MUST ONLY AUGMENT DIAGONAL ENTRIES!!!
// //             S(i,j) -= ( Y(i,j).ref()*transpose(W(i,j).ref()) + (U(j).ref() * (1.0f+lambda)) );
// //           else  
// //             S(i,j) -= ( Y(i,j).ref()*transpose(W(i,j).ref()) );
          
// //           ea(j) += epsilon_a(j).ref() - (Y(i,j).ref() * epsilon_b(i).ref());
// //           ++i;
// //         }
// //       }
      
//       //       mapped_vector<double> delta_a = inverse(S) * ea;

//       // Third pass
// //       i = 0;
// //       for (bundles_type::iterator iter = points.begin(); iter != points.end(); ++iter) {
// //         for (bundle_type::iterator point_iter = *iter.begin(); point_iter != *iter.end(); ++point_iter) {
// //           eb(j) += epsilon_b(i).ref() - (transpose(W(i,j).ref()) * delta_a(j).ref());
// //           ++i;
// //         }
// //       }

//       /// MUST ONLY AUGMENT DIAGONAL ENTRIES!!!
//       //      mapped_vector<double> delta_b = inverse( V[i] * (1.0f+lambda) ) * eb;
//     }

  };




  
}} // namespace vw::math

#endif // __VW_CAMERA_BUNDLE_ADJUST_H__
