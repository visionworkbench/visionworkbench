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
#include <vw/Core/Debugging.h>

// Boost 
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace vw {
namespace math {

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
  
}} // namespace vw::math

#endif // __VW_CAMERA_BUNDLE_ADJUST_H__
