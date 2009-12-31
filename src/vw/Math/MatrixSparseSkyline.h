// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file SparseSkylineMatrix.h
///
/// Provides sparse matrix support. This is made in mind for support
/// of Bundle Adjustment
///

#ifndef __VW_MATH_SPARSE_SKYLINE_MATRIX_H__
#define __VW_MATH_SPARSE_SKYLINE_MATRIX_H__

// Standard
#include <vector>

// Vision Workbench
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

// Boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>

namespace vw {
namespace math {

  //------------------------------------------------------------------
  //                 Sparse Skyline Matrix
  //
  // The hessian in the sparse LM algorithm takes the form of a Sparse
  // skyline matrix.  Under this sparseness condition, efficient
  // solutions to a linear system can be efficiently computed.  The
  // code below defines a class for a sparse skyline matrix as well as
  // methods for it's decomposition and for linear system solving.
  // -----------------------------------------------------------------

  // An sparse matrix class that wraps around a boost generalized
  // vector of compressed vectors, but keeps track of the first non-zero
  // element in each row (i.e. the "skyline") for more efficient processing
  // in some algorithms.

  // Skyline Structure
  // #..#..              #        0
  // .#.#.#              .#       1
  // ..##.#              ..#      2
  // ###### is stored as ####   & 0
  // ...###              ...##    3
  // .#####              .#####   1
  // "fingers toward zero"

  template <class ElemT>
  class MatrixSparseSkyline : public MatrixBase<MatrixSparseSkyline<ElemT> > {

    typedef boost::numeric::ublas::compressed_vector<ElemT> compressed_vec;
    typedef boost::numeric::ublas::vector<compressed_vec> inner_type;
    typedef boost::numeric::ublas::row_major order_type;
    typedef boost::numeric::ublas::generalized_vector_of_vector<ElemT, order_type, inner_type > sparse_matrix_type;

    sparse_matrix_type m_matrix;
    Vector<unsigned> m_skyline;

  public:
    /// Common Definitions
    typedef ElemT value_type;
    typedef typename sparse_matrix_type::reference reference_type;
    typedef typename sparse_matrix_type::const_reference const_reference_type;
    typedef ElemT* iterator;
    typedef ElemT const* const_iterator;
    typedef typename sparse_matrix_type::iterator1 sparse_iterator1;
    typedef typename sparse_matrix_type::iterator2 sparse_iterator2;
    typedef typename sparse_matrix_type::const_iterator1 const_sparse_iterator1;
    typedef typename sparse_matrix_type::const_iterator2 const_sparse_iterator2;

    /// Constructors
    MatrixSparseSkyline() :
    m_matrix(0,0), m_skyline(0) {
    }

    MatrixSparseSkyline( unsigned cols, unsigned rows ) :
    m_matrix(cols,rows), m_skyline(cols) {
      VW_ASSERT( cols == rows,
                 ArgumentErr() << "MatrixSparseSkyline must be square and symmetric.\n");
      // Setting non-zero to the identity point;
      for (unsigned i = 0; i < cols; ++i)
        m_skyline[i] = i;
    }

    /// Assignment operator specialized for sparse matrices
    template <class ElemT2>
    MatrixSparseSkyline& operator=( MatrixSparseSkyline<ElemT2> const& m ) {
      VW_ASSERT( m.rows() == m_matrix.size1() &&
                 m.cols() == m_matrix.size2(), ArgumentErr() << "Matrix must have dimensions "
                 << m_matrix.size1() << "x" << m_matrix.size2() << "." );
      // Iterate through non-zero elements
      for ( typename MatrixSparseSkyline<ElemT>::const_sparse_iterator1 it1 = m.sparse_begin();
            it1 != m.sparse_end(); it1++ ) {
        for ( typename MatrixSparseSkyline<ElemT>::const_sparse_iterator2 it2 = it1.begin();
              it2 != it1.end(); it2++ ) {
          (*this)(it2.index1(),it2.index2()) = *it2;
        }
      }
    }

    /// Assignment operator for generic (DESTROYS SPARSITY)
    template <class T>
    MatrixSparseSkyline& operator=( MatrixBase<T> const& m ) {
      VW_ASSERT( m.impl().rows() == m_matrix.size1() &&
                 m.impl().cols() == m_matrix.size2(), ArgumentErr() << "Matrix must have dimensions "
                 << m_matrix.size1() << "x" << m_matrix.size2() << "." );
      vw_out(WarningMessage, "math") << "Sparsity destroyed in generic assignment to MatrixSparseSkyline.\n";
      for ( unsigned i = 0; i < m.rows(); i++ ) {
        for ( unsigned j = 0; j < m.cols(); j++ ) {
          (*this)(i,j) = m(i,j);
        }
      }
    }

    /// Simple Access
    unsigned rows() const { return m_matrix.size1(); }
    unsigned cols() const { return m_matrix.size2(); }

    /// Change the size of the matrix
    void set_size( unsigned new_rows, unsigned new_cols, bool preserve = false ) {
      vw_throw( NoImplErr() << "MatrixSparseSkyline::set_size, code has not been written yet." );
    }

    /// Element Access
    reference_type operator()( unsigned i, unsigned j ) {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( i < rows() && j < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      if ( j > i )
        std::swap(i,j);
      if ( j < m_skyline[i] )
        m_skyline(i) = j;
      return m_matrix(i,j);
    }

    /// Element Access
    const_reference_type operator()( unsigned i, unsigned j ) const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( i < rows() && j < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      if ( j > i )
        return m_matrix(j,i);
      return m_matrix(i,j);
    }

    /// Pointer Access
    value_type *data() {
      vw_throw( NoImplErr() << "MatrixSparseSkyline::data can not provide direct pointer access.");
      return &(operator()(0,0));
    }
    iterator begin() {
      vw_throw( NoImplErr() << "MatrixSparseSkyline::begin can not provide direct pointer access.");
      return &(operator()(0,0));
    }
    iterator end() {
      vw_throw( NoImplErr() << "MatrixSparseSkyline::end can not provide direct pointer access.");
      return &(operator()(0,0));
    }
    sparse_iterator1 sparse_begin() { return m_matrix.begin1(); }
    const_sparse_iterator1 sparse_begin() const { return m_matrix.begin1(); }
    sparse_iterator1 sparse_end() { return m_matrix.end1(); }
    const_sparse_iterator1 sparse_end() const { return m_matrix.end1(); }

    // Special Access to Skyline Vector
    const Vector<unsigned>& skyline() const { return m_skyline; }

  };

  /// Dump for MatrixSparseSkyline
  template <class ElemT>
  inline std::ostream& operator<<( std::ostream& os, MatrixSparseSkyline<ElemT> const& m ) {
    // Need to first determine the sampling rate. In bundle adjustment
    // it doesn't make sense to show every element as camera variables will
    // come in blocks of 6.
    uint smallest_sampling_rate = m.cols();
    uint curr_sampling_rate = m.cols();
    unsigned last_value = 10000;
    Vector<unsigned> skyline = m.skyline();
    for ( uint i = 0; i < skyline.size(); i++ ) {
      if ( last_value != skyline(i) ) {
        if ( smallest_sampling_rate  > curr_sampling_rate )
          smallest_sampling_rate = curr_sampling_rate;
        curr_sampling_rate = 0;
      }
      curr_sampling_rate++;
      last_value = skyline(i);
    }

    os << "MatrixSparseSkyline" << m.rows() << "x" << m.cols() << "@" << smallest_sampling_rate << "\n";
    for ( unsigned i = 0, ii=0; i < m.rows(); i += smallest_sampling_rate, ii++ ) {
      for ( unsigned j = 0, jj=0; j < m.cols(); j += smallest_sampling_rate, jj++ ) {
        ElemT e = m(i,j);
        if (e)
          os << "#";
        else
          if ( ii%2 == 0 && jj%2 == 0 )
            os << ".";
          else
            os << " ";
      }
      os << "\n";
    }
    return os;
  }

  /*
  template <class ElemT>
  class SparseSkylineMatrix {
    //typedef boost::numeric::ublas::compressed_matrix<ElemT> sparse_matrix_type;
    //typedef boost::numeric::ublas::generalized_vector_of_vector<ElemT, row_major, vector<coordinate_vector<ElemT> > > sparse_matrix_type;
    typedef boost::numeric::ublas::compressed_vector<ElemT> compressed_vec;
    typedef boost::numeric::ublas::vector<compressed_vec> inner_type;
    typedef boost::numeric::ublas::row_major order_type;
    typedef boost::numeric::ublas::generalized_vector_of_vector<ElemT, order_type, inner_type > sparse_matrix_type;

    sparse_matrix_type m_matrix;
    std::vector<vw::uint32> m_skyline;

  public:
    SparseSkylineMatrix(unsigned cols, unsigned rows) :
    m_matrix(cols,rows), m_skyline(cols) {
      VW_ASSERT(cols == rows,
                ArgumentErr() << "SparseSkylineMatrix must be square and symmetric.\n");
      // Setting non-zero to the identity point;
      for (unsigned i = 0; i < cols; ++i)
        m_skyline[i] = i;
    }


    // Returns the "skyline" of the sparse, symmetric matrix.  That is,
    // each index i in the returned vector contains the index of the
    // first valid entry in row i (or equivelently, column i) of the
    // skyline matrix.
    const std::vector<vw::uint32>& skyline() const { return m_skyline; }

    vw::uint32 rows() const { return m_matrix.size1(); }
    vw::uint32 cols() const { return m_matrix.size2(); }

    // Some boost sparse matrix types define this method...
    void push_back (vw::uint32 i, vw::uint32 j, ElemT const& val) {
      return m_matrix.push_back(i,j,val);
    }

    typename sparse_matrix_type::const_reference operator () (vw::uint32 i, vw::uint32 j) const {
      VW_DEBUG_ASSERT(i < 0 || i >= this->rows() || j < 0 || j >= this->cols(),
                      ArgumentErr() << "SparseSkylineMatrix: index " << i << " "
                      << j << " out of bounds.");

      // Force symmetry by reflecting all points to the lower left
      // triangle.
      if (j > i)
        return m_matrix(j,i);
      else
        return m_matrix(i,j);
    }

    // REWRITE!
    typename sparse_matrix_type::reference operator () (vw::uint32 i, vw::uint32 j) {
      VW_DEBUG_ASSERT(i < 0 || i >= this->rows() || j < 0 || j >= this->cols(),
                      ArgumentErr() << "SparseSkylineMatrix: index " << i << " "
                      << j << " out of bounds.");

      // Force symmetry by reflecting all points to the lower left
      // triangle.
      if (j > i)
        std::swap(i,j);

      if (j < m_skyline[i])
        m_skyline[i] = j;
      return m_matrix(i,j);
    }

    // Handy for debugging...
    void print_sparse_structure(int sample_rate = 1) {
      vw_out(0) << "SPARSE STRUCTURE: \n";
      for (unsigned i = 0, ii=0; i < this->rows(); i += sample_rate, ii++) {
        for (unsigned j = 0, jj=0; j < this->cols(); j += sample_rate, jj++) {
          ElemT e = (*this)(i,j);
          if (e)
            vw_out(0) << "#";
          else
            if ( ii%2 == 0 &&
                 jj%2 == 0 )
              vw_out(0) << ".";
            else
              vw_out(0) << " ";
        }
        vw_out(0) << "\n";
      }
    }

  };
  */

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
  void sparse_ldl_decomposition(MatrixSparseSkyline<ElemT>& A) {
    VW_ASSERT(A.cols() == A.rows(),
              ArgumentErr() << "ldl_decomposition: argument must be square and symmetric.\n");

    Vector<unsigned> const& skyline = A.skyline();

    for (unsigned j = 0; j < A.cols(); ++j) {

      // Compute v(1:j)
      std::vector<double> v(j+1);
      v[j] = A(j,j);
      for (unsigned i = skyline(j); i < j; ++i) {
        v[i] = A(j,i)*A(i,i);
        v[j] -= A(j,i)*v[i];
      }

      // Store d(j) and compute L(j+1:n,j)
      A(j,j) = v[j];
      for (unsigned i = j+1; i < A.cols(); ++i) {
        double row_sum = 0;
        for (unsigned jj = skyline(i); jj < j; ++jj)
          row_sum += A(i,jj)*v[jj];
        if (j >= skyline(i))
          A(i,j) = ( A(i,j)-row_sum ) / v[j];
      }
    }
  }

  // This version excepts an outside skyline matrix
  template <class MatrixT, class VectorT>
  void sparse_ldl_decomposition(MatrixBase<MatrixT>& A, VectorT& skyline ) {
    for (unsigned j = 0; j < A.impl().cols(); ++j) {

      // Compute v(1:j)
      std::vector<double> v(j+1);
      v[j] = A.impl()(j,j);
      for (unsigned i = skyline(j); i < j; ++i) {
        v[i] = A.impl()(j,i)*A.impl()(i,i);
        v[j] -= A.impl()(j,i)*v[i];
      }

      // Store d(j) and compute L(j+1:n,j)
      A.impl()(j,j) = v[j];
      for (unsigned i = j+1; i < A.impl().cols(); ++i) {
        double row_sum = 0;
        for (unsigned jj = skyline(i); jj < j; ++jj)
          row_sum += A.impl()(i,jj)*v[jj];
        if (j >= skyline(i))
          A.impl()(i,j) = ( A.impl()(i,j)-row_sum ) / v[j];
      }
    }
  }

  //--------------------------------------------------------------
  //            Solve Spare Skyline Linear System: Ax=b
  //--------------------------------------------------------------

  /// Perform L*D*L^T decomposition on a sparse skyline symmetric
  /// semi-definite matrix A (in place) and then call
  /// forward/backward solver

  /// WARNING: Modifies the contents of the matrix A.
  template <class ElemT, class VectorT>
  Vector<ElemT> sparse_solve(MatrixSparseSkyline<ElemT>& A, VectorT const& b ) {
    VW_ASSERT(A.cols() == A.rows(), ArgumentErr() << "sparse_solve: matrix must be square and symmetric.\n");

    // Compute the L*D*L^T decomposition of A
    sparse_ldl_decomposition(A);
    Vector<ElemT> x = sparse_solve_ldl(A, b);

    return x;
  }

  /// This version can use an outside skyline matrix
  template <class MatrixT, class VectorT, class VectorSkyT>
  Vector<typename MatrixT::value_type> sparse_solve(MatrixBase<MatrixT>& A, VectorT const& b, VectorSkyT const& sky ) {
    VW_ASSERT(A.impl().cols() == A.impl().rows(), ArgumentErr() << "sparse_solve: matrix must be square and symmetric.\n");

    // Compute the L*D*L^T decomposition of A
    sparse_ldl_decomposition(A, sky);
    Vector<typename MatrixT::value_type> x = sparse_solve_ldl(A, b, sky);

    return x;
  }

 //--------------------------------------------------------------
  //            Solve Spare Skyline Linear System: AX=B, where
  //            X and B are matrices
  //--------------------------------------------------------------

  /// Perform L*D*L^T decomposition on a sparse skyline symmetric
  /// semi-definite matrix A (in place) and then call
  /// forward/backward solver

  /*
  /// WARNING: Modifies the contents of the matrix A.
  template <class ElemT, class BMatrixT>
  Matrix<typename PromoteType<typename BMatrixT::value_type, typename BMatrixT::value_type>::type>
    multi_sparse_solve(SparseSkylineMatrix<ElemT>& A, BMatrixT & B ) {
    VW_ASSERT(A.cols() == A.rows(), ArgumentErr() << "multi_sparse_solve: matrix must be square and symmetric.\n");
    VW_ASSERT(A.rows() == B.rows(), ArgumentErr() << "multi_sparse_solve: AX=B means A, B have same # of rows.\n");

    Matrix<double> X(B.rows(), B.cols());

    Vector<double> current_col(A.cols());


    // Compute the L*D*L^T decomposition of A
    ldl_decomposition(A);

    for(unsigned i = 0; i < B.cols(); i++){
      current_col = select_col(B, i);
      select_col(X, i) = sparse_solve_ldl(A, current_col);
    }
    return X;
  }
  */


  //-------------------------------------------------------------------
  //            Solve Spare Skyline Linear System: Ax=b Given LDL^T
  //-------------------------------------------------------------------

  /// Perform L*D*L^T decomposition on a sparse skyline symmetric
  /// semi-definite matrix A (in place) and then
  /// Solves an equation of the
  /// form Ax=b using forward and backward substitution.
  ///
  /// Assumes it receives LDL^T form of A
  template <class ElemT, class VectorT>
  Vector<ElemT> sparse_solve_ldl(MatrixSparseSkyline<ElemT>& A, VectorT const& b ) {
    Vector<unsigned> skyline = A.skyline();
    return sparse_solve_ldl(A,b,skyline);
  }

  /// Version that excepts an outside skyline vector
  template <class MatrixT, class VectorT, class VectorSkyT>
  Vector<typename MatrixT::value_type> sparse_solve_ldl(MatrixBase<MatrixT>& A, VectorT const& b, VectorSkyT const& sky) {
    VW_ASSERT(A.impl().cols() == A.impl().rows(), ArgumentErr() << "sparse_solve: matrix must be square and symmetric.\n");

    //const std::vector<vw::uint32>& skyline = A.skyline();
    Vector<unsigned> inverse_sky(sky.size());

    // Construct the inverse skyline, which is used to optimize the final
    // back substitution step below.
    for (unsigned j = 0; j < inverse_sky.size(); ++j) {
      inverse_sky(j) = 0;
      for (int i = sky.size()-1; i>=0; --i) {
        if (j < sky(i))
          ++(inverse_sky(j));
        else
          break; // Break out of the inner loop
      }
    }

    // Forward Substitution Step ( L*x'=b )
    Vector<typename MatrixT::value_type> x_prime(A.impl().cols());
    for (unsigned i = 0; i < x_prime.size(); ++i) {
      typename MatrixT::value_type sum = 0;
      for (unsigned j = sky(i); j < i; ++j)
        sum += A.impl()(i,j)*x_prime(j);
      x_prime(i) = b(i)-sum;
    }

    // Divide by D ( D*x''=x' )
    Vector<typename MatrixT::value_type> x_doubleprime(A.impl().cols());
    for (unsigned i = 0; i < x_doubleprime.size(); ++i)
      x_doubleprime(i) = x_prime(i)/A.impl()(i,i);

    // Back Substitution step ( L^T*x=x'' )
    Vector<typename MatrixT::value_type> x(A.impl().cols());
    for (int32 i = x.size()-1; i >= 0; --i) {
      typename MatrixT::value_type sum = 0;
      for (unsigned j = i+1; j < A.impl().cols()-inverse_sky[i]; ++j)
        sum += A.impl()(j,i)*x(j);
      x(i) = x_doubleprime(i) - sum;
    }

    return x;
  }


}} // namespace vw::math

#endif//__VW_MATH_SPARSE_SKYLINE_MATRIX_H__
