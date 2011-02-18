// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
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
#include <vw/Core/Log.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

// Boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

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
    Vector<size_t> m_skyline;

  public:
    /// Common Definitions
    typedef ElemT value_type;
    typedef typename sparse_matrix_type::reference reference_type;
    typedef typename sparse_matrix_type::const_reference const_reference_type;
    typedef IndexingMatrixIterator<MatrixSparseSkyline<ElemT> > iterator;
    typedef IndexingMatrixIterator<const MatrixSparseSkyline<ElemT> > const_iterator;
    typedef typename sparse_matrix_type::iterator1 sparse_iterator1;
    typedef typename sparse_matrix_type::iterator2 sparse_iterator2;
    typedef typename sparse_matrix_type::const_iterator1 const_sparse_iterator1;
    typedef typename sparse_matrix_type::const_iterator2 const_sparse_iterator2;

    /// Constructors
    MatrixSparseSkyline() :
    m_matrix(0,0), m_skyline(0) {
    }

    MatrixSparseSkyline( size_t size ) :
    m_matrix(size,size), m_skyline(size) {
      for ( size_t i = 0; i < size; ++i )
        m_skyline[i] = i;
    }

    MatrixSparseSkyline( size_t rows, size_t cols ) :
    m_matrix(rows,cols), m_skyline(rows) {
      VW_ASSERT( cols == rows,
                 ArgumentErr() << "MatrixSparseSkyline must be square and symmetric.\n");
      // Setting non-zero to the identity point;
      for (size_t i = 0; i < rows; ++i)
        m_skyline[i] = i;
    }

    /// Assignment operator specialized for sparse matrices
    template <class ElemT2>
    MatrixSparseSkyline& operator=( MatrixSparseSkyline<ElemT2> const& m ) {
      VW_ASSERT( m.rows() == m_matrix.size1() &&
                 m.cols() == m_matrix.size2(), ArgumentErr() << "Matrix must have dimensions "
                 << m_matrix.size1() << "x" << m_matrix.size2() << "." );
      m_matrix.clear();
      VectorClearImpl<Vector<size_t> >::clear(m_skyline);
      // Iterate through non-zero elements
      for ( typename MatrixSparseSkyline<ElemT>::const_sparse_iterator1 it1 = m.sparse_begin();
            it1 != m.sparse_end(); it1++ ) {
        for ( typename MatrixSparseSkyline<ElemT>::const_sparse_iterator2 it2 = it1.begin();
              it2 != it1.end(); it2++ ) {
          (*this)(it2.index1(),it2.index2()) = *it2;
        }
      }
      return *this;
    }

    /// Assignment operator for generic (DESTROYS SPARSITY)
    template <class T>
    MatrixSparseSkyline& operator=( MatrixBase<T> const& m ) {
      VW_ASSERT( m.impl().rows() == m_matrix.size1() &&
                 m.impl().cols() == m_matrix.size2(), ArgumentErr() << "Matrix must have dimensions "
                 << m_matrix.size1() << "x" << m_matrix.size2() << "." );
      vw_out(vw::WarningMessage, "math") << "Sparsity destroyed in generic assignment to MatrixSparseSkyline.\n";
      for ( size_t i = 0; i < m.rows(); i++ ) {
        for ( size_t j = 0; j < m.cols(); j++ ) {
          (*this)(i,j) = m(i,j);
        }
      }
      return *this;
    }

    /// Simple Access
    size_t rows() const { return m_matrix.size1(); }
    size_t cols() const { return m_matrix.size2(); }

    /// Change the size of the matrix
    void set_size( size_t new_rows, size_t new_cols, bool preserve = false ) {
      vw_throw( NoImplErr() << "MatrixSparseSkyline::set_size, code has not been written yet." );
    }

    /// Element Access
    reference_type operator()( size_t row, size_t col ) {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(),
                 LogicErr() << "operator() ran off end of matrix" );
#endif
      if ( col > row )
        std::swap(row,col);
      if ( col < m_skyline[row] )
        m_skyline[row] = col;
      return m_matrix(row,col);
    }

    /// Element Access
    const_reference_type operator()( size_t row, size_t col ) const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(),
                 LogicErr() << "operator() ran off end of matrix" );
#endif
      if ( col > row )
        return m_matrix(col,row);
      return m_matrix(row,col);
    }

    /// Pointer Access
    value_type *data() {
      vw_throw( NoImplErr() << "MatrixSparseSkyline::data can not provide direct pointer access.");
      return &(operator()(0,0));
    }
    iterator begin() { return iterator(*this,0,0); }
    const_iterator begin() const { return const_iterator(*this,0,0); }
    iterator end() { return iterator(*this,rows(),0); }
    const_iterator end() const { return const_iterator(*this,rows(),0); }
    sparse_iterator1 sparse_begin() { return m_matrix.begin1(); }
    const_sparse_iterator1 sparse_begin() const { return m_matrix.begin1(); }
    sparse_iterator1 sparse_end() { return m_matrix.end1(); }
    const_sparse_iterator1 sparse_end() const { return m_matrix.end1(); }

    // Special Access to Skyline Vector
    const Vector<size_t>& skyline() const { return m_skyline; }

  };

  /// Dump for MatrixSparseSkyline
  template <class ElemT>
  inline std::ostream& operator<<( std::ostream& os, MatrixSparseSkyline<ElemT> const& m ) {
    // Need to first determine the sampling rate. In bundle adjustment
    // it doesn't make sense to show every element as camera variables will
    // come in blocks of 6.
    size_t smallest_sampling_rate = m.cols();
    size_t curr_sampling_rate = m.cols();
    size_t last_value = 10000;
    const Vector<size_t>& skyline = m.skyline();
    for ( size_t i = 0; i < skyline.size(); i++ ) {
      if ( last_value != skyline(i) ) {
        if ( smallest_sampling_rate  > curr_sampling_rate )
          smallest_sampling_rate = curr_sampling_rate;
        curr_sampling_rate = 0;
      }
      curr_sampling_rate++;
      last_value = skyline(i);
    }

    os << "MatrixSparseSkyline" << m.rows() << "x" << m.cols() << "@" << smallest_sampling_rate << "\n";
    for ( size_t i = 0, ii=0; i < m.rows(); i += smallest_sampling_rate, ii++ ) {
      for ( size_t j = 0, jj=0; j < m.cols(); j += smallest_sampling_rate, jj++ ) {
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

    Vector<size_t> const& skyline = A.skyline();

    for (size_t j = 0; j < A.cols(); ++j) {

      // Compute v(1:j)
      std::vector<double> v(j+1);
      v[j] = A(j,j);
      for (size_t i = skyline(j); i < j; ++i) {
        v[i] = A(j,i)*A(i,i);
        v[j] -= A(j,i)*v[i];
      }

      // Store d(j) and compute L(j+1:n,j)
      A(j,j) = v[j];
      for (size_t i = j+1; i < A.cols(); ++i) {
        double row_sum = 0;
        for (size_t jj = skyline(i); jj < j; ++jj)
          row_sum += A(i,jj)*v[jj];
        if (j >= skyline(i))
          A(i,j) = ( A(i,j)-row_sum ) / v[j];
      }
    }
  }

  // This version excepts an outside skyline matrix
  template <class MatrixT, class VectorT>
  void sparse_ldl_decomposition(MatrixBase<MatrixT>& A, const VectorT& skyline ) {
    for (size_t j = 0; j < A.impl().cols(); ++j) {

      // Compute v(1:j)
      std::vector<double> v(j+1);
      v[j] = A.impl()(j,j);
      for (size_t i = skyline(j); i < j; ++i) {
        v[i] = A.impl()(j,i)*A.impl()(i,i);
        v[j] -= A.impl()(j,i)*v[i];
      }

      // Store d(j) and compute L(j+1:n,j)
      A.impl()(j,j) = v[j];
      for (size_t i = j+1; i < A.impl().cols(); ++i) {
        double row_sum = 0;
        for (size_t jj = skyline(i); jj < j; ++jj)
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
    VW_ASSERT(A.cols() == A.rows(),
              ArgumentErr() << "sparse_solve: matrix must be square and symmetric.\n");

    // Compute the L*D*L^T decomposition of A
    sparse_ldl_decomposition(A);
    Vector<ElemT> x = sparse_solve_ldl(A, b);

    return x;
  }

  /// This version can use an outside skyline matrix
  template <class MatrixT, class VectorT, class VectorSkyT>
  Vector<typename MatrixT::value_type> sparse_solve(MatrixBase<MatrixT>& A,
                                                    VectorT const& b,
                                                    VectorSkyT const& sky ) {
    VW_ASSERT(A.impl().cols() == A.impl().rows(),
              ArgumentErr() << "sparse_solve: matrix must be square and symmetric.\n");

    // Compute the L*D*L^T decomposition of A
    sparse_ldl_decomposition(A.impl(), sky);
    Vector<typename MatrixT::value_type> x = sparse_solve_ldl(A.impl(), b, sky);

    return x;
  }

  //--------------------------------------------------------------
  //            Solve Spare Skyline Linear System: AX=B, where
  //            X and B are matrices
  //--------------------------------------------------------------

  /// Perform L*D*L^T decomposition on a sparse skyline symmetric
  /// semi-definite matrix A (in place) and then call
  /// forward/backward solver

  /// WARNING: Modifies the contents of the matrix A.
  template <class ElemT, class BMatrixT>
  Matrix<typename PromoteType<typename BMatrixT::value_type, typename BMatrixT::value_type>::type>
  multi_sparse_solve(MatrixSparseSkyline<ElemT>& A, BMatrixT & B ) {
    return multi_sparse_solve(A, B, A.skyline());
  }

  template <class AMatrixT, class BMatrixT, class VectorT>
  Matrix<typename PromoteType<typename BMatrixT::value_type, typename BMatrixT::value_type>::type>
    multi_sparse_solve(AMatrixT & A, BMatrixT & B, const VectorT& skyline ) {
    VW_ASSERT(A.cols() == A.rows(), ArgumentErr() << "multi_sparse_solve: matrix must be square and symmetric.\n");
    VW_ASSERT(A.rows() == B.rows(), ArgumentErr() << "multi_sparse_solve: AX=B means A, B have same # of rows.\n");

    Matrix<double> X(B.rows(), B.cols());

    Vector<double> current_col(A.cols());


    // Compute the L*D*L^T decomposition of A
    sparse_ldl_decomposition(A, skyline);

    for(size_t i = 0; i < B.cols(); i++){
      current_col = select_col(B, i);
      select_col(X, i) = sparse_solve_ldl(A, current_col, skyline);
    }
    return X;
  }

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
    return sparse_solve_ldl(A,b,A.skyline());
  }

  /// Version that excepts an outside skyline vector
  template <class MatrixT, class VectorT, class VectorSkyT>
  Vector<typename MatrixT::value_type> sparse_solve_ldl(MatrixBase<MatrixT>& A,
                                                        VectorT const& b,
                                                        VectorSkyT const& sky) {
    MatrixT const& ar = A.impl();

    VW_ASSERT(ar.cols() == ar.rows(),
              ArgumentErr() << "sparse_solve: matrix must be square and symmetric.\n");

    Vector<size_t> inverse_sky(sky.size());

    // Construct the inverse skyline, which is used to optimize the final
    // back substitution step below.
    for (size_t j = 0; j < inverse_sky.size(); ++j) {
      inverse_sky(j) = 0;
      for (ssize_t i = sky.size()-1; i>=0; --i) {
        if (j < sky(i))
          ++(inverse_sky(j));
        else
          break; // Break out of the inner loop
      }
    }

    // Forward Substitution Step ( L*x'=b )
    Vector<typename MatrixT::value_type> x_prime(ar.cols());
    for (size_t i = 0; i < x_prime.size(); ++i) {
      typename MatrixT::value_type sum = 0;
      for (size_t j = sky(i); j < i; ++j)
        sum += ar(i,j)*x_prime(j);
      x_prime(i) = b(i)-sum;
    }

    // Divide by D ( D*x''=x' )
    Vector<typename MatrixT::value_type> x_doubleprime(A.impl().cols());
    for (size_t i = 0; i < x_doubleprime.size(); ++i)
      x_doubleprime(i) = x_prime(i)/ar(i,i);

    // Back Substitution step ( L^T*x=x'' )
    Vector<typename MatrixT::value_type> x(ar.cols());
    for (ssize_t i = x.size()-1; i >= 0; --i) {
      typename MatrixT::value_type sum = 0;
      for (size_t j = i+1; j < ar.cols()-inverse_sky[i]; ++j)
        sum += ar(j,i)*x(j);
      x(i) = x_doubleprime(i) - sum;
    }

    return x;
  }

  //-------------------------------------------------------------------
  // Reorganize Vector and Matrix Class
  //
  // Use to reorganize the order of the elements to reduce the bandwidth
  // of a problem. This is good for sparse matrices in that a solve
  // could be done immensely.
  //-------------------------------------------------------------------

  template <class VectorT>
  class VectorReorganize : public VectorBase<VectorReorganize<VectorT> > {
    VectorT & m_vector;
    std::vector<size_t> m_lookup;

  public:
    typedef typename VectorT::value_type value_type;

    typedef typename VectorT::reference_type reference_type;
    typedef typename VectorT::const_reference_type const_reference_type;

    typedef IndexingVectorIterator<VectorReorganize<VectorT> > iterator;
    typedef IndexingVectorIterator<const VectorReorganize<VectorT> > const_iterator;

    // Constructor
    explicit VectorReorganize( VectorT& vector, std::vector<size_t> const& lookup ) : m_vector(vector), m_lookup(lookup) {
      VW_ASSERT( vector.size()==lookup.size(),
                 ArgumentErr() << "Input Vector and Lookup Chart must have same dimensions" );
    }

    // Access to internal types
    VectorT& child() { return m_vector; }
    VectorT const& child() const { return m_vector; }
    std::vector<size_t>& lookup() { return m_lookup; }
    std::vector<size_t> const& lookup() const { return m_lookup; }
    std::vector<size_t> inverse() const {
      std::vector<size_t> ilookup;
      ilookup.resize( m_lookup.size() );
      for ( size_t i = 0; i < m_lookup.size(); i++ )
        ilookup[m_lookup[i]] = i;
      return ilookup;
    }

    // Standard Properties
    size_t size() const { return m_vector.size(); }
    void set_size( size_t new_size, bool preserve=false ) {
      VW_ASSERT( new_size==size(),
                 ArgumentErr() << "Cannot resize vector reorganize." );
    }

    // Element Access
    reference_type operator()( size_t i ) { return child()[m_lookup[i]]; }
    const_reference_type operator()( size_t i ) const { return child()[m_lookup[i]]; }
    reference_type operator[]( size_t i ) { return child()[m_lookup[i]]; }
    const_reference_type operator[]( size_t i ) const { return child()[m_lookup[i]]; }

    // Pointer Access
    iterator begin() { return iterator(*this,0); }
    const_iterator begin() const { return const_iterator(*this,0); }
    iterator end() { return iterator(*this,size()); }
    const_iterator end() const { return const_iterator(*this,size()); }
  };

  template <class MatrixT>
  class MatrixReorganize : public MatrixBase<MatrixReorganize<MatrixT> > {
    MatrixT & m_matrix;
    std::vector<size_t> m_lookup;

  public:
    typedef typename MatrixT::value_type value_type;

    typedef typename MatrixT::reference_type reference_type;
    typedef typename MatrixT::const_reference_type const_reference_type;

    typedef IndexingMatrixIterator<MatrixReorganize<MatrixT> > iterator;
    typedef IndexingMatrixIterator<const MatrixReorganize<MatrixT> > const_iterator;

    // Constructor
    explicit MatrixReorganize( MatrixT& matrix, std::vector<size_t> const& lookup ) : m_matrix(matrix), m_lookup(lookup) {
      VW_ASSERT( matrix.cols()==lookup.size() &&
                 matrix.rows()==lookup.size(),
                 ArgumentErr() << "Input Matrix must be square, and Lookup Chart must have same dimensions" );
    }

    // Copy assignment operator
    MatrixReorganize& operator=( MatrixReorganize const& m ) {
      VW_ASSERT( m.rows()==rows() && m.cols()==cols(), ArgumentErr() << "Matrix must have dimensions " << rows() << "x" << cols() << "." );
      Matrix<value_type> tmp( m );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class T>
    MatrixReorganize& operator=( MatrixBase<T> const& m ) {
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(), ArgumentErr() << "Matrix must have dimensions " << rows() << "x" << cols() << "." );
      Matrix<value_type> tmp( m );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW matrix expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class T>
    MatrixReorganize& operator=( MatrixNoTmp<T> const& m ) {
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(), ArgumentErr() << "Matrix must have dimensions " << rows() << "x" << cols() << "." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    // Access to internal types
    MatrixT& child() { return m_matrix; }
    MatrixT const& child() const { return m_matrix; }
    std::vector<size_t>& lookup() { return m_lookup; }
    std::vector<size_t> const& lookup() const { return m_lookup; }
    std::vector<size_t> inverse() const {
      std::vector<size_t> ilookup;
      ilookup.resize( m_lookup.size() );
      for ( size_t i = 0; i < m_lookup.size(); i++ )
        ilookup[m_lookup[i]] = i;
      return ilookup;
    }

    // Standard properties
    size_t rows() const { return m_matrix.rows(); }
    size_t cols() const { return m_matrix.cols(); }
    void set_size( size_t new_rows, size_t new_cols, bool preserve=false ) {
      VW_ASSERT( new_rows==rows() && new_cols==cols(),
                 ArgumentErr() << "Cannot resize matrix reorganize." );
    }

    // Element Access
    reference_type operator()( size_t row, size_t col ) {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      return child()(m_lookup[row],m_lookup[col]);
    }
    const_reference_type operator()( size_t row, size_t col ) const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      return child()(m_lookup[row],m_lookup[col]);
    }

    // Pointer Access
    IndexingMatrixIterator<MatrixReorganize<MatrixT> > begin() { return iterator(*this,0,0); }
    IndexingMatrixIterator<const MatrixReorganize<MatrixT> > begin() const { return const_iterator(*this,0,0); }
    IndexingMatrixIterator<MatrixReorganize<MatrixT> > end() { return iterator(*this,rows(),0); }
    IndexingMatrixIterator<const MatrixReorganize<MatrixT> > end() const { return const_iterator(*this,rows(),0); }
  };

  // User ease functions
  template <class VectorT>
  inline VectorReorganize<VectorT> reorganize( VectorBase<VectorT>& v,
                                               std::vector<size_t>& lookup ) {
    return VectorReorganize<VectorT>(v.impl(), lookup);
  }

  template <class VectorT>
  inline VectorReorganize<const VectorT> reorganize( VectorBase<VectorT> const& v,
                                                     std::vector<size_t> const& lookup ) {
    return VectorReorganize<const VectorT>(v.impl(), lookup);
  }

  template <class MatrixT>
  inline MatrixReorganize<MatrixT> reorganize( MatrixBase<MatrixT>& m,
                                               std::vector<size_t>& lookup ) {
    return MatrixReorganize<MatrixT>(m.impl(), lookup);
  }

  template <class MatrixT>
  inline MatrixReorganize<const MatrixT> reorganize( MatrixBase<MatrixT> const& m,
                                                     std::vector<size_t> const& lookup ) {
    return MatrixReorganize<const MatrixT>(m.impl(), lookup);
  }

  //------------------------------------------------------------------
  // Cuthill McKee Organization Solver
  //
  // Use to solve for a organization of the sparse matrix that best
  // reduces the bandwidth of the matrix. Reducing the bandwidth
  // can greatly increase the speed at which a sparse LDL can be solved
  // for.
  //------------------------------------------------------------------

  template <class ElemT>
  std::vector<size_t> cuthill_mckee_ordering(MatrixSparseSkyline<ElemT>& A,
                                             size_t sampling_rate ) {
    // First Working out the Sampling Rate (cheat to save time in
    // Bundle Adjustment where sampling rate is the number of camera
    // parameters)
    const Vector<size_t>& skyline = A.skyline();

    // Boost Graph Definitions
    typedef boost::adjacency_list<boost::vecS,boost::vecS, boost::undirectedS,
      boost::property<boost::vertex_color_t, boost::default_color_type,
      boost::property<boost::vertex_degree_t, size_t> > > Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::vertices_size_type size_type;

    // Finding Connections
    Graph G( A.cols() );
    for ( size_t i = 0; i < A.rows(); i += sampling_rate )
      for ( size_t j = skyline[i]; j < i; j += sampling_rate )
        if ( A(i,j) != 0 )
          for ( size_t ii = 0; ii < sampling_rate; ii++ )
            for ( size_t jj = 0; jj < sampling_rate; jj++ )
              boost::add_edge( i+ii, j+jj, G);

    boost::property_map<Graph,boost::vertex_index_t>::type
      index_map = get(boost::vertex_index, G);

    vw_out(DebugMessage,"math") << "-> CutHill McKee starting B: "
                                << boost::bandwidth(G) << "\n";

    // Solving for CutHill McKee
    std::vector<Vertex> inv_perm(boost::num_vertices(G));
    std::vector<size_type> perm(boost::num_vertices(G));
    cuthill_mckee_ordering(G, inv_perm.rbegin(), get(boost::vertex_color,G),
                           make_degree_map(G));

    // Building new lookup chart
    std::vector<size_t> lookup_chart( A.cols() );
    for ( size_t i = 0; i < inv_perm.size(); i++ )
      lookup_chart[i] = index_map[inv_perm[i]];

    // Finding new bandwidth for debug purposes
    for (size_type c = 0; c != inv_perm.size(); ++c )
      perm[index_map[inv_perm[c]]] = c;
    vw_out(DebugMessage,"math") << "-> CutHill McKee ending B: "
                                << boost::bandwidth(G, make_iterator_property_map(&perm[0], index_map, perm[0])) << "\n";

    return lookup_chart;
  }

  // --------------------------------------------------------------
  // Function to aid in solving for new skyline vector after a
  // reordering has happened from CutHill Mckee
  // --------------------------------------------------------------

  template <class MatrixT>
  Vector<size_t> solve_for_skyline(MatrixBase<MatrixT> const& A) {
    MatrixT const& ar = A.impl();
    size_t rows = ar.rows();
    Vector<size_t> skyline(rows);
    for ( size_t i = 0; i < rows; i++ ) {
      size_t j = 0;
      while ( j < i && ar(i,j) == 0 )
        j++;
      skyline[i] = j;
    }
    return skyline;
  }

}} // namespace vw::math

#endif//__VW_MATH_SPARSE_SKYLINE_MATRIX_H__
