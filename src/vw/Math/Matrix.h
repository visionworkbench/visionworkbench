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

/// \file Matrix.h
///
/// Provides the core Matrix<> mathematical matrix template classes.  
///
/// Currently I believe we support:
///   Element access via m(i,j) and m[i][j]
///   Seamless conversion between compatibile matrix types
///   Matrix addition, subtraction, and negation
///   Scalar multiplication and division
///   Matrix*matrix, matrix*vector, and vector*matrix products
///   Sum of elements via sum()
///   Trace via trace()
///   Transpose via transpose()
///   Inverse via inverse(), not particularly robust to singularity
///   Matrix norms via norm_1(), norm_inf(), and norm_frobenius()
///   Scalar addition and subtraction via elem_sum() and elem_dif()
///   Elementwise multiplication and division via elem_prod() and elem_quot()
///
/// We intentionally do *not* define matrix division using inverse(),
/// since users should probably avoid using inverse() anyway unless
/// they know what they are doing.  It's easy enough and preferable
/// for them to write an extra line of code to be explicit about what
/// they want, using e.g. the SVD.
///
#ifndef __VW_MATH_MATRIX_H__
#define __VW_MATH_MATRIX_H__

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/result_of.hpp>

#include <vw/Core/Exception.h>
#include <vw/Math/Vector.h>

namespace vw {
namespace math {

  namespace bnu = boost::numeric::ublas;

  /// A type function to compute the number of rows of a matrix
  /// expression at compile time (or zero for dynamically-sized
  /// vectors).
  template <class MatrixT>
  struct MatrixRows {
    const static int value = 0;
  };

  /// A type function to compute the number of columns of a matrix
  /// expression at compile time (or zero for dynamically-sized
  /// vectors).
  template <class MatrixT>
  struct MatrixCols {
    const static int value = 0;
  };


  // *******************************************************************
  // class MatrixBase<MatrixT>
  // A CRTP matrix base class.
  // *******************************************************************

  /// A CRTP base class for matrices and matrix expressions.  
  /// Provides a mechanism for restricting function arguments to 
  /// matrices, and provides the various arithmetic assignment 
  /// operators and other member functions.
  template <class MatrixT>
  struct MatrixBase {

    /// Returns the derived implementation type.
    MatrixT& impl() { return *static_cast<MatrixT*>(this); }

    /// Returns the derived implementation type.
    MatrixT const& impl() const { return *static_cast<MatrixT const*>(this); }

    /// Sum-assignment operator
    template <class T>
    MatrixT& operator+=( T const& m ) {
      return impl() = impl() + m;
    }

    /// Difference-assignment operator
    template <class T>
    MatrixT& operator-=( T const& m ) {
      return impl() = impl() - m;
    }

    /// Product-assignment operator
    template <class T>
    MatrixT& operator*=( T s ) {
      return impl() = impl() * s;
    }

    /// Quotient-assignment operator
    template <class T>
    MatrixT& operator/=( T s ) {
      return impl() = impl() / s;
    }

    /// Set the matrix to the identity matrix with the same dimensions.
    void set_identity() {
      VW_ASSERT( impl().rows()==impl().cols(), LogicErr() << "Only square matrices can be identity matrices." );
      int n=impl().rows();
      for( int i=0; i<n; ++i )
        for( int j=0; j<n; ++j )
          impl()(i,j)=(i==j)?(typename MatrixT::value_type(1)):(typename MatrixT::value_type(0));
    }

    /// Set the matrix to the identity matrix with the given dimensions.
    void set_identity( unsigned size ) {
      impl().set_size( size, size );
      set_identity();
    }

  protected:
    /// Okay, this is pretty clearly the wrong place for this, but
    /// here it is for now.
    template<class MT>
    static inline bnu::matrix<typename MT::value_type> bnu_matrix_inverse( bnu::matrix_expression<MT> const& m ) {
      typedef typename MT::value_type value_type;
      unsigned size = m().size1();
      bnu::matrix<value_type> buf( m() );
      bnu::permutation_matrix<unsigned> pm(size);
      lu_factorize(buf,pm);
      bnu::matrix<value_type> inverse = bnu::identity_matrix<value_type>(size);
      lu_substitute(buf, pm, inverse);
      return inverse;
    }
  };


  // *******************************************************************
  // class IndexingMatrixIterator<MatrixT>
  // A general-purpose matrix iterator type.
  // *******************************************************************

  /// An iterator for an arbitrary matrix type that iterates over the
  /// elements of the matrix in the standard (row-major) order.  It
  /// keeps track of the element indices, dereferencing via the 
  /// function call operator.
  template <class MatrixT>
  class IndexingMatrixIterator : public boost::iterator_facade<IndexingMatrixIterator<MatrixT>, 
                                                               typename boost::mpl::if_<boost::is_const<MatrixT>,
                                                                                        const typename MatrixT::value_type,
                                                                                        typename MatrixT::value_type>::type,
                                                               boost::random_access_traversal_tag,
                                                               typename boost::mpl::if_<boost::is_const<MatrixT>,
                                                                                        typename MatrixT::const_reference_type,
                                                                                        typename MatrixT::reference_type>::type>
  {
    friend class boost::iterator_core_access;
      
    MatrixT& m_matrix;
    unsigned m_row, m_col;
      
    bool equal( IndexingMatrixIterator const& iter ) const {
      return m_row==iter.m_row && m_col==iter.m_col;
    }

    ptrdiff_t distance_to( IndexingMatrixIterator const &iter ) const {
      ptrdiff_t coldiff = (iter.m_col>m_col) ? ptrdiff_t(iter.m_col-m_col) : -ptrdiff_t(m_col-iter.m_col);
      ptrdiff_t rowdiff = (iter.m_row>m_row) ? ptrdiff_t(iter.m_row-m_row) : -ptrdiff_t(m_row-iter.m_row);
      return coldiff + rowdiff * m_matrix.cols();
    }

    void increment() {
      if( ++m_col == m_matrix.cols() ) {
        m_col=0; ++m_row;
      }
    }

    void decrement() {
      if( m_col==0 ) {
        m_col=m_matrix.cols()-1;
        --m_row;
      }
      else {
        --m_col;
      }
    }
    
    void advance( ptrdiff_t n ) {
      if( n < 0 ) {
        ptrdiff_t rowdiff = 1 + (-n)/m_matrix.cols();
        m_row -= rowdiff;
        n += rowdiff * m_matrix.cols();
      }
      m_col += n;
      m_row += m_col / m_matrix.cols();
      m_col %= m_matrix.cols();
    }

    typename IndexingMatrixIterator::reference dereference() const {
      return m_matrix(m_row,m_col);
    }

  public:
    IndexingMatrixIterator( MatrixT& matrix, ptrdiff_t row, ptrdiff_t col ) : 
      m_matrix(matrix), m_row(row), m_col(col) {}
  };


  // *******************************************************************
  // class Matrix<ElemT,RowsN,ColsN>
  // A statically-allocated fixed-dimension matrix class.
  // *******************************************************************

  /// A fixed-dimension mathematical matrix class.
  template <class ElemT, int RowsN=0, int ColsN=0>
  class Matrix : public MatrixBase<Matrix<ElemT,RowsN,ColsN> >
  {
    typedef bnu::matrix<ElemT,bnu::row_major,FixedArray<ElemT,RowsN*ColsN> > core_type;
    core_type core_;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef ElemT* iterator;
    typedef ElemT const* const_iterator;

    /// Constructs a matrix of zeroes.
    Matrix() : core_(RowsN,ColsN) {
      core_.clear();
    }

    /// Constructs a matrix whose first element is as given.
    Matrix( ElemT e1 ) : core_(RowsN,ColsN) {
      BOOST_STATIC_ASSERT( RowsN*ColsN >= 1 );
      iterator i=begin();
      *(i++)=e1;
      for( ; i!=end(); ++i ) *i = ElemT();
    }

    /// Constructs a matrix whose first two elements are as given.
    Matrix( ElemT e1, ElemT e2 ) : core_(RowsN,ColsN) {
      BOOST_STATIC_ASSERT( RowsN*ColsN >= 2 );
      iterator i=begin();
      *(i++)=e1; *(i++)=e2;
      for( ; i!=end(); ++i ) *i = ElemT();
    }

    /// Constructs a matrix whose first three elements are as given.
    Matrix( ElemT e1, ElemT e2, ElemT e3 ) : core_(RowsN,ColsN) {
      BOOST_STATIC_ASSERT( RowsN*ColsN >= 3 );
      iterator i=begin();
      *(i++)=e1; *(i++)=e2; *(i++)=e3;
      for( ; i!=end(); ++i ) *i = ElemT();
    }

    /// Constructs a matrix whose first four elements are as given.
    Matrix( ElemT e1, ElemT e2, ElemT e3, ElemT e4 ) : core_(RowsN,ColsN) {
      BOOST_STATIC_ASSERT( RowsN*ColsN >= 4 );
      iterator i=begin();
      *(i++)=e1; *(i++)=e2; *(i++)=e3; *(i++)=e4;
      for( ; i!=end(); ++i ) *i = ElemT();
    }

    /// Constructs a matrix from given densely-packed row-mjor data.
    /// This constructor copies the data.  If you wish to make a
    /// shallow proxy object instead, see vw::MatrixProxy.
    Matrix( const ElemT data[RowsN*ColsN] ) : core_(RowsN,ColsN) {
      std::copy( data, data+RowsN*ColsN, core_.data().begin() );
    }

    /// Standard copy constructor.
    Matrix( Matrix const& m ) : core_( m.core_ ) {}

    /// Generalized copy constructor, from arbitrary VW matrix expressions.
    template <class T>
    Matrix( MatrixBase<T> const& m ) : core_(RowsN,ColsN) {
      VW_ASSERT( m.impl().rows()==RowsN && m.impl().cols()==ColsN, ArgumentErr() << "Matrix must have dimensions " << RowsN << "x" << ColsN << "." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
    }

    /// Generalized copy constructor, from arbitrary uBLAS matrix expressions.
    template <class T>
    Matrix( bnu::matrix_expression<T> const& m ) : core_(RowsN,ColsN) {
      VW_ASSERT( m().size1()==RowsN && m().size2()==ColsN, ArgumentErr() << "Matrix must have dimensions " << RowsN << "x" << ColsN << "." );
      core_.assign(m);
    }

    /// Standard copy assignment operator.
    Matrix& operator=( Matrix const& m ) {
      core_ = m.core_;
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class T>
    Matrix& operator=( MatrixBase<T> const& m ) { 
      VW_ASSERT( m.impl().rows()==RowsN && m.impl().cols()==ColsN, ArgumentErr() << "Matrix must have dimensions " << RowsN << "x" << ColsN << "." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary uBLAS matrix expressions.
    template <class T>
    Matrix& operator=( bnu::matrix_expression<T> const& m ) { 
      VW_ASSERT( m().size1()==RowsN && m().size2()==ColsN, ArgumentErr() << "Matrix must have dimensions " << RowsN << "x" << ColsN << "." );
      core_.assign( m );
      return *this;
    }

    /// Returns the number of rows in the matrix.
    unsigned rows() const { return RowsN; }

    /// Returns the number of columns in the matrix.
    unsigned cols() const { return ColsN; }

    /// Change the size of the matrix.  Elements in memory are preserved when specified.
    void set_size( unsigned new_rows, unsigned new_cols, bool preserve = false ) {
      VW_ASSERT( new_rows==rows() && new_cols==cols(),
                 ArgumentErr() << "Cannot resize a fixed-size Matrix." );
    }

    /// Access an element
    value_type& operator()( unsigned row, unsigned col ) {
      return core_(row,col);
    }

    /// Access an element
    value_type const& operator()( unsigned row, unsigned col ) const {
      return core_(row,col);
    }

    /// Access an individual matrix row, for further access using a second operator[].
    bnu::matrix_row<core_type> operator[]( unsigned row ) {
      return bnu::matrix_row<core_type>( core_, row );
    }

    /// Access an individual matrix row, for further access using a second operator[].
    bnu::matrix_row<const core_type> operator[]( unsigned row ) const {
      return bnu::matrix_row<const core_type>( core_, row );
    }

    value_type *data() {
      return &(operator()(0,0));
    }

    const value_type *data() const {
      return &(operator()(0,0));
    }

    iterator begin() { 
      return &(core_(0,0));
    }

    const_iterator begin() const { 
      return &(core_(0,0));
    }

    iterator end() {
      return &(core_(0,0)) + RowsN*ColsN;
    }

    const_iterator end() const {
      return &(core_(0,0)) + RowsN*ColsN;
    }

    Matrix inverse() const {
      return bnu_matrix_inverse( core_ );
    }
  };

  template <class ElemT, int RowsN, int ColsN>
  struct MatrixRows<Matrix<ElemT,RowsN,ColsN> > {
    static const int value = RowsN;
  };

  template <class ElemT, int RowsN, int ColsN>
  struct MatrixCols<Matrix<ElemT,RowsN,ColsN> > {
    static const int value = ColsN;
  };


  // *******************************************************************
  // class Matrix<ElemT>
  // A dynamically-allocated arbitrary-dimension matrix class.
  // *******************************************************************

  /// An arbitrary-dimension mathematical matrix class.
  template <class ElemT>
  class Matrix<ElemT,0,0> : public MatrixBase<Matrix<ElemT> > {
    typedef bnu::matrix<ElemT> core_type;
    core_type core_;
    friend class MatrixImplementation;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef ElemT* iterator;
    typedef ElemT const* const_iterator;

    /// Constructs a matrix with zero size.
    Matrix() {}

    /// Constructs a zero matrix of the given size.
    Matrix( unsigned rows, unsigned cols ) : core_(rows,cols) {
      core_.clear();
    }

    /// Constructs a matrix of the given size from given
    /// densely-packed row-mjor data.  This constructor copies the
    /// data.  If you wish to make a shallow proxy object instead, 
    /// see vw::MatrixProxy.
    Matrix( unsigned rows, unsigned cols, const ElemT *data ) : core_(rows,cols) {
      std::copy( data, data+rows*cols, core_.data().begin() );
    }

    /// Standard copy constructor.
    Matrix( Matrix const& m ) : core_( m.core_ ) {}

    /// Generalized copy constructor, from arbitrary VW matrix expressions.
    template <class T>
    Matrix( MatrixBase<T> const& m ) : core_(m.impl().rows(),m.impl().cols()) {
      std::copy( m.impl().begin(), m.impl().end(), begin() );
    }

    /// Generalized copy constructor, from arbitrary uBLAS matrix expressions.
    template <class T>
    Matrix( bnu::matrix_expression<T> const& m ) : core_(m) {}

    /// Standard copy assignment operator.
    Matrix& operator=( Matrix const& m ) {
      core_ = m.core_;
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class T>
    Matrix& operator=( MatrixBase<T> const& m ) {
      set_size( m.impl().rows(), m.impl().cols() );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary uBLAS matrix expressions.
    template <class T>
    Matrix& operator=( bnu::matrix_expression<T> const& m ) {
      core_ = m;
      return *this;
    } 

    /// Returns the number of rows in the matrix.
    unsigned rows() const { return core_.size1(); }

    /// Returns the number of columns in the matrix.
    unsigned cols() const { return core_.size2(); }

    /// Change the size of the matrix.  Elements in memory are preserved when specified.
    void set_size( unsigned new_rows, unsigned new_cols, bool preserve = false ) {
      core_.resize(new_rows, new_cols, preserve);
    }

    /// Access an element
    value_type& operator()( unsigned row, unsigned col ) {
      return core_(row,col);
    }

    /// Access an element
    value_type const& operator()( unsigned row, unsigned col ) const {
      return core_(row,col);
    }

    /// Access an individual matrix row, for further access using a second operator[].
    bnu::matrix_row<core_type> operator[]( unsigned row ) {
      return bnu::matrix_row<core_type>( core_, row );
    }

    /// Access an individual matrix row, for further access using a second operator[].
    bnu::matrix_row<const core_type> operator[]( unsigned row ) const {
      return bnu::matrix_row<const core_type>( core_, row );
    }

    value_type *data() {
      return &(operator()(0,0));
    }

    const value_type *data() const {
      return &(operator()(0,0));
    }

    iterator begin() { 
      return &(core_(0,0));
    }

    const_iterator begin() const { 
      return &(core_(0,0));
    }

    iterator end() {
      return &(core_(0,0)) + rows()*cols();
    }

    const_iterator end() const {
      return &(core_(0,0)) + rows()*cols();
    }

    Matrix inverse() const {
      return bnu_matrix_inverse( core_ );
    }
  };


  // *******************************************************************
  // class MatrixProxy<ElemT,RowsN,ColsN>
  // A fixed-dimension matrix proxy class, treating an arbitrary block 
  // of memory as a Matrix in packed row-major format.
  // *******************************************************************

  /// A fixed-dimension mathematical matrix class.
  template <class ElemT, int RowsN=0, int ColsN=0>
  class MatrixProxy : public MatrixBase<MatrixProxy<ElemT,RowsN,ColsN> >
  {
    ElemT *m_ptr;
    friend class MatrixImplementation;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef ElemT* iterator;
    typedef ElemT const* const_iterator;

    /// Constructs a matrix proxy
    MatrixProxy( ElemT* ptr ) : m_ptr(ptr) {}

    /// Standard copy assignment operator.
    MatrixProxy& operator=( MatrixProxy const& m ) {
      VW_ASSERT( m.rows()==rows() && m.cols()==cols(),
                 ArgumentErr() << "Matrix must have dimensions " << rows() 
                 << "x" << cols() << " in matrix proxy assignment." );
      std::copy( m.begin(), m.end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class T>
    MatrixProxy& operator=( MatrixBase<T> const& m ) { 
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(),
                 ArgumentErr() << "Matrix must have dimensions " << rows() 
                 << "x" << cols() << " in matrix proxy assignment." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary uBLAS matrix expressions.
    template <class T>
    MatrixProxy& operator=( bnu::matrix_expression<T> const& m ) {
      VW_ASSERT( m.impl().size1()==rows() && m.impl().size2()==cols(),
                 ArgumentErr() << "Matrix must have dimensions " << rows() 
                 << "x" << cols() << " in matrix proxy assignment." );
      std::copy( m.begin(), m.end(), begin() );
      return *this;
    } 

    /// Returns the number of rows in the matrix.
    unsigned rows() const { return RowsN; }

    /// Returns the number of columns in the matrix.
    unsigned cols() const { return ColsN; }

    /// Change the size of the matrix.
    /// Elements in memory are preserved when specified.
    void set_size( unsigned new_rows, unsigned new_cols, bool preserve = false ) {
      VW_ASSERT( new_rows==rows() && new_cols==cols(),
                 ArgumentErr() << "Cannot resize matrix proxy." );
    }

    /// Access an element
    value_type& operator()( unsigned row, unsigned col ) {
      return m_ptr[row*cols()+col];
    }

    /// Access an element
    value_type const& operator()( unsigned row, unsigned col ) const {
      return m_ptr[row*cols()+col];
    }

    value_type *data() {
      return &(operator()(0,0));
    }

    const value_type *data() const {
      return &(operator()(0,0));
    }

    iterator begin() { 
      return m_ptr;
    }

    const_iterator begin() const { 
      return m_ptr;
    }

    iterator end() {
      return m_ptr+rows()*cols();
    }

    const_iterator end() const {
      return m_ptr+rows()*cols();
    }

  };

  template <class ElemT, int RowsN, int ColsN>
  struct MatrixRows<MatrixProxy<ElemT,RowsN,ColsN> > {
    static const int value = RowsN;
  };

  template <class ElemT, int RowsN, int ColsN>
  struct MatrixCols<MatrixProxy<ElemT,RowsN,ColsN> > {
    static const int value = ColsN;
  };


  // *******************************************************************
  // class MatrixProxy<ElemT>
  // A arbitrary-dimension matrix proxy class, treating an arbitrary  
  // block of memory as a Matrix in packed row-major format.
  // *******************************************************************

  /// An arbitrary-dimension matrix proxy class.
  template <class ElemT>
  class MatrixProxy<ElemT,0,0> : public MatrixBase<MatrixProxy<ElemT> > {
    ElemT *m_ptr;
    unsigned m_rows, m_cols;
    friend class MatrixImplementation;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef ElemT* iterator;
    typedef ElemT const* const_iterator;

    /// Constructs a matrix with zero size.
    MatrixProxy( ElemT* ptr, unsigned rows, unsigned cols )
      : m_ptr(ptr), m_rows(rows), m_cols(cols) {}

    template <class ContainerT>
    MatrixProxy( ContainerT const& container )
      : m_ptr(container.data()), m_rows(container.rows()), m_cols(container.cols()) {}

    /// Standard copy assignment operator.
    MatrixProxy& operator=( MatrixProxy const& m ) {
      VW_ASSERT( m.rows()==rows() && m.cols()==cols(),
                 ArgumentErr() << "Matrix must have dimensions " << rows() 
                 << "x" << cols() << " in matrix proxy assignment." );
      std::copy( m.m_ptr, m.m_ptr+m_rows*m_cols, m_ptr );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class T>
    MatrixProxy& operator=( MatrixBase<T> const& m ) { 
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(),
                 ArgumentErr() << "Matrix must have dimensions " << rows() 
                 << "x" << cols() << " in matrix proxy assignment." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary uBLAS matrix expressions.
    template <class T>
    MatrixProxy& operator=( bnu::matrix_expression<T> const& m ) {
      VW_ASSERT( m.impl().size1()==rows() && m.impl().size2()==cols(),
                 ArgumentErr() << "Matrix must have dimensions " << rows() 
                 << "x" << cols() << " in matrix proxy assignment." );
      std::copy( m.begin(), m.end(), begin() );
      return *this;
    } 

    /// Returns the number of rows in the matrix.
    unsigned rows() const { return m_rows; }

    /// Returns the number of columns in the matrix.
    unsigned cols() const { return m_cols; }

    /// Change the size of the matrix.
    /// Elements in memory are preserved when specified.
    void set_size( unsigned new_rows, unsigned new_cols, bool preserve = false ) {
      VW_ASSERT( new_rows==rows() && new_cols==cols(),
                 ArgumentErr() << "Cannot resize matrix proxy." );
    }

    /// Access an element
    value_type& operator()( unsigned row, unsigned col ) {
      return m_ptr[row*m_cols+col];
    }

    /// Access an element
    value_type const& operator()( unsigned row, unsigned col ) const {
      return m_ptr[row*m_cols+col];
    }

    value_type *data() {
      return &(operator()(0,0));
    }

    const value_type *data() const {
      return &(operator()(0,0));
    }

    iterator begin() { 
      return m_ptr;
    }

    const_iterator begin() const { 
      return m_ptr;
    }

    iterator end() {
      return m_ptr+rows()*cols();
    }

    const_iterator end() const {
      return m_ptr+rows()*cols();
    }

  };

  /// Shallow view of an image as a matrix.  Returns a MatrixProxy
  /// object with the same element type as the channel type in the
  /// original image.
  template <class ContainerT>
  MatrixProxy<typename ContainerT::value_type> 
  matrix_proxy( ContainerT const& container ) {
    return MatrixProxy<typename ContainerT::value_type>( container );
  }

  /// Shallow view of an image as a matrix.  Returns a MatrixProxy
  /// object with the same element type as the channel type in the
  /// original image.
  template <class DataT>
  MatrixProxy<DataT> 
  matrix_proxy( DataT* data_ptr, int rows, int cols) {
    return MatrixProxy<DataT>( data_ptr, rows, cols );
  }


  // *******************************************************************
  // class VectorTranspose<VectorT>
  // A transposed matrix class.
  // *******************************************************************

  /// A matrix transpose class.
  template <class MatrixT>
  class MatrixTranspose : public MatrixBase<MatrixTranspose<MatrixT> > {
    MatrixT &m_matrix;
    friend class MatrixImplementation;

  public:
    typedef typename MatrixT::value_type value_type;

    typedef typename MatrixT::reference_type reference_type;
    typedef typename MatrixT::const_reference_type const_reference_type;

    typedef IndexingMatrixIterator<MatrixTranspose<MatrixT> > iterator;
    typedef IndexingMatrixIterator<const MatrixTranspose<MatrixT> > const_iterator;

    /// Constructs a matrix transpose.
    MatrixTranspose( MatrixT& matrix ) : m_matrix(matrix) {}

    /// Standard assignment operator.
    MatrixTranspose& operator=( MatrixTranspose const& m ) {
      VW_ASSERT( m.rows()==rows() && m.cols()==cols(),
                 ArgumentErr() << "Matrix must have dimensions " << rows() 
                 << "x" << cols() << " in matrix transpose assignment." );
      std::copy( m.m_matrix.begin(), m.m_matrix.end(), m_matrix.begin() );
      return *this;
    }

    /// Generalized assignment operator.
    template <class T>
    MatrixTranspose& operator=( MatrixBase<T> const& m ) { 
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(),
                 ArgumentErr() << "Matrix must have dimensions " << rows() 
                 << "x" << cols() << " in matrix transpose assignment." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    /// Returns the underlying non-transposed matrix.
    MatrixT& child() { return m_matrix; }
    
    /// Returns the underlying non-transposed matrix (const overload).
    MatrixT const& child() const { return m_matrix; }

    /// Returns the number of rows in the matrix.
    unsigned rows() const { return m_matrix.cols(); }

    /// Returns the number of columns in the matrix.
    unsigned cols() const { return m_matrix.rows(); }

    /// Change the size of the matrix.
    /// Elements in memory are preserved when specified.
    void set_size( unsigned new_rows, unsigned new_cols, bool preserve = false ) {
      VW_ASSERT( new_rows==rows() && new_cols==cols(),
                 ArgumentErr() << "Cannot resize matrix transpose." );
    }

    /// Access an element
    reference_type operator()( unsigned row, unsigned col ) {
      return m_matrix(col,row);
    }

    /// Access an element
    const_reference_type operator()( unsigned row, unsigned col ) const {
      return m_matrix(col,row);
    }

    iterator begin() { 
      return iterator(*this,0,0);
    }

    const_iterator begin() const { 
      return const_iterator(*this,0,0);
    }

    iterator end() {
      return iterator(*this,rows(),0);
    }

    const_iterator end() const {
      return const_iterator(*this,rows(),0);
    }
  };

  /// Matrix transpose.
  template <class MatrixT>
  inline MatrixTranspose<MatrixT> transpose( MatrixBase<MatrixT>& m ) {
    return MatrixTranspose<MatrixT>( m.impl() );
  }

  /// Matrix transpose (const overload).
  template <class MatrixT>
  inline MatrixTranspose<const MatrixT> transpose( MatrixBase<MatrixT> const& m ) {
    return MatrixTranspose<const MatrixT>( m.impl() );
  }

  /// Matrix transpose (transpose overload).
  template <class MatrixT>
  inline MatrixT& transpose( MatrixTranspose<MatrixT>& m ) {
    return m.child();
  }

  /// Matrix transpose (const transpose overload).
  template <class MatrixT>
  inline MatrixT const& transpose( MatrixTranspose<MatrixT> const& m ) {
    return m.child();
  }


  // *******************************************************************
  // class MatrixRow<MatrixT>
  // A matrix row class with vector semantics.
  // *******************************************************************

  /// A matrix row reference object.
  template <class MatrixT>
  class MatrixRow : public VectorBase<MatrixRow<MatrixT> > {
    MatrixT& m;
    unsigned row;
  public:
    typedef typename MatrixT::value_type value_type;

    typedef typename boost::mpl::if_<boost::is_const<MatrixT>,
                                     typename MatrixT::const_reference_type,
                                     typename MatrixT::reference_type>::type reference_type;
    typedef typename MatrixT::const_reference_type const_reference_type;

    typedef typename boost::mpl::if_<boost::is_const<MatrixT>,
                                     typename MatrixT::const_iterator,
                                     typename MatrixT::iterator>::type iterator;
    typedef typename MatrixT::const_iterator const_iterator;

    MatrixRow( MatrixT& m, unsigned row ) : m(m), row(row) {}
    
    MatrixRow& operator=( MatrixRow const& v ) {
      VW_ASSERT( v.size()==size(), 
                 ArgumentErr() << "Vectors must have same size in matrix row assignment" );
      std::copy( v.begin(), v.end(), begin() );
      return *this;
    }

    template <class OtherT>
    MatrixRow& operator=( VectorBase<OtherT> const& v ) {
      VW_ASSERT( v.impl().size()==size(),
                 ArgumentErr() << "Vectors must have same size in matrix row assignment" );
      std::copy( v.impl().begin(), v.impl().end(), begin() );
      return *this;
    }

    unsigned size() const {
      return m.cols();
    }

    reference_type operator()( int i ) {
      return m(row,i);
    }

    const_reference_type operator()( int i ) const {
      return m(row,i);
    }

    iterator begin() {
      return m.begin() + row*m.cols();
    }

    const_iterator begin() const {
      return m.begin() + row*m.cols();
    }

    iterator end() {
      return m.begin() + (row+1)*m.cols();
    }

    const_iterator end() const {
      return m.begin() + (row+1)*m.cols();
    }
  };

  template <class MatrixT>
  struct VectorSize<MatrixRow<MatrixT> > {
    static const int value = MatrixCols<MatrixT>::value;
  };

  /// Extract a row of a matrix as a vector.
  template <class MatrixT>
  inline MatrixRow<MatrixT> select_row( MatrixBase<MatrixT>& matrix, unsigned row ) {
    return MatrixRow<MatrixT>( matrix.impl(), row );
  }

  /// Extract a row of a matrix as a vector (const overload).
  template <class MatrixT>
  inline MatrixRow<const MatrixT> select_row( MatrixBase<MatrixT> const& matrix, unsigned row ) {
    return MatrixRow<const MatrixT>( matrix.impl(), row );
  }


  // *******************************************************************
  // class MatrixRow<MatrixT>
  // A matrix column class with vector semantics.
  // *******************************************************************

  /// A matrix column reference object.
  template <class MatrixT>
  class MatrixCol : public VectorBase<MatrixCol<MatrixT> >
  {
    MatrixT& m;
    unsigned col;

    template <class IterT>
    class Iterator : public boost::iterator_facade<Iterator<IterT>,
                                                   typename std::iterator_traits<IterT>::value_type,
                                                   boost::random_access_traversal_tag,
                                                   typename std::iterator_traits<IterT>::reference,
                                                   typename std::iterator_traits<IterT>::difference_type>
    {
      friend class boost::iterator_core_access;
      
      IterT i;
      typename Iterator::difference_type stride;
      
      bool equal( Iterator const& iter ) const { return i==iter.i; }
      typename Iterator::difference_type distance_to( Iterator const &iter ) const { return (iter.i - i) / stride; }
      void increment() { i += stride; }
      void decrement() { i -= stride; }
      void advance( ptrdiff_t n ) { i += n*stride; }
      typename Iterator::reference dereference() const { return *i; }
    public:
      Iterator( IterT const& i, typename Iterator::difference_type stride ) : i(i), stride(stride) {}
    };

  public:
    typedef typename MatrixT::value_type value_type;

    typedef typename boost::mpl::if_<boost::is_const<MatrixT>,
                                     typename MatrixT::const_reference_type,
                                     typename MatrixT::reference_type>::type reference_type;
    typedef typename MatrixT::const_reference_type const_reference_type;

    typedef typename boost::mpl::if_<boost::is_const<MatrixT>,
                                     Iterator<typename MatrixT::const_iterator>,
                                     Iterator<typename MatrixT::iterator> >::type iterator;
    typedef Iterator<typename MatrixT::const_iterator> const_iterator;

    MatrixCol( MatrixT& m, unsigned col ) : m(m), col(col) {}
    
    MatrixCol& operator=( MatrixCol const& v ) {
      VW_ASSERT( v.size()==size(), 
                 ArgumentErr() << "Vectors must have same size in matrix column assignment" );
      std::copy( v.begin(), v.end(), begin() );
      return *this;
    }

    template <class OtherT>
    MatrixCol& operator=( VectorBase<OtherT> const& v ) {
      VW_ASSERT( v.impl().size()==size(),
                 ArgumentErr() << "Vectors must have same size in matrix column assignment" );
      std::copy( v.impl().begin(), v.impl().end(), begin() );
      return *this;
    }

    unsigned size() const {
      return m.rows();
    }

    reference_type operator()( int i ) {
      return m(i,col);
    }

    const_reference_type operator()( int i ) const {
      return m(i,col);
    }

    iterator begin() {
      return iterator( m.begin() + col, m.cols() );
    }

    const_iterator begin() const {
      return const_iterator( m.begin() + col, m.cols() );
    }

    iterator end() {
      return begin() + size();
    }

    const_iterator end() const {
      return begin() + size();
    }

  };

  template <class MatrixT>
  struct VectorSize<MatrixCol<MatrixT> > {
    static const int value = MatrixRows<MatrixT>::value;
  };

  /// Extract a column of a matrix as a vector.
  template <class MatrixT>
  inline MatrixCol<MatrixT> select_col( MatrixBase<MatrixT>& matrix, unsigned col ) {
    return MatrixCol<MatrixT>( matrix.impl(), col );
  }

  /// Extract a column of a matrix as a vector (const overload).
  template <class MatrixT>
  inline MatrixCol<const MatrixT> select_col( MatrixBase<MatrixT> const& matrix, unsigned col ) {
    return MatrixCol<const MatrixT>( matrix.impl(), col );
  }


  // *******************************************************************
  // class SubMatrix<MatrixT>
  // An dynamically-sized submatrix (i.e. matrix block) class.
  // *******************************************************************

  /// A submatrix (matrix block) reference object.
  template <class MatrixT>
  class SubMatrix : public MatrixBase<SubMatrix<MatrixT> >
  {
    MatrixT& m_matrix;
    unsigned m_row, m_col;
    unsigned m_rows, m_cols;

  public:
    typedef typename MatrixT::value_type value_type;

    typedef typename boost::mpl::if_<boost::is_const<MatrixT>,
                                     typename MatrixT::const_reference_type,
                                     typename MatrixT::reference_type>::type reference_type;
    typedef typename MatrixT::const_reference_type const_reference_type;

    typedef IndexingMatrixIterator<typename boost::mpl::if_<boost::is_const<MatrixT>, const SubMatrix, SubMatrix>::type> iterator;
    typedef IndexingMatrixIterator<const SubMatrix> const_iterator;

    SubMatrix( MatrixT& m, unsigned row, unsigned col, unsigned rows, unsigned cols ) : 
      m_matrix(m), m_row(row), m_col(col), m_rows(rows), m_cols(cols) {}
    
    SubMatrix& operator=( SubMatrix const& m ) {
      VW_ASSERT( m.rows()==rows() && m.cols()==cols(),
                 ArgumentErr() << "Matrices must have same size in submatrix assignment" );
      std::copy( m.begin(), m.end(), begin() );
      return *this;
    }

    template <class OtherT>
    SubMatrix& operator=( MatrixBase<OtherT> const& m ) {
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(), 
                 ArgumentErr() << "Matrices must have same size in submatrix assignment" );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    unsigned rows() const {
      return m_rows;
    }

    unsigned cols() const {
      return m_cols;
    }

    reference_type operator()( int row, int col ) {
      return m_matrix( row+m_row, col+m_col );
    }

    const_reference_type operator()( int row, int col ) const {
      return m_matrix( row+m_row, col+m_col );
    }

    iterator begin() {
      return iterator( *this, 0, 0 );
    }

    const_iterator begin() const {
      return const_iterator( *this, 0, 0 );
    }

    iterator end() {
      return iterator( *this, rows(), 0 );
    }

    const_iterator end() const {
      return const_iterator( *this, rows(), 0 );
    }

  };

  /// Extract a submatrix, i.e. a matrix block.
  template <class MatrixT>
  inline SubMatrix<MatrixT> submatrix( MatrixBase<MatrixT>& matrix, unsigned row, unsigned col, unsigned rows, unsigned cols ) {
    return SubMatrix<MatrixT>( matrix.impl(), row, col, rows, cols );
  }

  /// Extract a submatrix, i.e. a matrix block (const overlaod).
  template <class MatrixT>
  inline SubMatrix<const MatrixT> submatrix( MatrixBase<MatrixT> const& matrix, unsigned row, unsigned col, unsigned rows, unsigned cols ) {
    return SubMatrix<const MatrixT>( matrix.impl(), row, col, rows, cols );
  }


  // *******************************************************************
  // class MatrixUnaryFunc<MatrixT,FuncT>
  // An unary elementwise matrix function class.
  // *******************************************************************

  template <class MatrixT, class FuncT>
  class MatrixUnaryFunc : public MatrixBase<MatrixUnaryFunc<MatrixT,FuncT> > {
    MatrixT const& m;
    FuncT func;
  public:
    typedef typename boost::result_of<FuncT(typename MatrixT::value_type)>::type value_type;
    
    typedef value_type reference_type;
    typedef value_type const_reference_type;

    MatrixUnaryFunc( MatrixT const& m ) : m(m) {}

    template <class Arg1>
    MatrixUnaryFunc( MatrixT const& m, Arg1 a1 ) : m(m), func(a1) {}

    unsigned rows() const {
      return m.rows();
    }

    unsigned cols() const {
      return m.cols();
    }

    reference_type operator()( int i, int j ) const {
      return func(m(i,j));
    }

    class iterator : public boost::iterator_facade<iterator, value_type, boost::random_access_traversal_tag, value_type> {
      friend class boost::iterator_core_access;

      typename MatrixT::const_iterator i;
      FuncT func;

      bool equal( iterator const& iter ) const { return i==iter.i; }
      ptrdiff_t distance_to( iterator const &iter ) const { return iter.i - i; }
      void increment() { ++i; }
      void decrement() { --i; }
      void advance( ptrdiff_t n ) { i+=n; }
      typename iterator::reference dereference() const { return func(*i); }
    public:
      iterator(typename MatrixT::const_iterator const& i,
               FuncT const& func) : i(i), func(func) {}
    };
    typedef iterator const_iterator;

    iterator begin() const { return iterator(m.begin(),func); }
    iterator end() const { return iterator(m.end(),func); }
  };

  template <class MatrixT, class FuncT>
  struct MatrixRows<MatrixUnaryFunc<MatrixT,FuncT> > {
    static const int value = MatrixRows<MatrixT>::value;
  };

  template <class MatrixT, class FuncT>
  struct MatrixCols<MatrixUnaryFunc<MatrixT,FuncT> > {
    static const int value = MatrixCols<MatrixT>::value;
  };


  // *******************************************************************
  // class MatrixBinaryFunc<Matrix1T,Matrix2T,FuncT>
  // A binary elementwise matrix function class.
  // *******************************************************************

  template <class Matrix1T, class Matrix2T, class FuncT>
  class MatrixBinaryFunc : public MatrixBase<MatrixBinaryFunc<Matrix1T,Matrix2T,FuncT> > {
    Matrix1T const& m1;
    Matrix2T const& m2;
    FuncT func;
  public:
    typedef typename boost::result_of<FuncT(typename Matrix1T::value_type, typename Matrix2T::value_type)>::type value_type;
    
    typedef value_type reference_type;
    typedef value_type const_reference_type;

    MatrixBinaryFunc( Matrix1T const& m1, Matrix2T const& m2 ) : m1(m1), m2(m2) {
      VW_ASSERT( m1.rows() == m2.rows() && m1.cols() == m2.cols(), ArgumentErr() << "Matrices must have same size in MatrixBinaryFunc" );
    }

    template <class Arg1>
    MatrixBinaryFunc( Matrix1T const& m1, Matrix2T const& m2, Arg1 a1 ) : m1(m1), m2(m2), func(a1) {
      VW_ASSERT( m1.rows() == m2.rows() && m1.cols() == m2.cols(), ArgumentErr() << "Matrices must have same size in MatrixBinaryFunc" );
    }

    unsigned rows() const {
      return m1.rows();
    }

    unsigned cols() const {
      return m1.cols();
    }

    reference_type operator()( int i, int j ) const {
      return func(m1(i,j),m2(i,j));
    }

    class iterator : public boost::iterator_facade<iterator, value_type, boost::random_access_traversal_tag, value_type> {
      friend class boost::iterator_core_access;

      typename Matrix1T::const_iterator i1;
      typename Matrix2T::const_iterator i2;
      FuncT func;

      bool equal( iterator const& iter ) const { return (i1==iter.i1) && (i2==iter.i2); }
      ptrdiff_t distance_to( iterator const &iter ) const { return iter.i1 - i1; }
      void increment() { ++i1; ++i2; }
      void decrement() { --i1; --i2; }
      void advance( ptrdiff_t n ) { i1+=n; i2+=n; }
      typename iterator::reference dereference() const { return func(*i1,*i2); }
    public:
      iterator(typename Matrix1T::const_iterator const& i1,
               typename Matrix2T::const_iterator const& i2,
               FuncT const& func) : i1(i1), i2(i2), func(func) {}
    };
    typedef iterator const_iterator;
    
    iterator begin() const { return iterator(m1.begin(),m2.begin(),func); }
    iterator end() const { return iterator(m1.end(),m2.end(),func); }
  };

  template <class Matrix1T, class Matrix2T, class FuncT>
  struct MatrixRows<MatrixBinaryFunc<Matrix1T,Matrix2T,FuncT> > {
    static const int value = (MatrixRows<Matrix1T>::value!=0)?(MatrixRows<Matrix1T>::value):(MatrixRows<Matrix2T>::value);
  };

  template <class Matrix1T, class Matrix2T, class FuncT>
  struct MatrixCols<MatrixBinaryFunc<Matrix1T,Matrix2T,FuncT> > {
    static const int value = (MatrixCols<Matrix1T>::value!=0)?(MatrixCols<Matrix1T>::value):(MatrixCols<Matrix2T>::value);
  };


  // *******************************************************************
  // Matrix iostream interface functions.
  // *******************************************************************

  /// Dumps a matrix to a std::ostream
  template <class MatrixT>
  inline std::ostream& operator<<( std::ostream& os, MatrixBase<MatrixT> const& m ) {
    MatrixT const& mr = m.impl();
    unsigned rows = mr.rows(), cols = mr.cols();
    os << "Matrix" << rows << 'x' << cols << '(';
    for( unsigned r=0; r<rows; ++r ) {
      os << '(' << mr(r,0);
      for( unsigned c=1; c<cols; ++c )
        os << ',' << mr(r,c);
      os << ')';
    }
    return os << ')';
  }


  // *******************************************************************
  // Explicit matrix expression evaluation functions.
  // *******************************************************************

  /// Forces evaluation of an arbitrary matrix expression to a Matrix
  /// object.
  template <class MatrixT>
  Matrix<typename MatrixT::value_type, MatrixRows<MatrixT>::value, MatrixCols<MatrixT>::value>
  inline eval( MatrixBase<MatrixT> const& m ) {
    return m;
  }

  /// Forwarding overload for plain Matrix objects.
  template <class ElemT, int RowsN, int ColsN>
  inline Matrix<ElemT,RowsN,ColsN> const& eval( Matrix<ElemT,RowsN,ColsN> const& m ) {
    return m;
  }


  // *******************************************************************
  // Matrix comparison operators and functions.
  // Note that only equality and inequality operators are provided,
  // in keeping with standard mathematical notation.  Users who want 
  // particular orderings can defined those operators appropriately.
  // *******************************************************************

  /// Equality of two matrices.  Two matrices are considered equal if
  /// they have the same dimensions and their elements are all
  /// equivalent with respect to the standard c++ operator==().
  template <class Matrix1T, class Matrix2T>
  inline bool operator==( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    if (m1.impl().rows() != m2.impl().rows() || 
        m1.impl().cols() != m2.impl().cols()) { return false; }

    typename Matrix1T::const_iterator iter1 = m1.impl().begin();
    typename Matrix2T::const_iterator iter2 = m2.impl().begin();
    for (; iter1 != m1.impl().end(); ++iter1, ++iter2)
      if (*iter1 != *iter2) { return false; }
    return true;
  }

  /// Equality of two matrices measured to within epsilon.  Two
  /// matrices are considered equal if they have the same dimensions
  /// and their elements are equal to within the specified tolerance.
  template <class Matrix1T, class Matrix2T>
  inline bool equal( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2, double epsilon = 0 ) {
    if (m1.impl().rows() != m2.impl().rows() || 
        m1.impl().cols() != m2.impl().cols()) { return false; }

    typename Matrix1T::const_iterator iter1 = m1.impl().begin();
    typename Matrix2T::const_iterator iter2 = m2.impl().begin();
    for (; iter1 != m1.impl().end(); ++iter1, ++iter2)
      if (fabs(*iter1 - *iter2) > epsilon) { return false; }
    return true;
  }

  /// Inequality of two matrices.  Two matrices are considered equal
  /// only if they have the same dimensions and their elements are all
  /// return true when compared with the standard c++ operator==().
  template <class Matrix1T, class Matrix2T>
  inline bool operator!=( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return ! (m1 == m2);
  }

  /// Inequality of two matrices measured to within epsilon.  Two
  /// matrices are considered equal only if they have the same
  /// dimensions and their elements are equal to within the specified
  /// tolerance.
  template <class Matrix1T, class Matrix2T>
  inline bool not_equal( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2, double epsilon = 0 ) {
    return ! equal(m1, m2, epsilon);
  }


  // *******************************************************************
  // Basic elementwise mathematical matrix operators and functions.
  // *******************************************************************

  /// Negation of a matrix.
  template <class MatrixT>
  MatrixUnaryFunc<MatrixT, ArgNegationFunctor >
  inline operator-( MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ArgNegationFunctor>( m.impl() );
  }


  /// Elementwise sum of two matrices.
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgSumFunctor>
  inline elem_sum( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgSumFunctor>( m1.impl(), m2.impl() );
  }

  /// Sum of a matrix and a matrix (same as elem_sum).
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgSumFunctor>
  inline operator+( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return elem_sum( m1, m2 );
  }

  /// Elementwise sum of a scalar and a matrix.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ValArgSumFunctor<ScalarT> > >::type
  inline elem_sum( ScalarT s, MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ValArgSumFunctor<ScalarT> >( m.impl(), s );
  }

  /// Elementwise sum of a matrix and a scalar.
  template <class MatrixT, class ScalarT>
  typename boost::enable_if< IsScalar<ScalarT>, 
                             MatrixUnaryFunc<MatrixT, ArgValSumFunctor<ScalarT> > >::type
  inline elem_sum( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return MatrixUnaryFunc<MatrixT, ArgValSumFunctor<ScalarT> >( m.impl(), s );
  }


  /// Elementwise difference of two matrices.
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgDifferenceFunctor>
  inline elem_diff( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgDifferenceFunctor>( m1.impl(), m2.impl() );
  }

  /// Difference of two matrices (same as elem_diff).
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgDifferenceFunctor>
  inline operator-( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return elem_diff( m1, m2 );
  }

  /// Elementwise difference of a scalar and a matrix.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ValArgDifferenceFunctor<ScalarT> > >::type
  inline elem_diff( ScalarT s, MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ValArgDifferenceFunctor<ScalarT> >( m.impl(), s );
  }

  /// Elementwise difference of a matrix and a scalar.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ArgValDifferenceFunctor<ScalarT> > >::type
  inline elem_diff( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return MatrixUnaryFunc<MatrixT, ArgValDifferenceFunctor<ScalarT> >( m.impl(), s );
  }


  /// Elementwise product of two matrices.
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgProductFunctor>
  inline elem_prod( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgProductFunctor>( m1.impl(), m2.impl() );
  }

  /// Elementwise product of a scalar and a matrix.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ValArgProductFunctor<ScalarT> > >::type
  inline elem_prod( ScalarT s, MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ValArgProductFunctor<ScalarT> >( m.impl(), s );
  }

  /// Product of a scalar and a matrix (same as elem_prod)
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ValArgProductFunctor<ScalarT> > >::type
  inline operator*( ScalarT s, MatrixBase<MatrixT> const& m ) {
    return elem_prod( s, m );
  }

  /// Elementwise product of a matrix and a scalar.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ArgValProductFunctor<ScalarT> > >::type
  inline elem_prod( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return MatrixUnaryFunc<MatrixT, ArgValProductFunctor<ScalarT> >( m.impl(), s );
  }

  /// Product of a matrix and a scalar (same as elem_prod).
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ArgValProductFunctor<ScalarT> > >::type
  inline operator*( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return elem_prod( m, s );
  }


  /// Elementwise quotient of two matrices.
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgQuotientFunctor>
  inline elem_quot( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgQuotientFunctor>( m1.impl(), m2.impl() );
  }

  /// Elementwise quotient of a scalar and a matrix.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ValArgQuotientFunctor<ScalarT> > >::type
  inline elem_quot( ScalarT s, MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ValArgQuotientFunctor<ScalarT> >( m.impl(), s );
  }

  /// Elementwise quotient of a matrix and a scalar.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ArgValQuotientFunctor<ScalarT> > >::type
  inline elem_quot( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return MatrixUnaryFunc<MatrixT, ArgValQuotientFunctor<ScalarT> >( m.impl(), s );
  }

  /// Quotient of a matrix and a scalar (same as elem_quot).
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ArgValQuotientFunctor<ScalarT> > >::type
  inline operator/( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return elem_quot( m, s );
  }


  // *******************************************************************
  // Elementwise vector comparison functions.
  // *******************************************************************

  /// Elementwise equality of two matrices.
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgEqualityFunctor>
  inline elem_eq( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgEqualityFunctor>( m1.impl(), m2.impl() );
  }

  /// Elementwise equality of a scalar and a matrix.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ValArgEqualityFunctor<ScalarT> > >::type
  inline elem_eq( ScalarT s, MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ValArgEqualityFunctor<ScalarT> >( m.impl(), s );
  }

  /// Elementwise equality of a matrix and a scalar.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ArgValEqualityFunctor<ScalarT> > >::type
  inline elem_eq( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return MatrixUnaryFunc<MatrixT, ArgValEqualityFunctor<ScalarT> >( m.impl(), s );
  }


  /// Elementwise inequality of two matrices.
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgInequalityFunctor>
  inline elem_neq( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgInequalityFunctor>( m1.impl(), m2.impl() );
  }

  /// Elementwise inequality of a scalar and a matrix.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ValArgInequalityFunctor<ScalarT> > >::type
  inline elem_neq( ScalarT s, MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ValArgInequalityFunctor<ScalarT> >( m.impl(), s );
  }

  /// Elementwise inequality of a matrix and a scalar.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ArgValInequalityFunctor<ScalarT> > >::type
  inline elem_neq( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return MatrixUnaryFunc<MatrixT, ArgValInequalityFunctor<ScalarT> >( m.impl(), s );
  }


  /// Elementwise less-than of two matrices.
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgLessThanFunctor>
  inline elem_lt( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgLessThanFunctor>( m1.impl(), m2.impl() );
  }

  /// Elementwise less-than of a scalar and a matrix.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ValArgLessThanFunctor<ScalarT> > >::type
  inline elem_lt( ScalarT s, MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ValArgLessThanFunctor<ScalarT> >( m.impl(), s );
  }

  /// Elementwise less-than of a matrix and a scalar.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ArgValLessThanFunctor<ScalarT> > >::type
  inline elem_lt( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return MatrixUnaryFunc<MatrixT, ArgValLessThanFunctor<ScalarT> >( m.impl(), s );
  }


  /// Elementwise less-than-or-equal-to of two matrices.
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgLessThanOrEqualFunctor>
  inline elem_lte( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgLessThanOrEqualFunctor>( m1.impl(), m2.impl() );
  }

  /// Elementwise less-than-or-equal-to of a scalar and a matrix.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ValArgLessThanOrEqualFunctor<ScalarT> > >::type
  inline elem_lte( ScalarT s, MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ValArgLessThanOrEqualFunctor<ScalarT> >( m.impl(), s );
  }

  /// Elementwise less-than-or-equal-to of a matrix and a scalar.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ArgValLessThanOrEqualFunctor<ScalarT> > >::type
  inline elem_lte( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return MatrixUnaryFunc<MatrixT, ArgValLessThanOrEqualFunctor<ScalarT> >( m.impl(), s );
  }


  /// Elementwise greater-than of two matrices.
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgGreaterThanFunctor>
  inline elem_gt( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgGreaterThanFunctor>( m1.impl(), m2.impl() );
  }

  /// Elementwise greater-than of a scalar and a matrix.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ValArgGreaterThanFunctor<ScalarT> > >::type
  inline elem_gt( ScalarT s, MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ValArgGreaterThanFunctor<ScalarT> >( m.impl(), s );
  }

  /// Elementwise greater-than of a matrix and a scalar.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ArgValGreaterThanFunctor<ScalarT> > >::type
  inline elem_gt( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return MatrixUnaryFunc<MatrixT, ArgValGreaterThanFunctor<ScalarT> >( m.impl(), s );
  }


  /// Elementwise greater-than-or-equal-to of two matrices.
  template <class Matrix1T, class Matrix2T>
  MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgGreaterThanOrEqualFunctor>
  inline elem_gte( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return MatrixBinaryFunc<Matrix1T, Matrix2T, ArgArgGreaterThanOrEqualFunctor>( m1.impl(), m2.impl() );
  }

  /// Elementwise greater-than-or-equal-to of a scalar and a matrix.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ValArgGreaterThanOrEqualFunctor<ScalarT> > >::type
  inline elem_gte( ScalarT s, MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ValArgGreaterThanOrEqualFunctor<ScalarT> >( m.impl(), s );
  }

  /// Elementwise greater-than-or-equal-to of a matrix and a scalar.
  template <class ScalarT, class MatrixT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             MatrixUnaryFunc<MatrixT, ArgValGreaterThanOrEqualFunctor<ScalarT> > >::type
  inline elem_gte( MatrixBase<MatrixT> const& m, ScalarT s ) {
    return MatrixUnaryFunc<MatrixT, ArgValGreaterThanOrEqualFunctor<ScalarT> >( m.impl(), s );
  }


  // *******************************************************************
  // Matrix norms and similar functions.
  // *******************************************************************

  /// Matrix 1-norm
  template <class MatrixT>
  inline double norm_1( MatrixBase<MatrixT> const& m ) {
    double result = 0.0;
    unsigned rows = m.impl().rows(), cols = m.impl().cols();
    for( unsigned j=0; j<cols; ++j ) {
      double sum = 0.0;
      for( unsigned i=0; i<rows; ++i )
        sum += fabs( m.impl()(i,j) );
      if( sum > result ) result = sum;
    }
    return result;
  }

  /// Matrix infinity-norm
  template <class MatrixT>
  inline double norm_inf( MatrixBase<MatrixT> const& m ) {
    double result = 0.0;
    unsigned rows = m.impl().rows(), cols = m.impl().cols();
    for( unsigned i=0; i<rows; ++i ) {
      double sum = 0.0;
      for( unsigned j=0; j<cols; ++j )
        sum += fabs( m.impl()(i,j) );
      if( sum > result ) result = sum;
    }
    return result;
  }

  /// Matrix Frobenius norm squared (i.e. sum of squares of the elements)
  template <class MatrixT>
  inline double norm_frobenius_sqr( MatrixBase<MatrixT> const& m ) {
    double result = 0.0;
    unsigned rows = m.impl().rows(), cols = m.impl().cols();
    for( unsigned i=0; i<rows; ++i ) {
      for( unsigned j=0; j<cols; ++j ) {
        typename MatrixT::value_type e = m.impl()(i,j);
        result += e * e;
      }
    }
    return result;
  }

  /// Matrix Frobenius norm
  template <class MatrixT>
  inline double norm_frobenius( MatrixBase<MatrixT> const& m ) {
    return sqrt(norm_frobenius_sqr(m));
  }

  /// Matrix element sum
  template <class MatrixT>
  inline typename MatrixT::value_type sum( MatrixBase<MatrixT> const& m ) {
    typename MatrixT::value_type result = typename MatrixT::value_type();
    unsigned rows = m.impl().rows(), cols = m.impl().cols();
    for( unsigned i=0; i<rows; ++i )
      for( unsigned j=0; j<cols; ++j )
        result += m.impl()(i,j);
    return result;
  }

  /// Matrix element product
  template <class MatrixT>
  inline typename MatrixT::value_type prod( MatrixBase<MatrixT> const& m ) {
    typename MatrixT::value_type result = typename MatrixT::value_type(1);
    unsigned rows = m.impl().rows(), cols = m.impl().cols();
    for( unsigned i=0; i<rows; ++i )
      for( unsigned j=0; j<cols; ++j )
        result *= m.impl()(i,j);
    return result;
  }

  /// Matrix trace
  template <class MatrixT>
  inline typename MatrixT::value_type trace( MatrixBase<MatrixT> const& m ) {
    typename MatrixT::value_type result = typename MatrixT::value_type();
    unsigned mindim = std::min( m.impl().rows(), m.impl().cols() );
    for( unsigned i=0; i<mindim; ++i )
      result += m.impl()(i,i);
    return result;
  }


  // *******************************************************************
  // Matrix vector product.
  // *******************************************************************

  /// A class representing a product of a matrix and a vector.
  template <class MatrixT, class VectorT, bool TransposeN>
  class MatrixVectorProduct : public VectorBase<MatrixVectorProduct<MatrixT,VectorT,TransposeN> > {

    template <class VecT> struct VectorClosure { typedef Vector<typename VecT::value_type, VectorSize<VecT>::value> type; };
    template <class ElemT, int SizeN> struct VectorClosure<Vector<ElemT,SizeN> > { typedef Vector<ElemT,SizeN> const& type; };
    template <class ElemT, int SizeN> struct VectorClosure<const Vector<ElemT,SizeN> > { typedef Vector<ElemT,SizeN> const& type; };
    template <class ElemT, int SizeN> struct VectorClosure<VectorProxy<ElemT,SizeN> > { typedef VectorProxy<ElemT,SizeN> const& type; };
    template <class ElemT, int SizeN> struct VectorClosure<const VectorProxy<ElemT,SizeN> > { typedef VectorProxy<ElemT,SizeN> const& type; };
    typename VectorClosure<VectorT>::type m_vector;

    template <class MatT> struct MatrixClosure { typedef Matrix<typename MatT::value_type> type; };
    template <class ElemT, int RowsN, int ColsN> struct MatrixClosure<Matrix<ElemT,RowsN,ColsN> > { typedef Matrix<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, int RowsN, int ColsN> struct MatrixClosure<const Matrix<ElemT,RowsN,ColsN> > { typedef Matrix<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, int RowsN, int ColsN> struct MatrixClosure<MatrixProxy<ElemT,RowsN,ColsN> > { typedef MatrixProxy<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, int RowsN, int ColsN> struct MatrixClosure<const MatrixProxy<ElemT,RowsN,ColsN> > { typedef MatrixProxy<ElemT,RowsN,ColsN> const& type; };
    typename MatrixClosure<MatrixT>::type m_matrix;

  public:
    typedef typename ProductType<typename MatrixT::value_type, typename VectorT::value_type>::type value_type;
    
    typedef value_type reference_type;
    typedef value_type const_reference_type;

    MatrixVectorProduct( MatrixT const& m, VectorT const& v ) : m_matrix(m), m_vector(v) {}

    unsigned size() const {
      return (TransposeN)?(m_matrix.cols()):(m_matrix.rows());
    }

    reference_type operator()( int i ) const {
      if( TransposeN ) return dot_prod( select_col(m_matrix,i), m_vector );
      else return dot_prod( select_row(m_matrix,i), m_vector );
    }

    class iterator : public boost::iterator_facade<iterator, value_type, boost::random_access_traversal_tag, value_type> {
      friend class boost::iterator_core_access;

      MatrixVectorProduct const& mvp;
      int i;

      bool equal( iterator const& iter ) const { return i==iter.i; }
      ptrdiff_t distance_to( iterator const &iter ) const { return iter.i - i; }
      void increment() { ++i; }
      void decrement() { --i; }
      void advance( ptrdiff_t n ) { i+=n; }
      typename iterator::reference dereference() const { return mvp(i); }
    public:
      iterator(MatrixVectorProduct const& mvp, int i) : mvp(mvp), i(i) {}
    };
    typedef iterator const_iterator;

    iterator begin() const { return iterator(*this,0); }
    iterator end() const { return iterator(*this,size()); }

  };

  /// Product of a matrix and a vector
  template <class MatrixT, class VectorT>
  MatrixVectorProduct<MatrixT,VectorT,false>
  inline operator*( MatrixBase<MatrixT> const& m, VectorBase<VectorT> const& v ) {
    return MatrixVectorProduct<MatrixT,VectorT,false>( m.impl(), v.impl() );
  }

  /// Product of a transposed matrix and a vector
  template <class MatrixT, class VectorT>
  MatrixVectorProduct<MatrixT,VectorT,true>
  inline operator*( MatrixTranspose<MatrixT> const& m, VectorBase<VectorT> const& v ) {
    return MatrixVectorProduct<MatrixT,VectorT,true>( m.child(), v.impl() );
  }

  /// Product of a transposed vector and a matrix
  template <class VectorT, class MatrixT>
  VectorTranspose<const MatrixVectorProduct<MatrixT,VectorT,true> >
  inline operator*( VectorTranspose<VectorT> const& v, MatrixBase<MatrixT> const& m ) {
    return transpose(MatrixVectorProduct<MatrixT,VectorT,true>( m.impl(), v.child() ));
  }


  // *******************************************************************
  // Matrix matrix product.
  // *******************************************************************

  /// A class representing a product of a matrix and a vector.
  template <class Matrix1T, class Matrix2T, bool Transpose1N, bool Transpose2N>
  class MatrixMatrixProduct : public MatrixBase<MatrixMatrixProduct<Matrix1T,Matrix2T,Transpose1N,Transpose2N> > {

    template <class MatT> struct MatrixClosure { typedef Matrix<typename MatT::value_type> type; };
    template <class ElemT, int RowsN, int ColsN> struct MatrixClosure<Matrix<ElemT,RowsN,ColsN> > { typedef Matrix<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, int RowsN, int ColsN> struct MatrixClosure<const Matrix<ElemT,RowsN,ColsN> > { typedef Matrix<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, int RowsN, int ColsN> struct MatrixClosure<MatrixProxy<ElemT,RowsN,ColsN> > { typedef MatrixProxy<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, int RowsN, int ColsN> struct MatrixClosure<const MatrixProxy<ElemT,RowsN,ColsN> > { typedef MatrixProxy<ElemT,RowsN,ColsN> const& type; };
    typename MatrixClosure<Matrix1T>::type m_matrix1;
    typename MatrixClosure<Matrix2T>::type m_matrix2;

  public:
    typedef typename ProductType<typename Matrix1T::value_type, typename Matrix2T::value_type>::type value_type;
    
    typedef value_type reference_type;
    typedef value_type const_reference_type;

    MatrixMatrixProduct( Matrix1T const& m1, Matrix2T const& m2 ) : m_matrix1(m1), m_matrix2(m2) {}

    unsigned rows() const {
      return (Transpose1N)?(m_matrix1.cols()):(m_matrix1.rows());
    }

    unsigned cols() const {
      return (Transpose2N)?(m_matrix2.rows()):(m_matrix2.cols());
    }

    reference_type operator()( int i, int j ) const {
      if     ( (!Transpose1N)&&(!Transpose2N) ) return dot_prod( select_row(m_matrix1,i), select_col(m_matrix2,j) );
      else if( (!Transpose1N)&&( Transpose2N) ) return dot_prod( select_row(m_matrix1,i), select_row(m_matrix2,j) );
      else if( ( Transpose1N)&&(!Transpose2N) ) return dot_prod( select_col(m_matrix1,i), select_col(m_matrix2,j) );
      else                                      return dot_prod( select_col(m_matrix1,i), select_row(m_matrix2,j) );
    }

    typedef IndexingMatrixIterator<const MatrixMatrixProduct> iterator;
    typedef iterator const_iterator;

    iterator begin() const { return iterator(*this,0,0); }
    iterator end() const { return iterator(*this,rows(),0); }
  };

  /// Product of two matrices.
  template <class Matrix1T, class Matrix2T>
  MatrixMatrixProduct<Matrix1T,Matrix2T,false,false>
  inline operator*( MatrixBase<Matrix1T> const& m1, MatrixBase<Matrix2T> const& m2 ) {
    return MatrixMatrixProduct<Matrix1T,Matrix2T,false,false>( m1.impl(), m2.impl() );
  }


  // *******************************************************************
  // Assorted mathematical matrix functions.
  // *******************************************************************

  /// Outer product of two vectors.
  template <class Vector1T, class Vector2T>
  Matrix<typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type,
         (VectorSize<Vector2T>::value==0)?(0):(VectorSize<Vector1T>::value),
         (VectorSize<Vector1T>::value==0)?(0):(VectorSize<Vector2T>::value)>
  inline operator*( VectorBase<Vector1T> const& v1, VectorTranspose<Vector2T> const& v2 ) {
    typedef Matrix<typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type,
      (VectorSize<Vector2T>::value==0)?(0):(VectorSize<Vector1T>::value),
      (VectorSize<Vector1T>::value==0)?(0):(VectorSize<Vector2T>::value)> result_type;
    result_type result;
    result.set_size( v1.impl().size(), v2.size() );
    for( int i=0; i<v1.impl().size(); ++i )
      for( int j=0; j<v2.size(); ++j )
        result(i,j) = v1.impl()(i) * v2(j);
    return result;
  }

  /// Matrix inversion
  template <class MatrixT>
  inline Matrix<typename MatrixT::value_type> inverse( MatrixBase<MatrixT> const& m ) {
    return eval(m.impl()).inverse();
  }

} // namespace math

  // Typedefs for commonly-used static vector types and using 
  // directives for backwards compatability.
  using math::Matrix;
  using math::MatrixBase;
  using math::MatrixProxy;
  typedef Matrix<float64,2,2> Matrix2x2;
  typedef Matrix<float64,3,3> Matrix3x3;
  typedef Matrix<float64,4,4> Matrix4x4;
  typedef Matrix<float32,2,2> Matrix2x2f;
  typedef Matrix<float32,3,3> Matrix3x3f;
  typedef Matrix<float32,4,4> Matrix4x4f;
  typedef Matrix<int32,2,2> Matrix2x2i;
  typedef Matrix<int32,3,3> Matrix3x3i;
  typedef Matrix<int32,4,4> Matrix4x4i;

} // namespace vw

#endif // __VW_MATH_MATRIX_H__
