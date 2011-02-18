// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
///   Determinant via det()
///   Transpose via transpose()
///   Inverse via inverse(), not particularly robust to singularity
///   Matrix norms via norm_1(), norm_inf(), and norm_frobenius()
///   Scalar addition and subtraction via elem_sum() and elem_dif()
///   Elementwise multiplication and division via elem_prod() and elem_quot()
///   One-line creation of identity matrices using identity_matrix()
///
/// We intentionally do *not* define matrix division using inverse(),
/// since users should probably avoid using inverse() anyway unless
/// they know what they are doing.  It's easy enough and preferable
/// for them to write an extra line of code to be explicit about what
/// they want, using e.g. the SVD.
///
#ifndef __VW_MATH_MATRIX_H__
#define __VW_MATH_MATRIX_H__

#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/result_of.hpp>

#include <stack>

#include <vw/Core/Exception.h>
#include <vw/Math/Vector.h>
#include <vw/config.h>

namespace vw {
namespace math {

  /// A type function to compute the number of rows of a matrix
  /// expression at compile time (or zero for dynamically-sized
  /// vectors).
  template <class MatrixT>
  struct MatrixRows {
    const static size_t value = 0;
  };

  /// A type function to compute the number of columns of a matrix
  /// expression at compile time (or zero for dynamically-sized
  /// vectors).
  template <class MatrixT>
  struct MatrixCols {
    const static size_t value = 0;
  };


  template <class MatrixT> class MatrixRow;
  template <class MatrixT> class MatrixCol;

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

    /// Access an individual matrix row, for further access using a second operator[].
    MatrixRow<MatrixT> operator[]( size_t row ) {
      return MatrixRow<MatrixT>( impl(), row );
    }

    /// Access an individual matrix row, for further access using a second operator[].
    MatrixRow<const MatrixT> operator[]( size_t row ) const {
      return MatrixRow<const MatrixT>( impl(), row );
    }

    /// Set the matrix to the identity matrix with the same dimensions.
    void set_identity() {
      VW_ASSERT( impl().rows()==impl().cols(), LogicErr() << "Only square matrices can be identity matrices." );
      size_t n=impl().rows();
      for( size_t i=0; i<n; ++i )
        for( size_t j=0; j<n; ++j )
          impl()(i,j)=(i==j)?(typename MatrixT::value_type(1)):(typename MatrixT::value_type(0));
    }

    /// Set the matrix to the identity matrix with the given dimensions.
    void set_identity( size_t size ) {
      impl().set_size( size, size );
      set_identity();
    }

    /// Set the matrix to entirely zeros.
    void set_zero() {
      for( size_t i=0; i<impl().rows(); ++i )
        for( size_t j=0; j<impl().cols(); ++j )
          impl()(i,j) = (typename MatrixT::value_type(0));
    }

    // Set the entire matrix to ones
    void set_ones() {
      for( size_t i=0; i<impl().rows(); ++i )
        for( size_t j=0; j<impl().cols(); ++j )
          impl()(i,j) = (typename MatrixT::value_type(1));
    }
  };


  // *******************************************************************
  // class MatrixNoTmp<MatrixT>
  // A matrix wrapper class that disables temporaries on assignment.
  // *******************************************************************

  /// A wrapper template class for matrices and matrix expressions.
  /// Provides a mechanism for disabling the use of temporary objects
  /// during matrix assignment in cases where the user deems it safe.
  template <class MatrixT>
  class MatrixNoTmp {
    MatrixT const& m_val;
  public:
    MatrixNoTmp( MatrixT const& val ) : m_val( val ) {}
    MatrixT const& impl() const { return m_val; }
  };

  /// A helper function that provides a mechanism for disabling the use
  /// of temporary objects during matrix assignment in cases where the
  /// user deems it safe.  Use with care.
  template <class MatrixT>
  MatrixNoTmp<MatrixT> no_tmp( MatrixBase<MatrixT> const& val ) {
    return MatrixNoTmp<MatrixT>( val.impl() );
  }


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
  public:
    typedef typename IndexingMatrixIterator::difference_type difference_type;

    IndexingMatrixIterator( MatrixT& matrix, difference_type row, difference_type col ) :
      m_matrix(&matrix), m_row(row), m_col(col) {}

  private:
    friend class boost::iterator_core_access;

    // This has to be a pointer and not a reference because we need to support
    // operator=, and references cannot be reseated.
    MatrixT* m_matrix;
    size_t m_row, m_col;

    bool equal( IndexingMatrixIterator const& iter ) const {
      return m_row==iter.m_row && m_col==iter.m_col;
    }

    difference_type distance_to( IndexingMatrixIterator const &iter ) const {
      difference_type coldiff = (iter.m_col>m_col) ? difference_type(iter.m_col-m_col) : -difference_type(m_col-iter.m_col);
      difference_type rowdiff = (iter.m_row>m_row) ? difference_type(iter.m_row-m_row) : -difference_type(m_row-iter.m_row);
      return coldiff + rowdiff * m_matrix->cols();
    }

    void increment() {
      if( ++m_col == m_matrix->cols() ) {
        m_col=0; ++m_row;
      }
    }

    void decrement() {
      if( m_col==0 ) {
        m_col=m_matrix->cols()-1;
        --m_row;
      }
      else {
        --m_col;
      }
    }

    void advance( difference_type n ) {
      // This safeguards against suprious division by zero troubles encountered
      // on some platforms when performing operations on degenerate matrices.
      if( m_matrix->cols() == 0 ) return;
      if( n < 0 ) {
        difference_type rowdiff = 1 + (-n)/m_matrix->cols();
        m_row -= rowdiff;
        n += rowdiff * m_matrix->cols();
      }
      m_col += n;
      m_row += m_col / m_matrix->cols();
      m_col %= m_matrix->cols();
    }

    typename IndexingMatrixIterator::reference dereference() const {
      return (*m_matrix)(m_row,m_col);
    }
  };


  // *******************************************************************
  // class Matrix<ElemT,RowsN,ColsN>
  // A statically-allocated fixed-dimension matrix class.
  // *******************************************************************

  /// A fixed-dimension mathematical matrix class.
  template <class ElemT, size_t RowsN=0, size_t ColsN=0>
  class Matrix : public MatrixBase<Matrix<ElemT,RowsN,ColsN> >
  {
    typedef boost::array<ElemT,RowsN*ColsN> core_type;
    core_type core_;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef typename core_type::iterator iterator;
    typedef typename core_type::const_iterator const_iterator;

    /// Constructs a matrix of zeroes.
    Matrix() {
      std::memset( core_.c_array(), 0, RowsN*ColsN*sizeof(ElemT) );
    }

    /// Constructs a matrix whose first element is as given.
    Matrix( ElemT e1 ) {
      BOOST_STATIC_ASSERT( RowsN*ColsN >= 1 );
      iterator i=begin();
      *(i++)=e1;
      for( ; i!=end(); ++i ) *i = ElemT();
    }

    /// Constructs a matrix whose first two elements are as given.
    Matrix( ElemT e1, ElemT e2 ) {
      BOOST_STATIC_ASSERT( RowsN*ColsN >= 2 );
      iterator i=begin();
      *(i++)=e1; *(i++)=e2;
      for( ; i!=end(); ++i ) *i = ElemT();
    }

    /// Constructs a matrix whose first three elements are as given.
    Matrix( ElemT e1, ElemT e2, ElemT e3 ) {
      BOOST_STATIC_ASSERT( RowsN*ColsN >= 3 );
      iterator i=begin();
      *(i++)=e1; *(i++)=e2; *(i++)=e3;
      for( ; i!=end(); ++i ) *i = ElemT();
    }

    /// Constructs a matrix whose first four elements are as given.
    Matrix( ElemT e1, ElemT e2, ElemT e3, ElemT e4 ) {
      BOOST_STATIC_ASSERT( RowsN*ColsN >= 4 );
      iterator i=begin();
      *(i++)=e1; *(i++)=e2; *(i++)=e3; *(i++)=e4;
      for( ; i!=end(); ++i ) *i = ElemT();
    }

    /// Constructs a matrix whose first 9 elements are as given.
    Matrix( ElemT e1, ElemT e2, ElemT e3, ElemT e4, ElemT e5,
            ElemT e6, ElemT e7, ElemT e8, ElemT e9 ) {
      BOOST_STATIC_ASSERT( RowsN*ColsN >= 4 );
      iterator i=begin();
      *(i++)=e1; *(i++)=e2; *(i++)=e3; *(i++)=e4;
      *(i++)=e5; *(i++)=e6; *(i++)=e7; *(i++)=e8; *(i++)=e9;
      for( ; i!=end(); ++i ) *i = ElemT();
    }

    /// Constructs a matrix from given densely-packed row-mjor data.
    /// This constructor copies the data.  If you wish to make a
    /// shallow proxy object instead, see vw::MatrixProxy.
    Matrix( const ElemT data[RowsN*ColsN] ) {
      std::copy( data, data+RowsN*ColsN, core_.begin() );
    }

    /// Standard copy constructor.
    Matrix( Matrix const& m ) : core_( m.core_ ) {}

    /// Generalized copy constructor, from arbitrary VW matrix expressions.
    template <class T>
    Matrix( MatrixBase<T> const& m ) {
      VW_ASSERT( m.impl().rows()==RowsN && m.impl().cols()==ColsN, ArgumentErr() << "Matrix must have dimensions " << RowsN << "x" << ColsN << "." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
    }

    /// Standard copy assignment operator.
    Matrix& operator=( Matrix const& m ) {
      Matrix tmp( m );
      core_ = tmp.core_;
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class T>
    Matrix& operator=( MatrixBase<T> const& m ) {
      VW_ASSERT( m.impl().rows()==RowsN && m.impl().cols()==ColsN, ArgumentErr() << "Matrix must have dimensions " << RowsN << "x" << ColsN << "." );
      Matrix tmp( m );
      core_ = tmp.core_;
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW matrix expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class T>
    Matrix& operator=( MatrixNoTmp<T> const& m ) {
      VW_ASSERT( m.impl().rows()==RowsN && m.impl().cols()==ColsN, ArgumentErr() << "Matrix must have dimensions " << RowsN << "x" << ColsN << "." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    /// Returns the number of rows in the matrix.
    size_t rows() const { return RowsN; }

    /// Returns the number of columns in the matrix.
    size_t cols() const { return ColsN; }

    /// Change the size of the matrix.  Elements in memory are preserved when specified.
    void set_size( size_t new_rows, size_t new_cols, bool /*preserve*/ = false ) {
      VW_ASSERT( new_rows==rows() && new_cols==cols(),
                 ArgumentErr() << "Cannot resize a fixed-size Matrix." );
    }

    /// Access an element
    value_type& operator()( size_t row, size_t col ) {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      return core_[row*ColsN+col];
    }

    /// Access an element
    value_type const& operator()( size_t row, size_t col ) const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      return core_[row*ColsN+col];
    }

    value_type *data() {
      return &(operator()(0,0));
    }

    const value_type *data() const {
      return &(operator()(0,0));
    }

    iterator begin() {
      return core_.begin();
    }

    const_iterator begin() const {
      return core_.begin();
    }

    iterator end() {
      return core_.end();
    }

    const_iterator end() const {
      return core_.end();
    }

  };

  template <class ElemT, size_t RowsN, size_t ColsN>
  struct MatrixRows<Matrix<ElemT,RowsN,ColsN> > {
    static const size_t value = RowsN;
  };

  template <class ElemT, size_t RowsN, size_t ColsN>
  struct MatrixCols<Matrix<ElemT,RowsN,ColsN> > {
    static const size_t value = ColsN;
  };


  // *******************************************************************
  // class Matrix<ElemT>
  // A dynamically-allocated arbitrary-dimension matrix class.
  // *******************************************************************

  /// An arbitrary-dimension mathematical matrix class.
  template <class ElemT>
  class Matrix<ElemT,0,0> : public MatrixBase<Matrix<ElemT> > {
    typedef VarArray<ElemT> core_type;
    core_type core_;
    size_t m_rows, m_cols;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef typename core_type::iterator iterator;
    typedef typename core_type::const_iterator const_iterator;

    /// Constructs a matrix with zero size.
    Matrix() : m_rows(0), m_cols(0) {}

    /// Constructs a zero matrix of the given size.
    Matrix( size_t rows, size_t cols )
      : core_(rows*cols), m_rows(rows), m_cols(cols) {}

    /// Constructs a matrix of the given size from given
    /// densely-packed row-mjor data.  This constructor copies the
    /// data.  If you wish to make a shallow proxy object instead,
    /// see vw::MatrixProxy.
    Matrix( size_t rows, size_t cols, const ElemT *data )
      : core_(data,data+rows*cols), m_rows(rows), m_cols(cols) {}

    /// Standard copy constructor.
    // FIXME Do we really need this?
    Matrix( Matrix const& m )
      : core_(m.core_), m_rows(m.m_rows), m_cols(m.m_cols) {}

    /// Generalized copy constructor, from arbitrary VW matrix expressions.
    template <class T>
    Matrix( MatrixBase<T> const& m )
      : core_(m.impl().begin(),m.impl().end()), m_rows(m.impl().rows()), m_cols(m.impl().cols()) {}

    /// Standard copy assignment operator.
    Matrix& operator=( Matrix const& m ) {
      Matrix tmp( m );
      m_rows = tmp.m_rows;
      m_cols = tmp.m_cols;
      core_ = tmp.core_;
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class T>
    Matrix& operator=( MatrixBase<T> const& m ) {
      Matrix tmp( m );
      m_rows = tmp.m_rows;
      m_cols = tmp.m_cols;
      core_ = tmp.core_;
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW matrix expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class T>
    Matrix& operator=( MatrixNoTmp<T> const& m ) {
      if( m.impl().rows()==rows() && m.impl().cols()==cols() ) {
        std::copy( m.impl().begin(), m.impl().end(), begin() );
        return *this;
      }
      else return *this = m.impl();
    }

    /// Returns the number of rows in the matrix.
    size_t rows() const { return m_rows; }

    /// Returns the number of columns in the matrix.
    size_t cols() const { return m_cols; }

    /// Change the size of the matrix.  Elements in memory are preserved when specified.
    void set_size( size_t rows, size_t cols, bool preserve = false ) {
      if( preserve ) {
        VarArray<ElemT> other(rows*cols);
        size_t mr = (std::min)(rows,m_rows);
        size_t mc = (std::min)(cols,m_cols);
        for( size_t r=0; r<mr; ++r )
          for( size_t c=0; c<mc; ++c )
            other[r*cols+c] = core_[r*m_cols+c];
        core_.swap( other );
      }
      else {
        core_.resize(rows*cols,false);
      }
      m_rows = rows;
      m_cols = cols;
    }

    /// Access an element
    value_type& operator()( size_t row, size_t col ) {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      return core_[row*m_cols+col];
    }

    /// Access an element
    value_type const& operator()( size_t row, size_t col ) const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      return core_[row*m_cols+col];
    }

    value_type *data() {
      return &(operator()(0,0));
    }

    const value_type *data() const {
      return &(operator()(0,0));
    }

    iterator begin() {
      return core_.begin();
    }

    const_iterator begin() const {
      return core_.begin();
    }

    iterator end() {
      return core_.end();
    }

    const_iterator end() const {
      return core_.end();
    }

  };


  // *******************************************************************
  // class MatrixProxy<ElemT,RowsN,ColsN>
  // A fixed-dimension matrix proxy class, treating an arbitrary block
  // of memory as a Matrix in packed row-major format.
  // *******************************************************************

  /// A fixed-dimension mathematical matrix class.
  template <class ElemT, size_t RowsN=0, size_t ColsN=0>
  class MatrixProxy : public MatrixBase<MatrixProxy<ElemT,RowsN,ColsN> >
  {
    ElemT *m_ptr;
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
      VW_ASSERT( m.rows()==RowsN && m.cols()==ColsN, ArgumentErr() << "Matrix must have dimensions " << RowsN << "x" << ColsN << "." );
      Matrix<ElemT,RowsN,ColsN> tmp( m );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class T>
    MatrixProxy& operator=( MatrixBase<T> const& m ) {
      VW_ASSERT( m.impl().rows()==RowsN && m.impl().cols()==ColsN, ArgumentErr() << "Matrix must have dimensions " << RowsN << "x" << ColsN << "." );
      Matrix<ElemT,RowsN,ColsN> tmp( m );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW matrix expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class T>
    MatrixProxy& operator=( MatrixNoTmp<T> const& m ) {
      VW_ASSERT( m.impl().rows()==RowsN && m.impl().cols()==ColsN, ArgumentErr() << "Matrix must have dimensions " << RowsN << "x" << ColsN << "." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    /// Returns the number of rows in the matrix.
    size_t rows() const { return RowsN; }

    /// Returns the number of columns in the matrix.
    size_t cols() const { return ColsN; }

    /// Change the size of the matrix.
    /// Elements in memory are preserved when specified.
    void set_size( size_t new_rows, size_t new_cols, bool /*preserve*/ = false ) {
      VW_ASSERT( new_rows==rows() && new_cols==cols(),
                 ArgumentErr() << "Cannot resize matrix proxy." );
    }

    /// Access an element
    value_type& operator()( size_t row, size_t col ) {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      return m_ptr[row*cols()+col];
    }

    /// Access an element
    value_type const& operator()( size_t row, size_t col ) const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
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

  template <class ElemT, size_t RowsN, size_t ColsN>
  struct MatrixRows<MatrixProxy<ElemT,RowsN,ColsN> > {
    static const size_t value = RowsN;
  };

  template <class ElemT, size_t RowsN, size_t ColsN>
  struct MatrixCols<MatrixProxy<ElemT,RowsN,ColsN> > {
    static const size_t value = ColsN;
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
    size_t m_rows, m_cols;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef ElemT* iterator;
    typedef ElemT const* const_iterator;

    /// Constructs a matrix with zero size.
    MatrixProxy( ElemT* ptr, size_t rows, size_t cols )
      : m_ptr(ptr), m_rows(rows), m_cols(cols) {}

    template <class ContainerT>
    MatrixProxy( ContainerT const& container )
      : m_ptr(container.data()), m_rows(container.rows()), m_cols(container.cols()) {}

    /// Standard copy assignment operator.
    MatrixProxy& operator=( MatrixProxy const& m ) {
      VW_ASSERT( m.rows()==rows() && m.cols()==cols(), ArgumentErr() << "Matrix must have dimensions " << rows() << "x" << cols() << "." );
      Matrix<ElemT> tmp( m );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class T>
    MatrixProxy& operator=( MatrixBase<T> const& m ) {
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(), ArgumentErr() << "Matrix must have dimensions " << rows() << "x" << cols() << "." );
      Matrix<ElemT> tmp( m );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW matrix expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class T>
    MatrixProxy& operator=( MatrixNoTmp<T> const& m ) {
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(), ArgumentErr() << "Matrix must have dimensions " << rows() << "x" << cols() << "." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    /// Returns the number of rows in the matrix.
    size_t rows() const { return m_rows; }

    /// Returns the number of columns in the matrix.
    size_t cols() const { return m_cols; }

    /// Change the size of the matrix.
    /// Elements in memory are preserved when specified.
    void set_size( size_t new_rows, size_t new_cols, bool /*preserve*/ = false ) {
      VW_ASSERT( new_rows==rows() && new_cols==cols(),
                 ArgumentErr() << "Cannot resize matrix proxy." );
    }

    /// Access an element
    value_type& operator()( size_t row, size_t col ) {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      return m_ptr[row*m_cols+col];
    }

    /// Access an element
    value_type const& operator()( size_t row, size_t col ) const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
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
  matrix_proxy( DataT* data_ptr, size_t rows, size_t cols) {
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

  public:
    typedef typename MatrixT::value_type value_type;

    typedef typename MatrixT::reference_type reference_type;
    typedef typename MatrixT::const_reference_type const_reference_type;

    typedef IndexingMatrixIterator<MatrixTranspose<MatrixT> > iterator;
    typedef IndexingMatrixIterator<const MatrixTranspose<MatrixT> > const_iterator;

    /// Constructs a matrix transpose.
    MatrixTranspose( MatrixT& matrix ) : m_matrix(matrix) {}

    /// Standard copy assignment operator.
    MatrixTranspose& operator=( MatrixTranspose const& m ) {
      VW_ASSERT( m.rows()==rows() && m.cols()==cols(), ArgumentErr() << "Matrix must have dimensions " << rows() << "x" << cols() << "." );
      Matrix<value_type> tmp( m );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class T>
    MatrixTranspose& operator=( MatrixBase<T> const& m ) {
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(), ArgumentErr() << "Matrix must have dimensions " << rows() << "x" << cols() << "." );
      Matrix<value_type> tmp( m );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW matrix expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class T>
    MatrixTranspose& operator=( MatrixNoTmp<T> const& m ) {
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(), ArgumentErr() << "Matrix must have dimensions " << rows() << "x" << cols() << "." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    /// Returns the underlying non-transposed matrix.
    MatrixT& child() { return m_matrix; }

    /// Returns the underlying non-transposed matrix (const overload).
    MatrixT const& child() const { return m_matrix; }

    /// Returns the number of rows in the matrix.
    size_t rows() const { return m_matrix.cols(); }

    /// Returns the number of columns in the matrix.
    size_t cols() const { return m_matrix.rows(); }

    /// Change the size of the matrix.
    /// Elements in memory are preserved when specified.
    void set_size( size_t new_rows, size_t new_cols, bool /*preserve*/ = false ) {
      VW_ASSERT( new_rows==rows() && new_cols==cols(),
                 ArgumentErr() << "Cannot resize matrix transpose." );
    }

    /// Access an element
    reference_type operator()( size_t row, size_t col ) {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      return child()(col,row);
    }

    /// Access an element
    const_reference_type operator()( size_t row, size_t col ) const {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
      VW_ASSERT( row < rows() && col < cols(), LogicErr() << "operator() ran off end of matrix" );
#endif
      return child()(col,row);
    }

    iterator begin() { return iterator(*this,0,0); }
    const_iterator begin() const { return const_iterator(*this,0,0); }
    iterator end() { return iterator(*this,rows(),0); }
    const_iterator end() const { return const_iterator(*this,rows(),0); }
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
    size_t row;
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

    MatrixRow( MatrixT& m, size_t row ) : m(m), row(row) {}

    /// Standard copy assignment operator.
    MatrixRow& operator=( MatrixRow const& v ) {
      VW_ASSERT( v.size()==size(), ArgumentErr() << "Vectors must have same size in matrix row assignment." );
      Vector<value_type> tmp( v );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class OtherT>
    MatrixRow& operator=( VectorBase<OtherT> const& v ) {
      VW_ASSERT( v.impl().size()==size(), ArgumentErr() << "Vectors must have same size in matrix row assignment." );
      Vector<value_type> tmp( v );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    // For the transposed vctor object
    template <class OtherT>
    MatrixRow& operator=( VectorTranspose<OtherT> const& v ) {
      VW_ASSERT( v.size() == size(), ArgumentErr() << "Vectors must have same size in matrix row assignment.");
      std::copy( v.begin(), v.end(), begin() );
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW vector expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class OtherT>
    MatrixRow& operator=( VectorNoTmp<OtherT> const& v ) {
      VW_ASSERT( v.impl().size()==size(), ArgumentErr() << "Vectors must have same size in matrix row assignment." );
      std::copy( v.impl().begin(), v.impl().end(), begin() );
      return *this;
    }

    MatrixT& child() { return m; }
    MatrixT const& child() const { return m; }

    size_t size() const {
      return child().cols();
    }

    reference_type operator()( size_t i ) { return child()(row,i); }
    const_reference_type operator()( size_t i ) const { return child()(row,i); }
    reference_type operator[]( size_t i ) { return child()(row,i); }
    const_reference_type operator[]( size_t i ) const { return child()(row,i); }

    iterator begin() { return child().begin() + row*child().cols(); }
    const_iterator begin() const { return child().begin() + row*child().cols(); }
    iterator end() { return child().begin() + (row+1)*child().cols(); }
    const_iterator end() const { return child().begin() + (row+1)*child().cols(); }
  };

  template <class MatrixT>
  struct VectorSize<MatrixRow<MatrixT> > {
    static const size_t value = MatrixCols<MatrixT>::value;
  };

  /// Extract a row of a matrix as a vector.
  template <class MatrixT>
  inline MatrixRow<MatrixT> select_row( MatrixBase<MatrixT>& matrix, size_t row ) {
    return MatrixRow<MatrixT>( matrix.impl(), row );
  }

  /// Extract a row of a matrix as a vector (const overload).
  template <class MatrixT>
  inline MatrixRow<const MatrixT> select_row( MatrixBase<MatrixT> const& matrix, size_t row ) {
    return MatrixRow<const MatrixT>( matrix.impl(), row );
  }


  // *******************************************************************
  // class MatrixCol<MatrixT>
  // A matrix column class with vector semantics.
  // *******************************************************************

  /// A matrix column reference object.
  template <class MatrixT>
  class MatrixCol : public VectorBase<MatrixCol<MatrixT> >
  {
    MatrixT & m;
    size_t col;

    template <class IterT>
    class Iterator : public boost::iterator_facade<Iterator<IterT>,
                                                   typename std::iterator_traits<IterT>::value_type,
                                                   boost::random_access_traversal_tag,
                                                   typename std::iterator_traits<IterT>::reference,
                                                   typename std::iterator_traits<IterT>::difference_type>
    {
      public:
        typedef typename Iterator::difference_type difference_type;

        Iterator( IterT const& i, difference_type stride ) : i(i), stride(stride) {}
      private:
        friend class boost::iterator_core_access;

        IterT i;
        difference_type stride;

        bool equal( Iterator const& iter ) const { return i==iter.i; }
        difference_type distance_to( Iterator const &iter ) const { return (iter.i - i) / stride; }
        void increment() { i += stride; }
        void decrement() { i -= stride; }
        void advance( difference_type n ) { i += n*stride; }
        typename Iterator::reference dereference() const { return *i; }
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

    MatrixCol( MatrixT& m, size_t col ) : m(m), col(col) {}

    /// Standard copy assignment operator.
    MatrixCol& operator=( MatrixCol const& v ) {
      VW_ASSERT( v.size()==size(), ArgumentErr() << "Vectors must have same size in matrix column assignment." );
      Vector<value_type> tmp( v );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class OtherT>
    MatrixCol& operator=( VectorBase<OtherT> const& v ) {
      VW_ASSERT( v.impl().size()==size(), ArgumentErr() << "Vectors must have same size in matrix column assignment." );
      Vector<value_type> tmp( v );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW vector expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class OtherT>
    MatrixCol& operator=( VectorNoTmp<OtherT> const& v ) {
      VW_ASSERT( v.impl().size()==size(), ArgumentErr() << "Vectors must have same size in matrix column assignment." );
      std::copy( v.impl().begin(), v.impl().end(), begin() );
      return *this;
    }

    MatrixT& child() { return m; }
    MatrixT const& child() const { return m; }

    size_t size() const {
      return child().rows();
    }

    reference_type operator()( size_t i ) { return child()(i,col); }
    const_reference_type operator()( size_t i ) const { return child()(i,col); }
    reference_type operator[]( size_t i ) { return child()(i,col); }
    const_reference_type operator[]( size_t i ) const { return child()(i,col); }

    iterator begin() { return iterator( child().begin() + col, child().cols() ); }
    const_iterator begin() const {
      return const_iterator( child().begin() + col, child().cols() );
    }
    iterator end() { return begin() + size(); }
    const_iterator end() const { return begin() + size(); }

  };

  template <class MatrixT>
  struct VectorSize<MatrixCol<MatrixT> > {
    static const size_t value = MatrixRows<MatrixT>::value;
  };

  /// Extract a column of a matrix as a vector.
  template <class MatrixT>
  inline MatrixCol<MatrixT> select_col( MatrixBase<MatrixT>& matrix, size_t col ) {
    return MatrixCol<MatrixT>( matrix.impl(), col );
  }

  /// Extract a column of a matrix as a vector (const overload).
  template <class MatrixT>
  inline MatrixCol<const MatrixT> select_col( MatrixBase<MatrixT> const& matrix, size_t col ) {
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
    size_t m_row, m_col;
    size_t m_rows, m_cols;

  public:
    typedef typename MatrixT::value_type value_type;

    typedef typename boost::mpl::if_<boost::is_const<MatrixT>,
                                     typename MatrixT::const_reference_type,
                                     typename MatrixT::reference_type>::type reference_type;
    typedef typename MatrixT::const_reference_type const_reference_type;

    typedef IndexingMatrixIterator<typename boost::mpl::if_<boost::is_const<MatrixT>, const SubMatrix, SubMatrix>::type> iterator;
    typedef IndexingMatrixIterator<const SubMatrix> const_iterator;

    SubMatrix( MatrixT& m, size_t row, size_t col, size_t rows, size_t cols ) :
      m_matrix(m), m_row(row), m_col(col), m_rows(rows), m_cols(cols) {}

    /// Standard copy assignment operator.
    SubMatrix& operator=( SubMatrix const& m ) {
      VW_ASSERT( m.rows()==rows() && m.cols()==cols(), ArgumentErr() << "Matrices must have same size in submatrix assignment." );
      Matrix<value_type> tmp( m );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW matrix expressions.
    template <class T>
    SubMatrix& operator=( MatrixBase<T> const& m ) {
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(), ArgumentErr() << "Matrices must have same size in submatrix assignment." );
      Matrix<value_type> tmp( m );
      std::copy( tmp.begin(), tmp.end(), begin() );
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW matrix expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class T>
    SubMatrix& operator=( MatrixNoTmp<T> const& m ) {
      VW_ASSERT( m.impl().rows()==rows() && m.impl().cols()==cols(), ArgumentErr() << "Matrices must have same size in submatrix assignment." );
      std::copy( m.impl().begin(), m.impl().end(), begin() );
      return *this;
    }

    MatrixT& child() { return m_matrix; }
    MatrixT const& child() const { return m_matrix; }

    size_t rows() const { return m_rows; }
    size_t cols() const { return m_cols; }

    reference_type operator()( size_t row, size_t col ) {
      return child()( row+m_row, col+m_col );
    }
    const_reference_type operator()( size_t row, size_t col ) const {
      return child()( row+m_row, col+m_col );
    }

    iterator begin() { return iterator( *this, 0, 0 ); }
    const_iterator begin() const { return const_iterator( *this, 0, 0 ); }
    iterator end() { return iterator( *this, rows(), 0 ); }
    const_iterator end() const { return const_iterator( *this, rows(), 0 ); }

  };

  /// Extract a submatrix, i.e. a matrix block.
  template <class MatrixT>
  inline SubMatrix<MatrixT> submatrix( MatrixBase<MatrixT>& matrix, size_t row, size_t col, size_t rows, size_t cols ) {
    return SubMatrix<MatrixT>( matrix.impl(), row, col, rows, cols );
  }

  /// Extract a submatrix, i.e. a matrix block (const overlaod).
  template <class MatrixT>
  inline SubMatrix<const MatrixT> submatrix( MatrixBase<MatrixT> const& matrix, size_t row, size_t col, size_t rows, size_t cols ) {
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

    MatrixT const& child() const { return m; }

    size_t rows() const { return m.rows(); }
    size_t cols() const { return m.cols(); }

    reference_type operator()( size_t i, size_t j ) const {
      return func(child()(i,j));
    }

    class iterator : public boost::iterator_facade<iterator, value_type, boost::random_access_traversal_tag, value_type> {
      friend class boost::iterator_core_access;

      typename MatrixT::const_iterator i;
      FuncT func;

      bool equal( iterator const& iter ) const { return i==iter.i; }
      typename iterator::difference_type distance_to( iterator const &iter ) const { return iter.i - i; }
      void increment() { ++i; }
      void decrement() { --i; }
      void advance( typename iterator::difference_type n ) { i+=n; }
      typename iterator::reference dereference() const { return func(*i); }
    public:
      iterator(typename MatrixT::const_iterator const& i,
               FuncT const& func) : i(i), func(func) {}
    };

    typedef iterator const_iterator;
    iterator begin() const { return iterator(child().begin(),func); }
    iterator end() const { return iterator(child().end(),func); }
  };

  template <class MatrixT, class FuncT>
  struct MatrixRows<MatrixUnaryFunc<MatrixT,FuncT> > {
    static const size_t value = MatrixRows<MatrixT>::value;
  };

  template <class MatrixT, class FuncT>
  struct MatrixCols<MatrixUnaryFunc<MatrixT,FuncT> > {
    static const size_t value = MatrixCols<MatrixT>::value;
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

    MatrixBinaryFunc( Matrix1T const& m1, Matrix2T const& m2 ) : m1(m1), m2(m2), func() {
      VW_ASSERT( m1.rows() == m2.rows() && m1.cols() == m2.cols(), ArgumentErr() << "Matrices must have same size in MatrixBinaryFunc" );
    }

    template <class Arg1>
    MatrixBinaryFunc( Matrix1T const& m1, Matrix2T const& m2, Arg1 a1 ) : m1(m1), m2(m2), func(a1) {
      VW_ASSERT( m1.rows() == m2.rows() && m1.cols() == m2.cols(), ArgumentErr() << "Matrices must have same size in MatrixBinaryFunc" );
    }

    Matrix1T const& child1() const { return m1; }
    Matrix2T const& child2() const { return m2; }

    size_t rows() const { return m1.rows(); }
    size_t cols() const { return m1.cols(); }

    reference_type operator()( size_t i, size_t j ) const {
      return func(child1()(i,j),child2()(i,j));
    }

    class iterator : public boost::iterator_facade<iterator, value_type, boost::random_access_traversal_tag, value_type> {
      friend class boost::iterator_core_access;

      typename Matrix1T::const_iterator i1;
      typename Matrix2T::const_iterator i2;
      FuncT func;

      bool equal( iterator const& iter ) const { return (i1==iter.i1) && (i2==iter.i2); }
      typename iterator::difference_type distance_to( iterator const &iter ) const { return iter.i1 - i1; }
      void increment() { ++i1; ++i2; }
      void decrement() { --i1; --i2; }
      void advance( typename iterator::difference_type n ) { i1+=n; i2+=n; }
      typename iterator::reference dereference() const { return func(*i1,*i2); }
    public:
      iterator(typename Matrix1T::const_iterator const& i1,
               typename Matrix2T::const_iterator const& i2,
               FuncT const& func) : i1(i1), i2(i2), func(func) {}
    };

    typedef iterator const_iterator;
    iterator begin() const { return iterator(child1().begin(),child2().begin(),func); }
    iterator end() const { return iterator(child1().end(),child2().end(),func); }
  };

  template <class Matrix1T, class Matrix2T, class FuncT>
  struct MatrixRows<MatrixBinaryFunc<Matrix1T,Matrix2T,FuncT> > {
    static const size_t value = (MatrixRows<Matrix1T>::value!=0)?(MatrixRows<Matrix1T>::value):(MatrixRows<Matrix2T>::value);
  };

  template <class Matrix1T, class Matrix2T, class FuncT>
  struct MatrixCols<MatrixBinaryFunc<Matrix1T,Matrix2T,FuncT> > {
    static const size_t value = (MatrixCols<Matrix1T>::value!=0)?(MatrixCols<Matrix1T>::value):(MatrixCols<Matrix2T>::value);
  };


  // *******************************************************************
  // Matrix iostream interface functions.
  // *******************************************************************

  /// Dumps a matrix to a std::ostream
  template <class MatrixT>
  inline std::ostream& operator<<( std::ostream& os, MatrixBase<MatrixT> const& m ) {
    MatrixT const& mr = m.impl();
    size_t rows = mr.rows(), cols = mr.cols();
    os << "Matrix" << rows << 'x' << cols << '(';
    for( size_t r=0; r<rows; ++r ) {
      os << '(' << mr(r,0);
      for( size_t c=1; c<cols; ++c )
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
  template <class ElemT, size_t RowsN, size_t ColsN>
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
  /// Note that if you have introduced std::equal into the global
  /// namespace then you will have to explicitly request this one as
  /// vw::math::equal().
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

  /// Absolute of a matrix
  template <class MatrixT>
  MatrixUnaryFunc<MatrixT, ArgAbsFunctor >
  inline abs( MatrixBase<MatrixT> const& m ) {
    return MatrixUnaryFunc<MatrixT, ArgAbsFunctor>( m.impl() );
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
    size_t rows = m.impl().rows(), cols = m.impl().cols();
    for( size_t j=0; j<cols; ++j ) {
      double sum = 0.0;
      for( size_t i=0; i<rows; ++i )
        sum += fabs( m.impl()(i,j) );
      if( sum > result ) result = sum;
    }
    return result;
  }

  /// Matrix infinity-norm
  template <class MatrixT>
  inline double norm_inf( MatrixBase<MatrixT> const& m ) {
    double result = 0.0;
    size_t rows = m.impl().rows(), cols = m.impl().cols();
    for( size_t i=0; i<rows; ++i ) {
      double sum = 0.0;
      for( size_t j=0; j<cols; ++j )
        sum += fabs( m.impl()(i,j) );
      if( sum > result ) result = sum;
    }
    return result;
  }

  /// Matrix Frobenius norm squared (i.e. sum of squares of the elements)
  template <class MatrixT>
  inline double norm_frobenius_sqr( MatrixBase<MatrixT> const& m ) {
    double result = 0.0;
    size_t rows = m.impl().rows(), cols = m.impl().cols();
    for( size_t i=0; i<rows; ++i ) {
      for( size_t j=0; j<cols; ++j ) {
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
    size_t rows = m.impl().rows(), cols = m.impl().cols();
    for( size_t i=0; i<rows; ++i )
      for( size_t j=0; j<cols; ++j )
        result += m.impl()(i,j);
    return result;
  }

  /// Matrix element product
  template <class MatrixT>
  inline typename MatrixT::value_type prod( MatrixBase<MatrixT> const& m ) {
    typename MatrixT::value_type result = typename MatrixT::value_type(1);
    size_t rows = m.impl().rows(), cols = m.impl().cols();
    for( size_t i=0; i<rows; ++i )
      for( size_t j=0; j<cols; ++j )
        result *= m.impl()(i,j);
    return result;
  }

  /// Matrix trace
  template <class MatrixT>
  inline typename MatrixT::value_type trace( MatrixBase<MatrixT> const& m ) {
    typename MatrixT::value_type result = typename MatrixT::value_type();
    size_t mindim = std::min( m.impl().rows(), m.impl().cols() );
    for( size_t i=0; i<mindim; ++i )
      result += m.impl()(i,i);
    return result;
  }

  /// Matrix determinant
  template <class MatrixT>
  typename MatrixT::value_type det( MatrixBase<MatrixT> const& m ) {
    typename MatrixT::value_type result = typename MatrixT::value_type();
    VW_ASSERT( m.impl().rows() == m.impl().cols(), vw::ArgumentErr() << "Can only compute determinant of a square matrix." );
    std::stack<std::pair<Matrix<typename MatrixT::value_type>,typename MatrixT::value_type> > s;
    s.push( std::make_pair( m.impl(), 1 ) );
    while( !s.empty() ) {
      Matrix<typename MatrixT::value_type> a = s.top().first;
      typename MatrixT::value_type scale = s.top().second;
      s.pop();
      VW_ASSERT( a.rows() == a.cols(), vw::LogicErr() << "Matrix has become non-square." );
      size_t dim = a.rows();
      Matrix<typename MatrixT::value_type> sub;
      switch( dim ) {
      case 0:
        break;
      case 1:
        result += scale*a(0,0);
        break;
      case 2:
        result += scale*(a(0,0)*a(1,1)-a(0,1)*a(1,0));
        break;
      default:
        {
          sub = submatrix( a, 1, 1, dim-1, dim-1 );
          s.push( std::make_pair( sub, scale*a(0,0) ) );
          scale *= -1;
        }
        for( size_t i=1; i<(dim-1); ++i ) {
          submatrix( sub, 0, 0, dim-1, i ) = submatrix( a, 1, 0, dim-1, i );
          submatrix( sub, 0, i, dim-1, dim-i-1 ) = submatrix( a, 1, i+1, dim-1, dim-i-1 );
          s.push( std::make_pair( sub, scale*a(0,i) ) );
          scale *= -1;
        }
        {
          sub = submatrix( a, 1, 0, dim-1, dim-1 );
          s.push( std::make_pair( sub, scale*a(0,dim-1) ) );
        }
        break;
      }
    }
    return result;
  }

  /// Pull out the maximium element
  template <class MatrixT>
  inline typename MatrixT::value_type max( MatrixBase<MatrixT> const& m ) {
    size_t rows = m.impl().rows(), cols = m.impl().cols();
    if ( rows == 0 ) return typename MatrixT::value_type(0);
    typename MatrixT::value_type result = m.impl()(0,0);
    for ( size_t i=0; i<rows; ++i )
      for ( size_t j=0; j<cols; ++j )
        if ( m.impl()(i,j) > result )
          result = m.impl()(i,j);
    return result;
  }

  /// Pull out the minimium element
  template <class MatrixT>
  inline typename MatrixT::value_type min( MatrixBase<MatrixT> const& m ) {
    size_t rows = m.impl().rows(), cols = m.impl().cols();
    if ( rows == 0 ) return typename MatrixT::value_type(0);
    typename MatrixT::value_type result = m.impl()(0,0);
    for ( size_t i=0; i<rows; ++i )
      for ( size_t j=0; j<cols; ++j )
        if ( m.impl()(i,j) < result )
          result = m.impl()(i,j);
    return result;
  }

  /// Fill the entire matrix with the specified value
  template <class MatrixT, class ValT>
  inline void fill( MatrixBase<MatrixT> &m, ValT const& val ) {
    size_t rows = m.impl().rows(), cols = m.impl().cols();
    for( size_t i=0; i<rows; ++i )
      for( size_t j=0; j<cols; ++j )
        m.impl()(i,j) = typename MatrixT::value_type(val);
  }


  // *******************************************************************
  // Matrix vector product.
  // *******************************************************************

  /// A class representing a product of a matrix and a vector.
  template <class MatrixT, class VectorT, bool TransposeN>
  class MatrixVectorProduct : public VectorBase<MatrixVectorProduct<MatrixT,VectorT,TransposeN> > {

    template <class MatT> struct MatrixClosure { typedef Matrix<typename MatT::value_type> type; };
    template <class ElemT, size_t RowsN, size_t ColsN> struct MatrixClosure<Matrix<ElemT,RowsN,ColsN> > { typedef Matrix<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, size_t RowsN, size_t ColsN> struct MatrixClosure<const Matrix<ElemT,RowsN,ColsN> > { typedef Matrix<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, size_t RowsN, size_t ColsN> struct MatrixClosure<MatrixProxy<ElemT,RowsN,ColsN> > { typedef MatrixProxy<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, size_t RowsN, size_t ColsN> struct MatrixClosure<const MatrixProxy<ElemT,RowsN,ColsN> > { typedef MatrixProxy<ElemT,RowsN,ColsN> const& type; };
    typename MatrixClosure<MatrixT>::type m_matrix;

    template <class VecT> struct VectorClosure { typedef Vector<typename VecT::value_type, VectorSize<VecT>::value> type; };
    template <class ElemT, size_t SizeN> struct VectorClosure<Vector<ElemT,SizeN> > { typedef Vector<ElemT,SizeN> const& type; };
    template <class ElemT, size_t SizeN> struct VectorClosure<const Vector<ElemT,SizeN> > { typedef Vector<ElemT,SizeN> const& type; };
    template <class ElemT, size_t SizeN> struct VectorClosure<VectorProxy<ElemT,SizeN> > { typedef VectorProxy<ElemT,SizeN> const& type; };
    template <class ElemT, size_t SizeN> struct VectorClosure<const VectorProxy<ElemT,SizeN> > { typedef VectorProxy<ElemT,SizeN> const& type; };
    typename VectorClosure<VectorT>::type m_vector;

  public:
    typedef typename ProductType<typename MatrixT::value_type, typename VectorT::value_type>::type value_type;

    typedef value_type reference_type;
    typedef value_type const_reference_type;

    MatrixVectorProduct( MatrixT const& m, VectorT const& v ) : m_matrix(m), m_vector(v) {}

    size_t size() const {
      return (TransposeN)?(m_matrix.cols()):(m_matrix.rows());
    }

    reference_type operator()( size_t i ) const {
      if( TransposeN ) return dot_prod( select_col(m_matrix,i), m_vector );
      else return dot_prod( select_row(m_matrix,i), m_vector );
    }

    class iterator : public boost::iterator_facade<iterator, value_type, boost::random_access_traversal_tag, value_type> {
      friend class boost::iterator_core_access;

      MatrixVectorProduct const& mvp;
      size_t i;

      bool equal( iterator const& iter ) const { return i==iter.i; }
      typename iterator::difference_type distance_to( iterator const &iter ) const { return iter.i - i; }
      void increment() { ++i; }
      void decrement() { --i; }
      void advance( typename iterator::difference_type n ) { i+=n; }
      typename iterator::reference dereference() const { return mvp(i); }
    public:
      iterator(MatrixVectorProduct const& mvp, size_t i) : mvp(mvp), i(i) {}
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
    template <class ElemT, size_t RowsN, size_t ColsN> struct MatrixClosure<Matrix<ElemT,RowsN,ColsN> > { typedef Matrix<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, size_t RowsN, size_t ColsN> struct MatrixClosure<const Matrix<ElemT,RowsN,ColsN> > { typedef Matrix<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, size_t RowsN, size_t ColsN> struct MatrixClosure<MatrixProxy<ElemT,RowsN,ColsN> > { typedef MatrixProxy<ElemT,RowsN,ColsN> const& type; };
    template <class ElemT, size_t RowsN, size_t ColsN> struct MatrixClosure<const MatrixProxy<ElemT,RowsN,ColsN> > { typedef MatrixProxy<ElemT,RowsN,ColsN> const& type; };
    typename MatrixClosure<Matrix1T>::type m_matrix1;
    typename MatrixClosure<Matrix2T>::type m_matrix2;

  public:
    typedef typename ProductType<typename Matrix1T::value_type, typename Matrix2T::value_type>::type value_type;

    typedef value_type reference_type;
    typedef value_type const_reference_type;

    MatrixMatrixProduct( Matrix1T const& m1, Matrix2T const& m2 ) : m_matrix1(m1), m_matrix2(m2) {}

    size_t rows() const {
      return (Transpose1N)?(m_matrix1.cols()):(m_matrix1.rows());
    }

    size_t cols() const {
      return (Transpose2N)?(m_matrix2.rows()):(m_matrix2.cols());
    }

    reference_type operator()( size_t i, size_t j ) const {
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
  // Convenience functions for returning a pre-made identity matrix in
  // one line of code.
  // *******************************************************************

  /// Create a square dynamic identity matrix.
  inline Matrix<double> identity_matrix(size_t size) {
    Matrix<double> id(size,size);
    id.set_identity();
    return id;
  }

  /// Create a square static identity matrix.
  template <size_t DimN>
  inline Matrix<double,DimN,DimN> identity_matrix() {
    Matrix<double,DimN,DimN> id;
    id.set_identity();
    return id;
  }


  // *******************************************************************
  // Convenience functions for returning a pre-made diagonal matrix in
  // one line of code.
  // *******************************************************************

  /// Create a square dynamic identity matrix.
  template <class VectorT>
  inline Matrix<typename VectorT::value_type> diagonal_matrix( VectorT const& diag ) {
    Matrix<typename VectorT::value_type> M(diag.size(),diag.size());
    for( size_t i=0; i<diag.size(); ++i )
      M(i,i) = diag(i);
    return M;
  }


  // *******************************************************************
  // Assorted mathematical matrix functions.
  // *******************************************************************
  template <class Vector1T, class Vector2T>
  Matrix<typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type,
         (VectorSize<Vector2T>::value?(VectorSize<Vector1T>::value):0),
         (VectorSize<Vector1T>::value?(VectorSize<Vector2T>::value):0)>
  inline operator*( VectorBase<Vector1T> const& v1, VectorTranspose<Vector2T> const& v2 ) {
    typedef Matrix<typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type,
                   (VectorSize<Vector2T>::value?(VectorSize<Vector1T>::value):0),
                   (VectorSize<Vector1T>::value?(VectorSize<Vector2T>::value):0)> result_type;
    result_type result;
    result.set_size( v1.impl().size(), v2.size() );
    for( size_t i=0; i<v1.impl().size(); ++i )
      for( size_t j=0; j<v2.size(); ++j )
        result(i,j) = v1.impl()(i) * v2(j);
    return result;
  }



  /// Outer product of two vectors.
  template <class Vector1T, class Vector2T>
  Matrix<typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type,
                              (VectorSize<Vector2T>::value?(VectorSize<Vector1T>::value):0),
                              (VectorSize<Vector1T>::value?(VectorSize<Vector2T>::value):0)>
  inline outer_prod( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    typedef Matrix<typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type,
                   (VectorSize<Vector2T>::value?(VectorSize<Vector1T>::value):0),
                   (VectorSize<Vector1T>::value?(VectorSize<Vector2T>::value):0)> result_type;
    result_type result;
    result.set_size( v1.impl().size(), v2.impl().size() );
    for( size_t i=0; i<v1.impl().size(); ++i )
      for( size_t j=0; j<v2.impl().size(); ++j )
        result(i,j) = v1.impl()(i) * v2.impl()(j);
    return result;
  }

  /// Matrix inversion
  template <class MatrixT>
  inline Matrix<typename MatrixT::value_type> inverse( MatrixBase<MatrixT> const& m ) {

    typedef typename MatrixT::value_type value_type;
    value_type zero = value_type();
    size_t size = m.impl().cols();
    Matrix<value_type> buf = m;

    // Initialize the permutation
    Vector<size_t> pm( size );
    for ( size_t i=0; i<size; ++i ) pm(i) = i;

    // Perform LU decomposition with partial pivoting
    for ( size_t i=0; i<size; ++i) {
      MatrixCol<Matrix<value_type> > mci(buf,i);
      MatrixRow<Matrix<value_type> > mri(buf,i);
      size_t i_norm_inf = i + index_norm_inf( subvector(mci,i,size-i) );
      if ( buf(i_norm_inf,i) == zero )
        vw_throw( MathErr() << "Matrix is singular in inverse()" );
      if ( i_norm_inf != i ) {
        size_t pbuf = pm(i);
        pm(i) = pm(i_norm_inf);
        pm(i_norm_inf) = pbuf;
        Vector<value_type> rowbuf = mri;
        mri = select_row(buf,i_norm_inf);
        select_row(buf,i_norm_inf) = rowbuf;
      }
      if ( i != size-1 ) {
        subvector(mci,i+1,size-i-1) /= buf(i,i);
        submatrix(buf, i+1, i+1, size-i-1, size-i-1) -= outer_prod( subvector(mci, i+1, size-i-1), subvector(mri,i+1,size-i-1) );
      }
    }

    // Build up a permuted identity matrix
    Matrix<value_type> inverse(size,size);
    for ( size_t i=0; i<size; ++i )
      inverse(i,pm(i)) = value_type(1);

    // Divide by the lower-triangular term
    for ( size_t i=0; i<size; ++i ) {
      for ( size_t j=0; j<size; ++j ) {
        value_type t = inverse(i,j);
        if ( t != zero ) {
          for ( size_t k=i+1; k<size; ++k ) {
            inverse(k,j) -= buf(k,i) * t;
          }
        }
      }
    }

    // Divide by the upper-triangular term
    for ( ssize_t i=size-1; i>=0; --i ) {
      for ( ssize_t j=size-1; j>=0; --j ) {
        value_type t = inverse(i,j) /= buf(i,i);
        if ( t != zero ) {
          for ( ssize_t k=i-1; k>=0; --k )
            inverse(k,j) -= buf(k,i) * t;
        }
      }
    }

    return inverse;
  }


} // namespace math

  // Typedefs for commonly-used static vector types and using
  // directives for backwards compatability.
  using math::Matrix;
  using math::MatrixBase;
  using math::MatrixProxy;
  using math::identity_matrix;
  using math::diagonal_matrix;
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
