// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

/// \file Vector.h
///
/// Provides the core Vector<> mathematical vector template classes.  
///
/// Currently I believe we support:
///   Dynamic (Vector<T>) and fixed (Vector<T,N>) vectors
///   Element access via vec(i) and vec[i]
///   Seamless conversion between compatibile vector types
///   Vector<T,N> construction with between zero and four initializer elements
///   Vector addition, subtraction, and negation
///   Scalar multiplication and division
///   Scalar addition and subtraction
///   Dot product via dot_prod()
///   Cross product of 3-element vectors via cross_prod()
///   Sum of elements via sum()
///   Norms via norm_1(), norm_2(), and norm_inf()
///   Elementwise multiplication and division via elem_prod() and elem_quot()
///
#ifndef __VW_MATH__VECTOR_H__
#define __VW_MATH__VECTOR_H__

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/static_assert.hpp>
#include <boost/array.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/utility/result_of.hpp>

#include <vw/Core/Exception.h>
#include <vw/Core/TypeDeduction.h>
#include <vw/Core/Functors.h>

namespace vw {
namespace math {
  namespace bnu = boost::numeric::ublas;

  /// \cond INTERNAL
  // A wrapper around boost::array that provides a constructor taking 
  // a size parameter, which the BOOST uBLAS containers expect to be 
  // available.  With debugging disabled this whole thing should be 
  // optimized out of existence.
  template <class ElemT, int SizeN>
  class FixedArray : public boost::array<ElemT,SizeN> {
  public:
    FixedArray() {}
    FixedArray( unsigned size ) {
      VW_ASSERT( size==SizeN, LogicErr() << "Internal size mismatch in FixedArray." );
    }
  };
  /// \endcond

  // Forward delcaration for friend statements
  class MatrixImplementation;


  template <class VectorT>
  struct VectorBase {

    /// Returns the derived implementation type.
    VectorT& impl() { return *static_cast<VectorT*>(this); }

    /// Returns the derived implementation type.
    VectorT const& impl() const { return *static_cast<VectorT const*>(this); }

    /// Scalar product-assignment operator
    template <class ScalarT>
    typename boost::enable_if< IsScalar<ScalarT>, VectorT& >::type
    operator*=( ScalarT s ) {
      return impl() = impl() * s;
    }

    /// Scalar quotient-assignment operator
    template <class ScalarT>
    typename boost::enable_if< IsScalar<ScalarT>, VectorT& >::type
    operator/=( ScalarT s ) {
      return impl() = impl() / s;
    }
    
    /// Vector sum-assignment operator
    template <class OtherT>
    VectorT& operator+=( VectorBase<OtherT> const& other ) {
      return impl() = impl() + other.impl();
    }

    /// Vector difference-assignment operator
    template <class OtherT>
    VectorT& operator-=( VectorBase<OtherT> const& other ) {
      return impl() = impl() - other.impl();
    }

  };


  /// A fixed-dimension mathematical vector class.
  template <class ElemT, int SizeN = 0>
  class Vector : public VectorBase<Vector<ElemT,SizeN> >
  {
    typedef bnu::vector<ElemT,FixedArray<ElemT,SizeN> > core_type;
    core_type core_;
    friend class MatrixImplementation;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef ElemT* iterator;
    typedef const ElemT* const_iterator;

    /// Constructs a vector of zeroes.
    Vector() : core_(SizeN) {
      core_.clear();
    }

    /// Constructs a vector whose first element is as given.
    Vector( ElemT e1 ) : core_(SizeN) {
      BOOST_STATIC_ASSERT( SizeN >= 1 );
      initialize_( e1 );
      (*this)[0] = e1;
      for( unsigned i=1; i<SizeN; ++i ) (*this)[i] = ElemT();
    }

    /// Constructs a vector whose first two elements are as given.
    Vector( ElemT e1, ElemT e2 ) : core_(SizeN) {
      BOOST_STATIC_ASSERT( SizeN >= 2 );
      (*this)[0] = e1; (*this)[1] = e2;
      for( unsigned i=2; i<SizeN; ++i ) (*this)[i] = ElemT();
    }

    /// Constructs a vector whose first three elements are as given.
    Vector( ElemT e1, ElemT e2, ElemT e3 ) : core_(SizeN) {
      BOOST_STATIC_ASSERT( SizeN >= 3 );
      (*this)[0] = e1; (*this)[1] = e2; (*this)[2] = e3;
      for( unsigned i=3; i<SizeN; ++i ) (*this)[i] = ElemT();
    }

    /// Constructs a vector whose first four elements are as given.
    Vector( ElemT e1, ElemT e2, ElemT e3, ElemT e4 ) : core_(SizeN) {
      BOOST_STATIC_ASSERT( SizeN >= 4 );
      (*this)[0] = e1; (*this)[1] = e2; (*this)[2] = e3; (*this)[3] = e4;
      for( unsigned i=4; i<SizeN; ++i ) (*this)[i] = ElemT();
    }

    /// Constructs a vector from given densely-packed data.  This
    /// constructor copies the data.  If you wish to make a shallow
    /// proxy object instead, see vw::VectorProxy.
    Vector( const ElemT data[SizeN] ) : core_(SizeN) {
      std::copy( data, data+SizeN, begin() );
    }

    /// Standard copy constructor.
    Vector( Vector const& v ) : core_( v.core_ ) {}

    /// Generalized copy constructor, from arbitrary VW vector expressions.
    template <class T>
    Vector( VectorBase<T> const& v ) : core_(SizeN) { 
      VW_ASSERT( v.impl().size()==SizeN, ArgumentErr() << "Vector must have dimension " << SizeN << "." );
      std::copy( v.impl().begin(), v.impl().end(), begin() );
    }

    /// Generalized copy constructor, from arbitrary uBLAS vector expressions.
    template <class T>
    Vector( bnu::vector_expression<T> const& v ) : core_(SizeN) {
      VW_ASSERT( v().size()==SizeN, ArgumentErr() << "Vector must have dimension " << SizeN << "." );
      core_.assign(v);
    }

    /// Standard copy assignment operator.
    Vector& operator=( Vector const& v ) {
      core_ = v.core_;
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW vector expressions.
    template <class T>
    Vector& operator=( VectorBase<T> const& v ) { 
      VW_ASSERT( v.impl().size()==SizeN, ArgumentErr() << "Vector must have dimension " << SizeN << "." );
      std::copy( v.impl().begin(), v.impl().end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary uBLAS vector expressions.
    template <class T>
    Vector& operator=( bnu::vector_expression<T> const& v ) {
      VW_ASSERT( v().size()==SizeN, ArgumentErr() << "Vector must have dimension " << SizeN << "." );
      core_.assign(v);
      return *this;
    }

    /// Returns the size of the vector.
    unsigned size() const { 
      return SizeN;
    }

    /// Change the size of the vector. Elements in memory are preserved when specified.
    void set_size( unsigned new_size, bool preserve = false ) {
      VW_ASSERT( new_size==size(), ArgumentErr() << "Cannot change size of fixed-size Vector." );
    }

    reference_type operator()( unsigned i ) {
      return core_[i];
    }

    const_reference_type operator()( unsigned i ) const {
      return core_[i];
    }

    reference_type operator[]( unsigned i ) {
      return core_[i];
    }

    const_reference_type operator[]( unsigned i ) const {
      return core_[i];
    }

    reference_type x() {
      BOOST_STATIC_ASSERT( SizeN >= 1 );
      return core_[0];
    }

    const_reference_type x() const {
      BOOST_STATIC_ASSERT( SizeN >= 1 );
      return core_[0];
    }

    reference_type y() {
      BOOST_STATIC_ASSERT( SizeN >= 2 );
      return core_[1];
    }

    const_reference_type y() const {
      BOOST_STATIC_ASSERT( SizeN >= 2 );
      return core_[1];
    }

    reference_type z() {
      BOOST_STATIC_ASSERT( SizeN >= 3 );
      return core_[2];
    }

    const_reference_type z() const {
      BOOST_STATIC_ASSERT( SizeN >= 3 );
      return core_[2];
    }

    iterator begin() { 
      return &(core_(0));
    }

    const_iterator begin() const { 
      return &(core_(0));
    }

    iterator end() {
      return &(core_(0)) + SizeN;
    }

    const_iterator end() const {
      return &(core_(0)) + SizeN;
    }

  };


  /// An arbitrary-dimension mathematical vector class.
  template <class ElemT>
  class Vector<ElemT,0> : public VectorBase<Vector<ElemT> > {
    typedef bnu::vector<ElemT> core_type;
    core_type core_;
    friend class MatrixImplementation;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef ElemT* iterator;
    typedef const ElemT* const_iterator;

    /// Constructs a vector with zero size.
    Vector() {}

    /// Constructs a zero vector of the given size.
    Vector( unsigned size ) : core_(size) {
      core_.clear();
    }

    /// Constructs a vector of the given size from given
    /// densely- packed data.  This constructor copies the data.
    /// If you wish to make a shallow proxy object instead, see
    /// vw::VectorProxy.
    Vector( unsigned size, const ElemT *data ) : core_(size) {
      std::copy( data, data+size, begin() );
    }

    /// Standard copy constructor.
    Vector( Vector const& v ) : core_( v.core_ ) {}

    /// Generalized copy constructor, from arbitrary VW vector expressions.
    template <class T>
    Vector( VectorBase<T> const& v ) {
      set_size( v.impl().size() );
      std::copy( v.impl().begin(), v.impl().end(), begin() );
    }

    /// Generalized copy constructor, from arbitrary uBLAS vector expressions.
    template <class T>
    Vector( bnu::vector_expression<T> const& v ) {
      set_size( v().size() );
      core_.assign(v);
    }

    /// Standard copy assignment operator.
    Vector& operator=( Vector const& v ) {
      core_ = v.core_;
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW vector expressions.
    template <class T>
    Vector& operator=( VectorBase<T> const& v ) {
      set_size( v.impl().size() );
      std::copy( v.impl().begin(), v.impl().end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary uBLAS vector expressions.
    template <class T>
    Vector& operator=( bnu::vector_expression<T> const& v ) {
      set_size( v().size() );
      core_.assign(v);
      return *this;
    }

    /// Returns the size of the vector.
    unsigned size() const { 
      return core_.size();
    }

    /// Change the size of the vector. Elements in memory are preserved when specified.
    void set_size( unsigned new_size, bool preserve = false ) {
      core_.resize(new_size, preserve);
    }

    reference_type operator()( unsigned i ) {
      return core_[i];
    }

    const_reference_type operator()( unsigned i ) const {
      return core_[i];
    }

    reference_type operator[]( unsigned i ) {
      return core_[i];
    }

    const_reference_type operator[]( unsigned i ) const {
      return core_[i];
    }

    iterator begin() { 
      return &(core_(0));
    }

    const_iterator begin() const { 
      return &(core_(0));
    }

    iterator end() {
      return &(core_(0)) + size();
    }

    const_iterator end() const {
      return &(core_(0)) + size();
    }

  };


  /// A fixed-dimension mathematical vector class.
  template <class ElemT, int SizeN = 0>
  class VectorProxy : public VectorBase<VectorProxy<ElemT,SizeN> >
  {
    ElemT *m_ptr;
    friend class MatrixImplementation;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef ElemT* iterator;
    typedef const ElemT* const_iterator;

    /// Constructs a vector proxy.
    VectorProxy( ElemT *ptr ) : m_ptr(ptr) {}

    /// Standard copy assignment operator.
    VectorProxy& operator=( VectorProxy const& v ) {
      VW_ASSERT( v.size()==size(), 
                 ArgumentErr() << "Vector must have dimension " 
                 << size() << " in vector proxy assignment." );
      std::copy( v.begin(), v.end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW vector expressions.
    template <class T>
    VectorProxy& operator=( VectorBase<T> const& v ) { 
      VW_ASSERT( v.impl().size()==size(), 
                 ArgumentErr() << "Vector must have dimension " 
                 << size() << " in vector proxy assignment." );
      std::copy( v.impl().begin(), v.impl().end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary uBLAS vector expressions.
    template <class T>
    VectorProxy& operator=( bnu::vector_expression<T> const& v ) {
      VW_ASSERT( v().size()==size(), 
                 ArgumentErr() << "Vector must have dimension " 
                 << size() << " in vector proxy assignment." );
      std::copy( v().begin(), v().end(), begin() );
      return *this;
    }

    /// Returns the size of the vector.
    unsigned size() const { 
      return SizeN;
    }

    /// Change the size of the vector. Elements in memory are preserved when specified.
    void set_size( unsigned new_size, bool preserve = false ) {
      VW_ASSERT( new_size==size(), ArgumentErr() << "Cannot resize a vector proxy." );
    }

    reference_type operator()( unsigned i ) {
      return m_ptr[i];
    }

    const_reference_type operator()( unsigned i ) const {
      return m_ptr[i];
    }

    reference_type operator[]( unsigned i ) {
      return m_ptr[i];
    }

    const_reference_type operator[]( unsigned i ) const {
      return m_ptr[i];
    }

    reference_type x() {
      BOOST_STATIC_ASSERT( SizeN >= 1 );
      return m_ptr[0];
    }

    const_reference_type x() const {
      BOOST_STATIC_ASSERT( SizeN >= 1 );
      return m_ptr[0];
    }

    reference_type y() {
      BOOST_STATIC_ASSERT( SizeN >= 2 );
      return m_ptr[1];
    }

    const_reference_type y() const {
      BOOST_STATIC_ASSERT( SizeN >= 2 );
      return m_ptr[1];
    }

    reference_type z() {
      BOOST_STATIC_ASSERT( SizeN >= 3 );
      return m_ptr[2];
    }

    const_reference_type z() const {
      BOOST_STATIC_ASSERT( SizeN >= 3 );
      return m_ptr[2];
    }

    iterator begin() { 
      return m_ptr;
    }

    const_iterator begin() const { 
      return m_ptr;
    }

    iterator end() {
      return m_ptr + size();
    }

    const_iterator end() const {
      return m_ptr + size();
    }

  };


  /// An arbitrary-dimension vector proxy class.
  template <class ElemT>
  class VectorProxy<ElemT,0> : public VectorBase<VectorProxy<ElemT> > {
    ElemT *m_ptr;
    unsigned m_size;
    friend class MatrixImplementation;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef ElemT* iterator;
    typedef const ElemT* const_iterator;

    /// Constructs a vector with zero size.
    VectorProxy( unsigned size, ElemT *ptr ) : m_ptr(ptr), m_size(size) {}

    VectorProxy& operator=( VectorProxy const& v ) {
      VW_ASSERT( v.size()==size(), 
                 ArgumentErr() << "Vector must have dimension " 
                 << size() << " in vector proxy assignment." );
      std::copy( v.begin(), v.end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW vector expressions.
    template <class T>
    VectorProxy& operator=( VectorBase<T> const& v ) { 
      VW_ASSERT( v.impl().size()==size(), 
                 ArgumentErr() << "Vector must have dimension " 
                 << size() << " in vector proxy assignment." );
      std::copy( v.impl().begin(), v.impl().end(), begin() );
      return *this;
    }

    /// Generalized assignment operator, from arbitrary uBLAS vector expressions.
    template <class T>
    VectorProxy& operator=( bnu::vector_expression<T> const& v ) {
      VW_ASSERT( v().size()==size(), 
                 ArgumentErr() << "Vector must have dimension " 
                 << size() << " in vector proxy assignment." );
      std::copy( v().begin(), v().end(), begin() );
      return *this;
    }

    /// Returns the size of the vector.
    unsigned size() const { 
      return m_size;
    }

    /// Change the size of the vector. Elements in memory are preserved when specified.
    void set_size( unsigned new_size, bool preserve = false ) {
      VW_ASSERT( new_size==size(), ArgumentErr() << "Cannot resize a vector proxy." );
    }

    reference_type operator()( unsigned i ) {
      return m_ptr[i];
    }

    const_reference_type operator()( unsigned i ) const {
      return m_ptr[i];
    }

    reference_type operator[]( unsigned i ) {
      return m_ptr[i];
    }

    const_reference_type operator[]( unsigned i ) const {
      return m_ptr[i];
    }

    iterator begin() { 
      return m_ptr;
    }

    const_iterator begin() const { 
      return m_ptr;
    }

    iterator end() {
      return m_ptr + size();
    }

    const_iterator end() const {
      return m_ptr + size();
    }

  };


  template <class VectorT>
  struct VectorSize {
    const static int value = 0;
  };

  template <class ElemT, int SizeN>
  struct VectorSize<Vector<ElemT,SizeN> > {
    const static int value = SizeN;
  };

  template <class ElemT, int SizeN>
  struct VectorSize<VectorProxy<ElemT,SizeN> > {
    const static int value = SizeN;
  };


  template <class VectorT>
  class SubVector : public VectorBase<SubVector<VectorT> > {
    VectorT& m_vector;
    unsigned m_pos, m_size;
  public:
    typedef typename VectorT::value_type value_type;

    typedef typename boost::mpl::if_<boost::is_const<VectorT>,
                                     typename VectorT::const_reference_type,
                                     typename VectorT::reference_type>::type reference_type;
    typedef typename VectorT::const_reference_type const_reference_type;

    typedef typename boost::mpl::if_<boost::is_const<VectorT>,
                                     typename VectorT::const_iterator,
                                     typename VectorT::iterator>::type iterator;
    typedef typename VectorT::const_iterator const_iterator;

    SubVector( VectorT& v, unsigned pos, unsigned size ) : m_vector(v), m_pos(pos), m_size(size) {}
    
    SubVector& operator=( SubVector const& v ) {
      VW_ASSERT( v.size()==size(), 
                 ArgumentErr() << "Vectors must have same size in subvector assignment" );
      std::copy( v.begin(), v.end(), begin() );
      return *this;
    }

    template <class OtherT>
    SubVector& operator=( VectorBase<OtherT> const& v ) {
      VW_ASSERT( v.impl().size()==m_size, ArgumentErr() << "Vectors must have same size in subvector assignment" );
      std::copy( v.impl().begin(), v.impl().end(), begin() );
      return *this;
    }

    unsigned size() const {
      return m_size;
    }

    reference_type operator()( int i ) {
      return m_vector(m_pos+i);
    }

    const_reference_type operator()( int i ) const {
      return m_vector(m_pos+i);
    }

    iterator begin() {
      return m_vector.begin() + m_pos;
    }

    const_iterator begin() const {
      return m_vector.begin() + m_pos;
    }

    iterator end() {
      return m_vector.begin() + m_pos + m_size;
    }

    const_iterator end() const {
      return m_vector.begin() + m_pos + m_size;
    }

  };

  template <class VectorT>
  inline SubVector<VectorT> subvector( VectorBase<VectorT>& vector, unsigned pos, unsigned size ) {
    return SubVector<VectorT>( vector.impl(), pos, size );
  }

  template <class VectorT>
  inline SubVector<const VectorT> subvector( VectorBase<VectorT> const& vector, unsigned pos, unsigned size ) {
    return SubVector<const VectorT>( vector.impl(), pos, size );
  }



  template <class VectorT>
  class TransposedVector {
    // We want to store Vector objects by reference so that we don't copy 
    // them, but we want to store everything else by value so that we can 
    // return transposed versions of various vector expressions.
    template <class T> struct IsPlainVector : boost::false_type {};
    template <class ElemT, int SizeN> struct IsPlainVector<Vector<ElemT,SizeN> > : boost::true_type{};
    template <class ElemT, int SizeN> struct IsPlainVector<const Vector<ElemT,SizeN> > : boost::true_type{};
    typename boost::mpl::if_<IsPlainVector<VectorT>,VectorT&,VectorT>::type m_vector;
  public:
    typedef typename VectorT::value_type value_type;

    typedef typename boost::mpl::if_<boost::is_const<VectorT>,
                                     typename VectorT::const_reference_type,
                                     typename VectorT::reference_type>::type reference_type;
    typedef typename VectorT::const_reference_type const_reference_type;

    typedef typename boost::mpl::if_<boost::is_const<VectorT>,
                                     typename VectorT::const_iterator,
                                     typename VectorT::iterator>::type iterator;
    typedef typename VectorT::const_iterator const_iterator;

    explicit TransposedVector( VectorT& v ) : m_vector(v) {}
    
    template <class OtherT>
    TransposedVector& operator=( TransposedVector<OtherT> const& v ) {
      VW_ASSERT( v.size()==size(), 
                 ArgumentErr() << "Vectors must have same size in transposed vector assignment" );
      m_vector = v.impl();
      return *this;
    }

    VectorT& impl() {
      return m_vector;
    }

    VectorT const& impl() const {
      return m_vector;
    }

    unsigned size() const {
      return m_vector.size();
    }

    reference_type operator()( int i ) {
      return m_vector(i);
    }

    const_reference_type operator()( int i ) const {
      return m_vector(i);
    }

    iterator begin() {
      return m_vector.begin();
    }

    const_iterator begin() const {
      return m_vector.begin();
    }

    iterator end() {
      return m_vector.begin();
    }

    const_iterator end() const {
      return m_vector.begin();
    }

  };

  template <class VectorT>
  inline TransposedVector<VectorT> transpose( VectorBase<VectorT>& vector ) {
    return TransposedVector<VectorT>( vector.impl() );
  }

  template <class VectorT>
  inline TransposedVector<const VectorT> transpose( VectorBase<VectorT> const& vector ) {
    return TransposedVector<const VectorT>( vector.impl() );
  }

  template <class VectorT>
  inline VectorT& transpose( TransposedVector<VectorT>& vector ) {
    return vector.impl();
  }

  template <class VectorT>
  inline VectorT const& transpose( TransposedVector<VectorT> const& vector ) {
    return vector.impl();
  }


  template <class VectorT, class FuncT>
  class VectorUnaryFunc : public VectorBase<VectorUnaryFunc<VectorT,FuncT> > {
    VectorT const& v;
    FuncT func;
  public:
    typedef typename boost::result_of<FuncT(typename VectorT::value_type)>::type value_type;
    
    typedef value_type reference_type;
    typedef value_type const_reference_type;

    VectorUnaryFunc( VectorT const& v ) : v(v) {}

    template <class Arg1>
    VectorUnaryFunc( VectorT const& v, Arg1 a1 ) : v(v), func(a1) {}

    unsigned size() const {
      return v.size();
    }

    reference_type operator()( int i ) const {
      return func(v(i));
    }

    class iterator : public boost::iterator_facade<iterator, value_type, boost::random_access_traversal_tag, value_type> {
      friend class boost::iterator_core_access;

      typename VectorT::const_iterator i;
      FuncT func;

      bool equal( iterator const& iter ) const { return i==iter.i; }
      ptrdiff_t distance_to( iterator const &iter ) const { return iter.i - i; }
      void increment() { ++i; }
      void decrement() { --i; }
      void advance( ptrdiff_t n ) { i+=n; }
      typename iterator::reference dereference() const { return func(*i); }
    public:
      iterator(typename VectorT::const_iterator const& i,
               FuncT const& func) : i(i), func(func) {}
    };
    typedef iterator const_iterator;

    iterator begin() const { return iterator(v.begin(),func); }
    iterator end() const { return iterator(v.end(),func); }

  };

  template <class VectorT, class FuncT>
  struct VectorSize<VectorUnaryFunc<VectorT,FuncT> > {
    static const int value = VectorSize<VectorT>::value;
  };


  template <class Vector1T, class Vector2T, class FuncT>
  class VectorBinaryFunc : public VectorBase<VectorBinaryFunc<Vector1T,Vector2T,FuncT> > {
    Vector1T const& v1;
    Vector2T const& v2;
    FuncT func;
  public:
    typedef typename boost::result_of<FuncT(typename Vector1T::value_type, typename Vector2T::value_type)>::type value_type;
    
    typedef value_type reference_type;
    typedef value_type const_reference_type;

    VectorBinaryFunc( Vector1T const& v1, Vector2T const& v2 ) : v1(v1), v2(v2) {
      VW_ASSERT( v1.size() == v2.size(), ArgumentErr() << "Vectors must have same size in VectorBinaryFunc" );
    }

    template <class Arg1>
    VectorBinaryFunc( Vector1T const& v1, Vector2T const& v2, Arg1 a1 ) : v1(v1), v2(v2), func(a1) {
      VW_ASSERT( v1.size() == v2.size(), ArgumentErr() << "Vectors must have same size in VectorBinaryFunc" );
    }

    unsigned size() const {
      return v1.size();
    }

    reference_type operator()( int i ) const {
      return func(v1(i),v2(i));
    }

    class iterator : public boost::iterator_facade<iterator, value_type, boost::random_access_traversal_tag, value_type> {
      friend class boost::iterator_core_access;

      typename Vector1T::const_iterator i1;
      typename Vector2T::const_iterator i2;
      FuncT func;

      bool equal( iterator const& iter ) const { return (i1==iter.i1) && (i2==iter.i2); }
      ptrdiff_t distance_to( iterator const &iter ) const { return iter.i1 - i1; }
      void increment() { ++i1; ++i2; }
      void decrement() { --i1; --i2; }
      void advance( ptrdiff_t n ) { i1+=n; i2+=n; }
      typename iterator::reference dereference() const { return func(*i1,*i2); }
    public:
      iterator(typename Vector1T::const_iterator const& i1,
               typename Vector2T::const_iterator const& i2,
               FuncT const& func) : i1(i1), i2(i2), func(func) {}
    };
    typedef iterator const_iterator;

    iterator begin() const { return iterator(v1.begin(),v2.begin(),func); }
    iterator end() const { return iterator(v1.end(),v2.end(),func); }
  };

  template <class Vector1T, class Vector2T, class FuncT>
  struct VectorSize<VectorBinaryFunc<Vector1T,Vector2T,FuncT> > {
    static const int value = (VectorSize<Vector1T>::value!=0)?(VectorSize<Vector1T>::value):(VectorSize<Vector2T>::value);
  };


  /// Dumps a vector to a std::ostream
  template <class VectorT>
  inline std::ostream& operator<<( std::ostream& os, VectorBase<VectorT> const& v ) {
    VectorT const& vr = v.impl();
    unsigned size = vr.size();
    os << '[' << size << "](";
    if( size > 0 ) os << vr(0);
    for( unsigned i=1; i<size; ++i ) os << ',' << vr(i);
    return os << ')';
  }

  /// Dumps a transposed vector to a std::ostream
  template <class VectorT>
  inline std::ostream& operator<<( std::ostream& os, TransposedVector<VectorT> const& v ) {
    VectorT const& vr = v.impl();
    unsigned size = vr.size();
    os << '[' << size << "'](";
    if( size > 0 ) os << vr(0);
    for( unsigned i=1; i<size; ++i ) os << ',' << vr(i);
    return os << ')';
  }


  /// Equality of two vectors.  Two vectors are considered equal if they have the 
  /// same dimensions and their elements are all equivalent with respect to the 
  /// standard c++ operator==(). 
  template <class Vector1T, class Vector2T>
  inline bool operator==( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    if (v1.impl().size() != v2.impl().size()) { return false; }

    typename Vector1T::const_iterator iter1 = v1.impl().begin();
    typename Vector2T::const_iterator iter2 = v2.impl().begin();
    for (; iter1 != v1.impl().end(); ++iter1, ++iter2)
      if (*iter1 != *iter2) { return false; }
    return true;
  }

  /// Equality of two vectors measured to within epsilon.  Two vectors
  /// are considered equal if they have the same dimensions and their
  /// elements are equal to within the specified tolerance. 
  template <class Vector1T, class Vector2T>
  inline bool equal( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2, double epsilon = 0 ) {
    if (v1.impl().size() != v2.impl().size()) { return false; }

    typename Vector1T::const_iterator iter1 = v1.impl().begin();
    typename Vector2T::const_iterator iter2 = v2.impl().begin();
    for (; iter1 != v1.impl().end(); ++iter1, ++iter2)
      if (fabs(*iter1 - *iter2) >= epsilon) { return false; }
    return true;
  }

  /// Inequality of two vectors.  Two vectors are considered equal
  /// only if they have the same dimensions and their elements are all
  /// return true when compared with the standard c++ operator==().
  template <class Vector1T, class Vector2T>
  inline bool operator!=( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return ! (v1 == v2);
  }

  /// Inequality of two vectors measured to within epsilon.  Two
  /// vectors are considered equal only if they have the same
  /// dimensions and their elements are equal to within the specified
  /// tolerance.
  template <class Vector1T, class Vector2T>
  inline bool not_equal( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2, double epsilon = 0 ) {
    return ! equal(v1, v2, epsilon);
  }


  /// Negation of a vector
  template <class VectorT>
  VectorUnaryFunc<VectorT, ArgNegationFunctor>
  inline operator-( VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ArgNegationFunctor>( v.impl() );
  }

  /// Negation of a transposed vector.
  template <class VectorT>
  TransposedVector<const VectorUnaryFunc<VectorT, ArgNegationFunctor> >
  inline operator-( TransposedVector<VectorT> const& v ) {
    return transpose(-v.impl());
  }


  /// Elementwise sum of two vectors.
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgSumFunctor>
  inline elem_sum( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgSumFunctor>( v1.impl(), v2.impl() );
  }

  /// Sum of two vectors (same as elem_sum).
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgSumFunctor>
  inline operator+( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return elem_sum( v1, v2 );
  }

  /// Sum of two transposed vectors.
  template <class Vector1T, class Vector2T>
  TransposedVector<const VectorBinaryFunc<Vector1T, Vector2T, ArgArgSumFunctor> >
  inline operator+( TransposedVector<Vector1T> const& v1, TransposedVector<Vector2T> const& v2 ) {
    return transpose(v1.impl()+v2.impl());
  }

  /// Elementwise sum of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValSumFunctor<ScalarT> > >::type
  inline elem_sum( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ArgValSumFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise sum of a vector and a scalar.
  template <class VectorT, class ScalarT>
  typename boost::enable_if< IsScalar<ScalarT>, 
                             VectorUnaryFunc<VectorT, ArgValSumFunctor<ScalarT> > >::type
  inline elem_sum( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValSumFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise difference of two vectors.
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgDifferenceFunctor>
  inline elem_diff( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgDifferenceFunctor>( v1.impl(), v2.impl() );
  }

  /// Difference of two vectors (same as elem_diff).
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgDifferenceFunctor>
  inline operator-( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return elem_diff( v1, v2 );
  }

  /// Difference of two transposed vectors.
  template <class Vector1T, class Vector2T>
  TransposedVector<const VectorBinaryFunc<Vector1T, Vector2T, ArgArgDifferenceFunctor> >
  inline operator-( TransposedVector<Vector1T> const& v1, TransposedVector<Vector2T> const& v2 ) {
    return transpose(v1.impl()-v2.impl());
  }

  /// Elementwise difference of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgDifferenceFunctor<ScalarT> > >::type
  inline elem_diff( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgDifferenceFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise difference of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValDifferenceFunctor<ScalarT> > >::type
  inline elem_diff( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValDifferenceFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise product of two vectors.
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgProductFunctor>
  inline elem_prod( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgProductFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise product of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgProductFunctor<ScalarT> > >::type
  inline elem_prod( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgProductFunctor<ScalarT> >( v.impl(), s );
  }

  /// Product of a scalar and a vector (same as elem_prod).
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgProductFunctor<ScalarT> > >::type
  inline operator*( ScalarT s, VectorBase<VectorT> const& v ) {
    return elem_prod( s, v );
  }

  /// Product of a scalar and a transposed vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             TransposedVector<const VectorUnaryFunc<VectorT, ValArgProductFunctor<ScalarT> > > >::type
  inline operator*( ScalarT s, TransposedVector<VectorT> const& v ) {
    return transpose(s*v.impl());
  }

  /// Elementwise product of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValProductFunctor<ScalarT> > >::type
  inline elem_prod( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValProductFunctor<ScalarT> >( v.impl(), s );
  }

  /// Product of a vector and a scalar (same as elem_prod).
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValProductFunctor<ScalarT> > >::type
  inline operator*( VectorBase<VectorT> const& v, ScalarT s ) {
    return elem_prod( v, s );
  }

  /// Product of a transposed vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             TransposedVector<const VectorUnaryFunc<VectorT, ArgValProductFunctor<ScalarT> > > >::type
  inline operator*( TransposedVector<VectorT> const& v, ScalarT s ) {
    return transpose(v.impl()*s);
  }


  /// Elementwise quotient of two vectors.
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgQuotientFunctor>
  inline elem_quot( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgQuotientFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise quotient of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgQuotientFunctor<ScalarT> > >::type
  inline elem_quot( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgQuotientFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise quotient of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValQuotientFunctor<ScalarT> > >::type
  inline elem_quot( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValQuotientFunctor<ScalarT> >( v.impl(), s );
  }

  /// Quotient of a vector and a scalar (same as elem_quot).
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValQuotientFunctor<ScalarT> > >::type
  inline operator/( VectorBase<VectorT> const& v, ScalarT s ) {
    return elem_quot( v, s );
  }

  /// Quotient of a transposed vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             TransposedVector<const VectorUnaryFunc<VectorT, ArgValQuotientFunctor<ScalarT> > > >::type
  inline operator/( TransposedVector<VectorT> const& v, ScalarT s ) {
    return transpose(v.impl()/s);
  }


  /// Elementwise equality operator
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgEqualityFunctor>
  inline elem_eq( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgEqualityFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise equality operator of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgEqualityFunctor<ScalarT> > >::type
  inline elem_eq( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgEqualityFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise equality operator of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValEqualityFunctor<ScalarT> > >::type
  inline elem_eq( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValEqualityFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise inequality operator
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgInequalityFunctor>
  inline elem_neq( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgInequalityFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise inequality operator of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgInequalityFunctor<ScalarT> > >::type
  inline elem_neq( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgInequalityFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise inequality operator of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValInequalityFunctor<ScalarT> > >::type
  inline elem_neq( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValInequalityFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise less than operator
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgLessThanFunctor>
  inline elem_lt( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgLessThanFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise less than operator of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgLessThanFunctor<ScalarT> > >::type
  inline elem_lt( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgLessThanFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise less than operator of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValLessThanFunctor<ScalarT> > >::type
  inline elem_lt( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValLessThanFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise less than or equal to operator
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgLessThanOrEqualFunctor>
  inline elem_lte( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgLessThanOrEqualFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise less than or equal to operator of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgLessThanOrEqualFunctor<ScalarT> > >::type
  inline elem_lte( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgLessThanOrEqualFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise less than or equal to operator of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValLessThanOrEqualFunctor<ScalarT> > >::type
  inline elem_lte( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValLessThanOrEqualFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise greater than operator
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgGreaterThanFunctor>
  inline elem_gt( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgGreaterThanFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise greater than operator of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgGreaterThanFunctor<ScalarT> > >::type
  inline elem_gt( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgGreaterThanFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise greater than operator of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValGreaterThanFunctor<ScalarT> > >::type
  inline elem_gt( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValGreaterThanFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise greater than or equal to operator
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgGreaterThanOrEqualFunctor>
  inline elem_gte( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgGreaterThanOrEqualFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise less than or equal to operator of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgGreaterThanOrEqualFunctor<ScalarT> > >::type
  inline elem_gte( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgGreaterThanOrEqualFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise less than or equal to operator of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValGreaterThanOrEqualFunctor<ScalarT> > >::type
  inline elem_gte( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValGreaterThanOrEqualFunctor<ScalarT> >( v.impl(), s );
  }


  /// Vector 1-norm
  template <class VectorT>
  inline typename VectorT::value_type norm_1( VectorBase<VectorT> const& v ) {
    double result = 0.0;
    typename VectorT::const_iterator i = v.impl().begin(), end = v.impl().end();
    for( ; i != end ; ++i ) result += fabs( *i );
    return result;
  }

  /// Square of vector 2-norm
  template <class VectorT>
  inline double norm_2_sqr( VectorBase<VectorT> const& v ) {
    double result = 0.0;
    typename VectorT::const_iterator i = v.impl().begin(), end = v.impl().end();
    for( ; i != end ; ++i ) result += (*i) * (*i);
    return result;
  }

  /// Vector 2-norm
  template <class VectorT>
  inline double norm_2( VectorBase<VectorT> const& v ) {
    return sqrt( norm_2_sqr(v) );
  }

  /// Vector infinity-norm
  template <class VectorT>
  inline typename VectorT::value_type norm_inf( VectorBase<VectorT> const& v ) {
    double result = 0.0;
    typename VectorT::const_iterator i = v.impl().begin(), end = v.impl().end();
    for( ; i != end ; ++i ) {
      double a = fabs( *i );
      if( a > result ) result = a;
    }
    return result;
  }

  /// Sum of elements
  template <class VectorT>
  inline typename VectorT::value_type sum( VectorBase<VectorT> const& v ) {
    typename VectorT::const_iterator i = v.impl().begin(), end = v.impl().end();
    if( i == end ) return typename VectorT::value_type();
    typename VectorT::value_type result = *(i++);
    for( ; i != end ; ++i ) result += *i;
    return result;
  }

  /// Product of elements
  template <class VectorT>
  inline typename VectorT::value_type prod( VectorBase<VectorT> const& v ) {
    typename VectorT::const_iterator i = v.impl().begin(), end = v.impl().end();
    if( i == end ) return typename VectorT::value_type(1);
    typename VectorT::value_type result = *(i++);
    for( ; i != end ; ++i ) result *= *i;
    return result;
  }

  /// Returns a normalized (unit) vector.
  template <class VectorT>
  VectorUnaryFunc<VectorT, ArgValQuotientFunctor<double> >
  inline normalize( VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ArgValQuotientFunctor<double> >( v.impl(), norm_2(v) );
  }

  /// Vector dot product
  template <class Vector1T, class Vector2T>
  typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type
  inline dot_prod( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type result = 
      typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type();
    typename Vector1T::const_iterator i1 = v1.impl().begin(), end1 = v1.impl().end();
    typename Vector2T::const_iterator i2 = v2.impl().begin();
    for( ; i1 != end1 ; ++i1, ++i2 ) result += *i1 * *i2;
    return result;
  }

  /// Vector inner product via transpose
  template <class Vector1T, class Vector2T>
  typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type
  inline operator*( TransposedVector<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return dot_prod( v1.impl(), v2 );
  }

  /// Vector cross product. (Only valid for 3-element vectors.)
  template <class Vector1T, class Vector2T>
  Vector<typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type, 3>
  inline cross_prod( VectorBase<Vector1T> const& v1_, VectorBase<Vector2T> const& v2_ ) {
    typedef Vector<typename ProductType<typename Vector1T::value_type, typename Vector2T::value_type>::type, 3> result_type;
    Vector1T const& v1 = v1_.impl();
    Vector2T const& v2 = v2_.impl();
    return result_type( v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0] );
  }

} // namespace math

  using math::Vector;
  using math::VectorBase;
  using math::VectorProxy;
  typedef Vector<double,2> Vector2;
  typedef Vector<double,3> Vector3;
  typedef Vector<double,4> Vector4;

} // namespace vw

#endif // __VW_MATH__VECTOR_H__
