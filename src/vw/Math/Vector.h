// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
///   Dynamic and fixed VectorProxy proxy objects
///   Subvector access via subvector()
///   Vector transposition via transpose()
///   Explicit expression evaluation via eval()
///   Printing of vectors to ostreams
///   Equality of vectors with or without an epsilon
///   Vector negation, addition, subtraction, abs
///   Scalar multiplication and division
///   Elementwise vector addition, subtraction, multiplication, and division
///   Elementwise scalar addition, subtraction, multiplication, and division
///   Vector negation, addition, subtraction, abs of transposed vectors
///   Scalar multiplication and division of transposed vectors
///   Elementwise comprison operations.
///   Norms via norm_1(), norm_2(), norm_2_sqr(), and norm_inf()
///   Sum and product of elements via sum() and prod()
///   Vector normalization via normalize()
///   Vector homogenization and dehomogenization via hom() and dehom()
///   Dot/inner product via dot_prod(v1,v2) or transpose(v1)*v2
///   Cross product of 3-element vectors via cross_prod()
///   Real and imaginary parts of complex vectors via real() and imag()
///
#ifndef __VW_MATH_VECTOR_H__
#define __VW_MATH_VECTOR_H__

#include <cstring> // for memset
#include <vector>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/static_assert.hpp>
#include <boost/array.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/utility/result_of.hpp>

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/VarArray.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Functors.h>
#include <vw/Math/Functors.h>

namespace vw {
namespace math {

  /// A type function to compute the dimension of a vector expression
  /// at compile time (or zero for dynamically-sized vectors).
  template <class VectorT>
  struct VectorSize {
    const static size_t value = 0;
  };


  // *******************************************************************
  // class VectorBase<VectorT>
  // A CRTP vector base class.
  // *******************************************************************

  /// A CRTP base class for vectors and vector expressions.
  /// Provides a mechanism for restricting function arguments to
  /// vectors, and provides the various arithmetic assignment
  /// operators.
  template <class VectorT>
  struct VectorBase {

    /// Returns the derived implementation type.
    VectorT& impl() { return *static_cast<VectorT*>(this); }

    /// Returns the derived implementation type.
    VectorT const& impl() const { return *static_cast<VectorT const*>(this); }

    /// Sum-assignment operator
    template <class T>
    VectorT& operator+=( T const& v ) {
      return impl() = impl() + v;
    }

    /// Difference-assignment operator
    template <class T>
    VectorT& operator-=( T const& v ) {
      return impl() = impl() - v;
    }

    /// Product-assignment operator
    template <class T>
    VectorT& operator*=( T s ) {
      return impl() = impl() * s;
    }

    /// Quotient-assignment operator
    template <class T>
    VectorT& operator/=( T s ) {
      return impl() = impl() / s;
    }

  };


  // *******************************************************************
  // class VectorNoTmp<VectorT>
  // A vector wrapper class that disables temporaries on assignment.
  // *******************************************************************

  /// A wrapper template class for vectors and vector expressions.
  /// Provides a mechanism for disabling the use of temporary objects
  /// during vector assignment in cases where the user deems it safe.
  template <class VectorT>
  class VectorNoTmp {
    VectorT const& m_val;
  public:
    VectorNoTmp( VectorT const& val ) : m_val( val ) {}
    VectorT const& impl() const { return m_val; }
  };

  /// A helper function that provides a mechanism for disabling the use
  /// of temporary objects during vector assignment in cases where the
  /// user deems it safe.  Use with care.
  template <class VectorT>
  VectorNoTmp<VectorT> no_tmp( VectorBase<VectorT> const& val ) {
    return VectorNoTmp<VectorT>( val.impl() );
  }


  /// This helper class allows overriding the basic vector assignment
  /// operations in specific cases for efficiency, using template
  /// specialization.
  template <class DstVecT, class SrcVecT>
  struct VectorAssignImpl {
    static void assign( DstVecT& dst, SrcVecT const& src ) {
      std::copy( src.begin(), src.end(), dst.begin() );
    }
  };

  /// This helper class allows overriding the basic vector clearing
  /// operation in specific cases for efficiency, using template
  /// specialization.
  template <class VectorT>
  struct VectorClearImpl {
    static void clear( VectorT& vec ) {
      std::fill( vec.begin(), vec.end(), typename VectorT::value_type() );
    }
  };

  // *******************************************************************
  // class IndexingVectorIterator<VectorT>
  // A general purpose vector iterator type.
  // *******************************************************************

  template <class VectorT>
  class IndexingVectorIterator : public boost::iterator_facade<IndexingVectorIterator<VectorT>,
    typename boost::mpl::if_<boost::is_const<VectorT>,
    const typename VectorT::value_type,
    typename VectorT::value_type>::type,
    boost::random_access_traversal_tag,
    typename boost::mpl::if_<boost::is_const<VectorT>,
    typename VectorT::const_reference_type,
    typename VectorT::reference_type>::type> {

  public:
    typedef typename IndexingVectorIterator::difference_type difference_type;

    IndexingVectorIterator( VectorT& vector, difference_type index ) :
      m_vector(&vector), m_index(index) {}
  private:
    friend class boost::iterator_core_access;

    // This has to be a pointer and not a reference because we need to support
    // operator=, and references cannot be reseated.
    VectorT* m_vector;
    difference_type m_index;

    bool equal( IndexingVectorIterator const& iter ) const {
      return m_index==iter.m_index;
    }

    difference_type distance_to( IndexingVectorIterator const& iter ) const {
      return difference_type(iter.m_index)-difference_type(m_index);
    }

    void increment() { ++m_index; }
    void decrement() { --m_index; }

    void advance( difference_type n ) {
      if ( m_vector->size() == 0 ) return;
      m_index += n;
    }

    typename IndexingVectorIterator::reference dereference() const {
      return (*m_vector)[m_index];
    }
  };

  // *******************************************************************
  // class Vector<ElemT,SizeN>
  // A statically-allocated fixed-dimension vector class.
  // *******************************************************************

  /// A fixed-dimension mathematical vector class.
  template <class ElemT, size_t SizeN = 0>
  class Vector : public VectorBase<Vector<ElemT,SizeN> >
  {
    typedef boost::array<ElemT,SizeN> core_type;
    core_type core_;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef typename core_type::iterator iterator;
    typedef typename core_type::const_iterator const_iterator;

    /// Constructs a vector of zeroes.
    Vector() {
      VectorClearImpl<Vector>::clear(*this);
    }

    /// Constructs a vector whose first element is as given.
    Vector( ElemT e1 ) {
      BOOST_STATIC_ASSERT( SizeN >= 1 );
      (*this)[0] = e1;
      for( size_t i=1; i<SizeN; ++i ) (*this)[i] = ElemT();
    }

    /// Constructs a vector whose first two elements are as given.
    Vector( ElemT e1, ElemT e2 ) {
      BOOST_STATIC_ASSERT( SizeN >= 2 );
      core_[0] = e1; core_[1] = e2;
      for( size_t i=2; i<SizeN; ++i ) (*this)[i] = ElemT();
    }

    /// Constructs a vector whose first three elements are as given.
    Vector( ElemT e1, ElemT e2, ElemT e3 ) {
      BOOST_STATIC_ASSERT( SizeN >= 3 );
      (*this)[0] = e1; (*this)[1] = e2; (*this)[2] = e3;
      for( size_t i=3; i<SizeN; ++i ) (*this)[i] = ElemT();
    }

    /// Constructs a vector whose first four elements are as given.
    Vector( ElemT e1, ElemT e2, ElemT e3, ElemT e4 ) {
      BOOST_STATIC_ASSERT( SizeN >= 4 );
      (*this)[0] = e1; (*this)[1] = e2; (*this)[2] = e3; (*this)[3] = e4;
      for( size_t i=4; i<SizeN; ++i ) (*this)[i] = ElemT();
    }

    /// Constructs a vector whose first five elements are as given.
    Vector( ElemT e1, ElemT e2, ElemT e3, ElemT e4, ElemT e5 ) {
      BOOST_STATIC_ASSERT( SizeN >= 5 );
      (*this)[0] = e1; (*this)[1] = e2; (*this)[2] = e3; (*this)[3] = e4; (*this)[4] = e5;
      for( size_t i=5; i<SizeN; ++i ) (*this)[i] = ElemT();
    }

    /// Constructs a vector from given densely-packed data.  This
    /// constructor copies the data.  If you wish to make a shallow
    /// proxy object instead, see vw::VectorProxy.
    Vector( const ElemT data[SizeN] ) {
      std::copy( data, data+SizeN, begin() );
    }

    /// Standard copy constructor.
    Vector( Vector const& v ) : core_( v.core_ ) {}

    /// Generalized copy constructor, from arbitrary VW vector expressions.
    template <class T>
    Vector( VectorBase<T> const& v ) {
      VW_ASSERT( v.impl().size()==SizeN, ArgumentErr() << "Vector must have dimension " << SizeN << "." );
      VectorAssignImpl<Vector,T>::assign(*this,v.impl());
    }

    /// Standard copy assignment operator.
    Vector& operator=( Vector const& v ) {
      Vector tmp( v );
      core_ = tmp.core_;
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW vector expressions.
    template <class T>
    Vector& operator=( VectorBase<T> const& v ) {
      VW_ASSERT( v.impl().size()==SizeN, ArgumentErr() << "Vector must have dimension " << SizeN << "." );
      Vector tmp( v );
      core_ = tmp.core_;
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW vector expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class T>
    Vector& operator=( VectorNoTmp<T> const& v ) {
      VW_ASSERT( v.impl().size()==SizeN, ArgumentErr() << "Vector must have dimension " << SizeN << "." );
      VectorAssignImpl<Vector,T>::assign(*this,v.impl());
      return *this;
    }

    /// Returns the size of the vector.
    size_t size() const {
      return SizeN;
    }

    /// Change the size of the vector. Elements in memory are preserved when specified.
    void set_size( size_t new_size, bool /*preserve*/ = false ) {
      VW_ASSERT( new_size==size(), ArgumentErr() << "Cannot change size of fixed-size Vector." );
    }

    reference_type operator()( size_t i ) {
      return core_[i];
    }

    const_reference_type operator()( size_t i ) const {
      return core_[i];
    }

    reference_type operator[]( size_t i ) {
      return core_[i];
    }

    const_reference_type operator[]( size_t i ) const {
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

  template <class ElemT, size_t SizeN>
  struct VectorSize<Vector<ElemT,SizeN> > {
    const static size_t value = SizeN;
  };

  template <class DstElemT, class SrcVecT>
  struct VectorAssignImpl<Vector<DstElemT,2>,SrcVecT> {
    static void assign( Vector<DstElemT,2>& dst, SrcVecT const& src ) {
      dst(0) = (DstElemT) src(0);
      dst(1) = (DstElemT) src(1);
    }
  };

  template <class DstElemT, class SrcVecT>
  struct VectorAssignImpl<Vector<DstElemT,3>,SrcVecT> {
    static void assign( Vector<DstElemT,3>& dst, SrcVecT const& src ) {
      dst(0) = (DstElemT) src(0);
      dst(1) = (DstElemT) src(1);
      dst(2) = (DstElemT) src(2);
    }
  };

  template <class DstElemT, class SrcVecT>
  struct VectorAssignImpl<Vector<DstElemT,4>,SrcVecT> {
    static void assign( Vector<DstElemT,4>& dst, SrcVecT const& src ) {
      dst(0) = (DstElemT) src(0);
      dst(1) = (DstElemT) src(1);
      dst(2) = (DstElemT) src(2);
      dst(3) = (DstElemT) src(3);
    }
  };

  template <class ElemT>
  struct VectorClearImpl<Vector<ElemT,2> > {
    static void clear( Vector<ElemT,2>& v ) {
      v(0) = v(1) = ElemT(0);
    }
  };

  template <class ElemT>
  struct VectorClearImpl<Vector<ElemT,3> > {
    static void clear( Vector<ElemT,3>& v ) {
      v(0) = v(1) = v(2) = ElemT();
    }
  };

  template <class ElemT>
  struct VectorClearImpl<Vector<ElemT,4> > {
    static void clear( Vector<ElemT,4>& v ) {
      v(0) = v(1) = v(2) = v(3) = ElemT();
    }
  };

  template <class ElemT, size_t N>
  struct VectorClearImpl<Vector<ElemT,N> > {
    static void clear( Vector<ElemT,N>& v ) {
      std::memset( &v(0), 0, N*sizeof(ElemT) );
    }
  };


  // *******************************************************************
  // class Vector<ElemT>
  // A dynamically-allocated arbitrary-dimension vector class.
  // *******************************************************************

  /// An arbitrary-dimension mathematical vector class.
  template <class ElemT>
  class Vector<ElemT,0> : public VectorBase<Vector<ElemT> > {
    typedef VarArray<ElemT> core_type;
    core_type core_;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef typename core_type::iterator iterator;
    typedef typename core_type::const_iterator const_iterator;

    /// Constructs a vector with zero size.
    Vector() {}

    /// Constructs a zero vector of the given size.
    Vector( size_t size ) : core_(size) {}

    /// Constructs a vector of the given size from given
    /// densely- packed data.  This constructor copies the data.
    /// If you wish to make a shallow proxy object instead, see
    /// vw::VectorProxy.
    Vector( size_t size, const ElemT *data ) : core_(data, data+size) {}

    /// Standard copy constructor.
    Vector( Vector const& v ) : core_( v.core_ ) {}

    /// Generalized copy constructor, from arbitrary VW vector expressions.
    template <class T>
    Vector( VectorBase<T> const& v ) {
      set_size( v.impl().size() );
      VectorAssignImpl<Vector,T>::assign(*this,v.impl());
    }

    /// Standard copy assignment operator.
    Vector& operator=( Vector const& v ) {
      Vector tmp( v );
      core_ = tmp.core_;
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW vector expressions.
    template <class T>
    Vector& operator=( VectorBase<T> const& v ) {
      Vector tmp( v );
      core_ = tmp.core_;
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW vector expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class T>
    Vector& operator=( VectorNoTmp<T> const& v ) {
      if( v.impl().size()==size() ) {
        VectorAssignImpl<Vector,T>::assign(*this,v.impl());
        return *this;
      }
      else return *this = v.impl();
    }

    /// Returns the size of the vector.
    size_t size() const {
      return core_.size();
    }

    /// Change the size of the vector. Elements in memory are preserved when specified.
    void set_size( size_t new_size, bool preserve = false ) {
      core_.resize(new_size, preserve);
    }

    reference_type operator()( size_t i ) {
      return core_[i];
    }

    const_reference_type operator()( size_t i ) const {
      return core_[i];
    }

    reference_type operator[]( size_t i ) {
      return core_[i];
    }

    const_reference_type operator[]( size_t i ) const {
      return core_[i];
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
  // class VectorProxy<ElemT,SizeN>
  // A fixed-dimension vector proxy class, treating an arbitrary block
  // of memory as a Vector.
  // *******************************************************************

  /// A fixed-dimension mathematical vector class.
  template <class ElemT, size_t SizeN = 0>
  class VectorProxy : public VectorBase<VectorProxy<ElemT,SizeN> >
  {
    ElemT *m_ptr;
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
      VW_ASSERT( v.size()==size(), ArgumentErr() << "Vector must have dimension " << size() << " in vector proxy assignment." );
      Vector<value_type> tmp( v );
      VectorAssignImpl<VectorProxy,Vector<value_type> >::assign(*this,tmp);
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW vector expressions.
    template <class T>
    VectorProxy& operator=( VectorBase<T> const& v ) {
      VW_ASSERT( v.impl().size()==size(), ArgumentErr() << "Vector must have dimension " << size() << " in vector proxy assignment." );
      Vector<value_type> tmp( v );
      VectorAssignImpl<VectorProxy,Vector<value_type> >::assign(*this,tmp);
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW vector expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class T>
    VectorProxy& operator=( VectorNoTmp<T> const& v ) {
      VW_ASSERT( v.impl().size()==size(), ArgumentErr() << "Vector must have dimension " << size() << " in vector proxy assignment." );
      VectorAssignImpl<VectorProxy,T>::assign(*this,v.impl());
      return *this;
    }

    /// Returns the size of the vector.
    size_t size() const {
      return SizeN;
    }

    /// Change the size of the vector. Elements in memory are preserved when specified.
    void set_size( size_t new_size, bool /*preserve*/ = false ) {
      VW_ASSERT( new_size==size(), ArgumentErr() << "Cannot resize a vector proxy." );
    }

    reference_type operator()( size_t i ) {
      return m_ptr[i];
    }

    const_reference_type operator()( size_t i ) const {
      return m_ptr[i];
    }

    reference_type operator[]( size_t i ) {
      return m_ptr[i];
    }

    const_reference_type operator[]( size_t i ) const {
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

  template <class ElemT, size_t SizeN>
  struct VectorSize<VectorProxy<ElemT,SizeN> > {
    const static size_t value = SizeN;
  };


  // *******************************************************************
  // class VectorProxy<ElemT>
  // An arbitrary-dimension vector proxy class, treating an arbitrary
  // block of memory as a Vector.
  // *******************************************************************

  /// An arbitrary-dimension vector proxy class.
  template <class ElemT>
  class VectorProxy<ElemT,0> : public VectorBase<VectorProxy<ElemT> > {
    ElemT *m_ptr;
    size_t m_size;
  public:
    typedef ElemT value_type;

    typedef ElemT& reference_type;
    typedef ElemT const& const_reference_type;

    typedef ElemT* iterator;
    typedef const ElemT* const_iterator;

    /// Constructs a vector with zero size.
    VectorProxy( size_t size, ElemT *ptr ) : m_ptr(ptr), m_size(size) {}

    /// Standard copy assignment operator.
    VectorProxy& operator=( VectorProxy const& v ) {
      VW_ASSERT( v.size()==size(), ArgumentErr() << "Vector must have dimension " << size() << " in vector proxy assignment." );
      Vector<value_type> tmp( v );
      VectorAssignImpl<VectorProxy,Vector<value_type> >::assign(*this,tmp);
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW vector expressions.
    template <class T>
    VectorProxy& operator=( VectorBase<T> const& v ) {
      VW_ASSERT( v.impl().size()==size(), ArgumentErr() << "Vector must have dimension " << size() << " in vector proxy assignment." );
      Vector<value_type> tmp( v );
      VectorAssignImpl<VectorProxy,Vector<value_type> >::assign(*this,tmp);
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW vector expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class T>
    VectorProxy& operator=( VectorNoTmp<T> const& v ) {
      VW_ASSERT( v.impl().size()==size(), ArgumentErr() << "Vector must have dimension " << size() << " in vector proxy assignment." );
      VectorAssignImpl<VectorProxy,T>::assign(*this,v.impl());
      return *this;
    }

    /// Returns the size of the vector.
    size_t size() const {
      return m_size;
    }

    /// Change the size of the vector. Elements in memory are preserved when specified.
    void set_size( size_t new_size, bool /*preserve*/ = false ) {
      VW_ASSERT( new_size==size(), ArgumentErr() << "Cannot resize a vector proxy." );
    }

    reference_type operator()( size_t i ) {
      return m_ptr[i];
    }

    const_reference_type operator()( size_t i ) const {
      return m_ptr[i];
    }

    reference_type operator[]( size_t i ) {
      return m_ptr[i];
    }

    const_reference_type operator[]( size_t i ) const {
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

  /// Shallow proxy view of an block of memory as a vector.
  template <class DataT>
  VectorProxy<DataT>
  vector_proxy( DataT* data, size_t size) {
    return VectorProxy<DataT>( data, size );
  }


  // *******************************************************************
  // class VectorTranspose<VectorT>
  // A transposed vector wrapper class.
  // *******************************************************************

  /// A transposed vector class.  This class represents the transposed
  /// version of a vector.  Note that unlike a transposed matrix, a
  /// transposed vector is not simply another vector.  For this reason
  /// VectorTransposed does not derive from VectorBase<>, since that
  /// would permit a variety of mathematically meaninless expressions
  /// to be created with unexpected results.  An unfortunate side
  /// effect of this is that explicit overloads must be provided for
  /// all operations on transposed vectors.  For the moment only a
  /// limited set of operations are supported, including primarily the
  /// mathematical operators.
  template <class VectorT>
  class VectorTranspose {
    // We want to store Vector objects by reference so that we don't copy
    // them, but we want to store everything else by value so that we can
    // return transposed versions of various vector expressions.
    template <class T> struct VectorClosure { typedef T type; };
    template <class ElemT, size_t SizeN> struct VectorClosure<Vector<ElemT,SizeN> > { typedef Vector<ElemT,SizeN>& type; };
    template <class ElemT, size_t SizeN> struct VectorClosure<const Vector<ElemT,SizeN> > { typedef Vector<ElemT,SizeN> const& type; };
    typename VectorClosure<VectorT>::type m_vector;
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

    explicit VectorTranspose( VectorT& v ) : m_vector(v) {}

    template <class OtherT>
    VectorTranspose& operator=( VectorTranspose<OtherT> const& v ) {
      VW_ASSERT( v.size()==size(), ArgumentErr() << "Vectors must have same size in transposed vector assignment" );
      m_vector = v.child();
      return *this;
    }

    VectorT& child() {
      return m_vector;
    }

    VectorT const& child() const {
      return m_vector;
    }

    size_t size() const {
      return m_vector.size();
    }

    reference_type operator()( size_t i ) {
      return m_vector(i);
    }

    const_reference_type operator()( size_t i ) const {
      return m_vector(i);
    }

    reference_type operator[]( size_t i ) {
      return m_vector[i];
    }

    const_reference_type operator[]( size_t i ) const {
      return m_vector[i];
    }

    iterator begin() {
      return m_vector.begin();
    }

    const_iterator begin() const {
      return m_vector.begin();
    }

    iterator end() {
      return m_vector.end();
    }

    const_iterator end() const {
      return m_vector.end();
    }

  };

  /// Vector transpose.
  template <class VectorT>
  inline VectorTranspose<VectorT> transpose( VectorBase<VectorT>& vector ) {
    return VectorTranspose<VectorT>( vector.impl() );
  }

  /// Vector transpose (const overload).
  template <class VectorT>
  inline VectorTranspose<const VectorT> transpose( VectorBase<VectorT> const& vector ) {
    return VectorTranspose<const VectorT>( vector.impl() );
  }

  /// Vector transpose (transpose overload).
  template <class VectorT>
  inline VectorT& transpose( VectorTranspose<VectorT>& vector ) {
    return vector.child();
  }

  /// Vector transpose (const transpose overload).
  template <class VectorT>
  inline VectorT const& transpose( VectorTranspose<VectorT> const& vector ) {
    return vector.child();
  }


  // *******************************************************************
  // class SubVector<VectorT>
  // An dynamically-sized subvector class.
  // *******************************************************************

  /// A subvector class.
  template <class VectorT>
  class SubVector : public VectorBase<SubVector<VectorT> > {
    VectorT& m_vector;
    size_t m_pos, m_size;
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

    SubVector( VectorT& v, size_t pos, size_t size ) : m_vector(v), m_pos(pos), m_size(size) {}

    /// Standard copy assignment operator.
    SubVector& operator=( SubVector const& v ) {
      VW_ASSERT( v.size()==size(), ArgumentErr() << "Vectors must have same size in subvector assignment" );
      VectorAssignImpl<SubVector,Vector<value_type> >::assign(*this,v.impl());
      return *this;
    }

    /// Generalized assignment operator, from arbitrary VW vector expressions.
    template <class OtherT>
    SubVector& operator=( VectorBase<OtherT> const& v ) {
      VW_ASSERT( v.impl().size()==m_size, ArgumentErr() << "Vectors must have same size in subvector assignment" );
      VectorAssignImpl<SubVector,OtherT >::assign(*this,v.impl());
      return *this;
    }

    /// Temporary-free generalized assignment operator, from arbitrary VW vector expressions.
    /// This is a performance-optimizing function to be used with caution!
    template <class OtherT>
    SubVector& operator=( VectorNoTmp<OtherT> const& v ) {
      VW_ASSERT( v.impl().size()==m_size, ArgumentErr() << "Vectors must have same size in subvector assignment" );
      VectorAssignImpl<SubVector,OtherT>::assign(*this,v.impl());
      return *this;
    }

    size_t size() const {
      return m_size;
    }

    reference_type operator()( size_t i ) {
      return m_vector(m_pos+i);
    }

    const_reference_type operator()( size_t i ) const {
      return m_vector(m_pos+i);
    }

    reference_type operator[]( size_t i ) {
      return m_vector[m_pos+i];
    }

    const_reference_type operator[]( size_t i ) const {
      return m_vector[m_pos+i];
    }

    iterator begin() {
      return m_vector.begin() + m_pos;
    }

    const_iterator begin() const {
      return const_cast<VectorT const&>(m_vector).begin() + m_pos;
    }

    iterator end() {
      return m_vector.begin() + m_pos + m_size;
    }

    const_iterator end() const {
      return const_cast<VectorT const&>(m_vector).begin() + m_pos + m_size;
    }

  };

  template <class VectorT>
  inline SubVector<VectorT> subvector( VectorBase<VectorT>& vector, size_t pos, size_t size ) {
    return SubVector<VectorT>( vector.impl(), pos, size );
  }

  template <class VectorT>
  inline SubVector<const VectorT> subvector( VectorBase<VectorT> const& vector, size_t pos, size_t size ) {
    return SubVector<const VectorT>( vector.impl(), pos, size );
  }


  // *******************************************************************
  // class VectorUnaryFunc<VectorT,FuncT>
  // An unary elementwise vector function class.
  // *******************************************************************

  /// An unary elementwise vector function class.
  template <class VectorT, class FuncT>
  class VectorUnaryFunc : public VectorBase<VectorUnaryFunc<VectorT,FuncT> >
  {
    VectorT const& v;
    FuncT func;
  public:
    typedef typename boost::result_of<FuncT(typename VectorT::value_type)>::type value_type;

    typedef value_type reference_type;
    typedef value_type const_reference_type;

    class iterator : public boost::iterator_facade<iterator,
                                                   value_type,
                                                   boost::random_access_traversal_tag,
                                                   reference_type>
    {
      friend class boost::iterator_core_access;

      typename VectorT::const_iterator i;
      FuncT func;

      bool equal( iterator const& iter ) const { return i==iter.i; }
      typename iterator::difference_type distance_to( iterator const &iter ) const { return iter.i - i; }
      void increment() { ++i; }
      void decrement() { --i; }
      void advance( typename iterator::difference_type n ) { i+=n; }
      typename iterator::reference dereference() const { return func(*i); }
    public:
      iterator(typename VectorT::const_iterator const& i,
               FuncT const& func) : i(i), func(func) {}
    };

    typedef iterator const_iterator;

    VectorUnaryFunc( VectorT const& v ) : v(v) {}

    template <class Arg1>
    VectorUnaryFunc( VectorT const& v, Arg1 a1 ) : v(v), func(a1) {}

    size_t size() const {
      return v.size();
    }

    reference_type operator()( size_t i ) const {
      return func(v(i));
    }

    reference_type operator[]( size_t i ) const {
      return func(v[i]);
    }

    iterator begin() const { return iterator(v.begin(),func); }
    iterator end() const { return iterator(v.end(),func); }
  };

  template <class VectorT, class FuncT>
  struct VectorSize<VectorUnaryFunc<VectorT,FuncT> > {
    static const size_t value = VectorSize<VectorT>::value;
  };


  // *******************************************************************
  // class VectorBinaryFunc<Vector1T,Vector2T,FuncT>
  // A binary elementwise vector function class.
  // *******************************************************************

  /// A binary elementwise vector function class.
  template <class Vector1T, class Vector2T, class FuncT>
  class VectorBinaryFunc : public VectorBase<VectorBinaryFunc<Vector1T,Vector2T,FuncT> > {
    Vector1T const& v1;
    Vector2T const& v2;
    FuncT func;
  public:
    typedef typename boost::result_of<FuncT(typename Vector1T::value_type, typename Vector2T::value_type)>::type value_type;

    typedef value_type reference_type;
    typedef value_type const_reference_type;

    class iterator : public boost::iterator_facade<iterator,
                                                   value_type,
                                                   boost::random_access_traversal_tag,
                                                   reference_type>
    {
      friend class boost::iterator_core_access;

      typename Vector1T::const_iterator i1;
      typename Vector2T::const_iterator i2;
      FuncT func;

      bool equal( iterator const& iter ) const { return (i1==iter.i1) && (i2==iter.i2); }
      typename iterator::difference_type distance_to( iterator const &iter ) const { return iter.i1 - i1; }
      void increment() { ++i1; ++i2; }
      void decrement() { --i1; --i2; }
      void advance( typename iterator::difference_type n ) { i1+=n; i2+=n; }
      typename iterator::reference dereference() const { return func(*i1,*i2); }
    public:
      iterator(typename Vector1T::const_iterator const& i1,
               typename Vector2T::const_iterator const& i2,
               FuncT const& func) : i1(i1), i2(i2), func(func) {}
    };

    typedef iterator const_iterator;

    VectorBinaryFunc( Vector1T const& v1, Vector2T const& v2 ) : v1(v1), v2(v2), func() {
      VW_ASSERT( v1.size() == v2.size(), ArgumentErr() << "Vectors must have same size in VectorBinaryFunc" );
    }

    template <class Arg1>
    VectorBinaryFunc( Vector1T const& v1, Vector2T const& v2, Arg1 a1 ) : v1(v1), v2(v2), func(a1) {
      VW_ASSERT( v1.size() == v2.size(), ArgumentErr() << "Vectors must have same size in VectorBinaryFunc" );
    }

    size_t size() const {
      return v1.size();
    }

    reference_type operator()( size_t i ) const {
      return func(v1(i),v2(i));
    }

    reference_type operator[]( size_t i ) const {
      return func(v1[i],v2[i]);
    }

    iterator begin() const { return iterator(v1.begin(),v2.begin(),func); }
    iterator end() const { return iterator(v1.end(),v2.end(),func); }
  };

  template <class Vector1T, class Vector2T, class FuncT>
  struct VectorSize<VectorBinaryFunc<Vector1T,Vector2T,FuncT> > {
    static const size_t value = (VectorSize<Vector1T>::value!=0)?(VectorSize<Vector1T>::value):(VectorSize<Vector2T>::value);
  };


  // *******************************************************************
  // Explicit vector expression evaluation functions.
  // *******************************************************************

  /// Forces evaluation of an arbitrary vector expression to a Vector
  /// object.
  template <class VectorT>
  Vector<typename VectorT::value_type, VectorSize<VectorT>::value>
  inline eval( VectorBase<VectorT> const& v ) {
    return v;
  }

  /// Forwarding overload for plain Vector objects.
  template <class ElemT, size_t SizeN>
  inline Vector<ElemT,SizeN> const& eval( Vector<ElemT,SizeN> const& v ) {
    return v;
  }


  // *******************************************************************
  // Vector iostream interface functions.
  // *******************************************************************

  /// Dumps a vector to a std::ostream
  template <class VectorT>
  inline std::ostream& operator<<( std::ostream& os, VectorBase<VectorT> const& v ) {
    VectorT const& vr = v.impl();
    size_t size = vr.size();
    os << "Vector" << size << '(';
    if( size > 0 ) os << vr(0);
    for( size_t i=1; i<size; ++i ) os << ',' << vr(i);
    return os << ')';
  }

  /// Dumps a transposed vector to a std::ostream
  template <class VectorT>
  inline std::ostream& operator<<( std::ostream& os, VectorTranspose<VectorT> const& v ) {
    VectorT const& vr = v.child();
    size_t size = vr.size();
    os << "Vector" << size << '(';
    if( size > 0 ) os << vr(0);
    for( size_t i=1; i<size; ++i ) os << ',' << vr(i);
    return os << ")'";
  }


  // *******************************************************************
  // Vector comparison operators and functions.
  // Note that only equality and inequality operators are provided,
  // in keeping with standard mathematical notation.  Users who want
  // particular orderings can defined those operators appropriately.
  // *******************************************************************

  /// Equality of two vectors.  Two vectors are considered equal if
  /// they have the same dimensions and their elements are all
  /// equivalent with respect to the standard c++ operator==().
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
  /// elements are equal to within the specified tolerance.  Note that
  /// if you have introduced std::equal into the global namespace then
  /// you will have to explicitly request this one as vw::math::equal().
  template <class Vector1T, class Vector2T>
  inline bool equal( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2, double epsilon = 0 ) {
    if (v1.impl().size() != v2.impl().size()) { return false; }

    typename Vector1T::const_iterator iter1 = v1.impl().begin();
    typename Vector2T::const_iterator iter2 = v2.impl().begin();
    for (; iter1 != v1.impl().end(); ++iter1, ++iter2)
      if (fabs(*iter1 - *iter2) > epsilon) { return false; }
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


  // *******************************************************************
  // Basic elementwise mathematical vector operators and functions.
  // *******************************************************************

  /// Negation of a vector.
  template <class VectorT>
  VectorUnaryFunc<VectorT, ArgNegationFunctor>
  inline operator-( VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ArgNegationFunctor>( v.impl() );
  }

  /// Negation of a transposed vector.
  template <class VectorT>
  VectorTranspose<const VectorUnaryFunc<VectorT, ArgNegationFunctor> >
  inline operator-( VectorTranspose<VectorT> const& v ) {
    return transpose(-v.child());
  }

  /// Absolute of an image.
  template <class VectorT>
  VectorUnaryFunc<VectorT, ArgAbsFunctor>
  inline abs( VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ArgAbsFunctor>( v.impl() );
  }

  /// Absolute of a transposed vector
  template <class VectorT>
  VectorTranspose<const VectorUnaryFunc<VectorT, ArgAbsFunctor> >
  inline abs( VectorTranspose<VectorT> const& v ) {
    return transpose(abs(v.child()));
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
  VectorTranspose<const VectorBinaryFunc<Vector1T, Vector2T, ArgArgSumFunctor> >
  inline operator+( VectorTranspose<Vector1T> const& v1, VectorTranspose<Vector2T> const& v2 ) {
    return transpose(v1.child()+v2.child());
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
  VectorTranspose<const VectorBinaryFunc<Vector1T, Vector2T, ArgArgDifferenceFunctor> >
  inline operator-( VectorTranspose<Vector1T> const& v1, VectorTranspose<Vector2T> const& v2 ) {
    return transpose(v1.child()-v2.child());
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
                             VectorTranspose<const VectorUnaryFunc<VectorT, ValArgProductFunctor<ScalarT> > > >::type
  inline operator*( ScalarT s, VectorTranspose<VectorT> const& v ) {
    return transpose(s*v.child());
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
                             VectorTranspose<const VectorUnaryFunc<VectorT, ArgValProductFunctor<ScalarT> > > >::type
  inline operator*( VectorTranspose<VectorT> const& v, ScalarT s ) {
    return transpose(v.child()*s);
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
                             VectorTranspose<const VectorUnaryFunc<VectorT, ArgValQuotientFunctor<ScalarT> > > >::type
  inline operator/( VectorTranspose<VectorT> const& v, ScalarT s ) {
    return transpose(v.child()/s);
  }


  // *******************************************************************
  // Elementwise vector comparison functions.
  // *******************************************************************

  /// Elementwise equality of two vectors.
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgEqualityFunctor>
  inline elem_eq( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgEqualityFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise equality of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgEqualityFunctor<ScalarT> > >::type
  inline elem_eq( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgEqualityFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise equality of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValEqualityFunctor<ScalarT> > >::type
  inline elem_eq( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValEqualityFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise inequality of two vectors.
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgInequalityFunctor>
  inline elem_neq( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgInequalityFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise inequality of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgInequalityFunctor<ScalarT> > >::type
  inline elem_neq( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgInequalityFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise inequality of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValInequalityFunctor<ScalarT> > >::type
  inline elem_neq( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValInequalityFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise less-than of two vectors.
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgLessThanFunctor>
  inline elem_lt( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgLessThanFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise less-than of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgLessThanFunctor<ScalarT> > >::type
  inline elem_lt( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgLessThanFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise less-than of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValLessThanFunctor<ScalarT> > >::type
  inline elem_lt( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValLessThanFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise less-than-or-equal-to of two vectors.
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgLessThanOrEqualFunctor>
  inline elem_lte( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgLessThanOrEqualFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise less-than-or-equal-to of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgLessThanOrEqualFunctor<ScalarT> > >::type
  inline elem_lte( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgLessThanOrEqualFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise less-than-or-equal-to of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValLessThanOrEqualFunctor<ScalarT> > >::type
  inline elem_lte( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValLessThanOrEqualFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise greater-than of two vectors.
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgGreaterThanFunctor>
  inline elem_gt( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgGreaterThanFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise greater-than of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgGreaterThanFunctor<ScalarT> > >::type
  inline elem_gt( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgGreaterThanFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise greater-than of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValGreaterThanFunctor<ScalarT> > >::type
  inline elem_gt( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValGreaterThanFunctor<ScalarT> >( v.impl(), s );
  }


  /// Elementwise greater-than-or-equal-to of two vectors.
  template <class Vector1T, class Vector2T>
  VectorBinaryFunc<Vector1T, Vector2T, ArgArgGreaterThanOrEqualFunctor>
  inline elem_gte( VectorBase<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return VectorBinaryFunc<Vector1T, Vector2T, ArgArgGreaterThanOrEqualFunctor>( v1.impl(), v2.impl() );
  }

  /// Elementwise greater-than-or-equal-to of a scalar and a vector.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ValArgGreaterThanOrEqualFunctor<ScalarT> > >::type
  inline elem_gte( ScalarT s, VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ValArgGreaterThanOrEqualFunctor<ScalarT> >( v.impl(), s );
  }

  /// Elementwise greater-than-or-equal-to of a vector and a scalar.
  template <class ScalarT, class VectorT>
  typename boost::enable_if< IsScalar<ScalarT>,
                             VectorUnaryFunc<VectorT, ArgValGreaterThanOrEqualFunctor<ScalarT> > >::type
  inline elem_gte( VectorBase<VectorT> const& v, ScalarT s ) {
    return VectorUnaryFunc<VectorT, ArgValGreaterThanOrEqualFunctor<ScalarT> >( v.impl(), s );
  }


  // *******************************************************************
  // Vector norms and similar functions.
  // *******************************************************************

  /// Vector 1-norm
  template <class VectorT>
  inline typename VectorT::value_type norm_1( VectorBase<VectorT> const& v ) {
    double result = 0.0;
    typename VectorT::const_iterator i = v.impl().begin(), end = v.impl().end();
    for( ; i != end ; ++i ) result += fabs( *i );
    return static_cast<typename VectorT::value_type>(result);
  }

  /// Square of vector 2-norm
  template <class VectorT>
  inline double norm_2_sqr( VectorBase<VectorT> const& v ) {
    double result = 0.0;
    typename VectorT::const_iterator i = v.impl().begin(), end = v.impl().end();
    for( ; i != end ; ++i ) result += (*i) * (*i);
    return static_cast<typename VectorT::value_type>(result);
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
    return static_cast<typename VectorT::value_type>(result);
  }

  /// Index of the element with the largest magnitude
  template <class VectorT>
  inline size_t index_norm_inf( VectorBase<VectorT> const& v ) {
    double maxval = -1;
    size_t index=0, result=0;
    typename VectorT::const_iterator i = v.impl().begin(), end = v.impl().end();
    for( ; i != end ; ++i, ++index ) {
      double a = fabs( *i );
      if( a > maxval ) {
        maxval = a;
        result = index;
      }
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

  /// Maximium of elements
  template <class VectorT>
  inline typename VectorT::value_type max( VectorBase<VectorT> const& v ) {
    typename VectorT::const_iterator i = v.impl().begin(), end = v.impl().end();
    if( i == end ) return typename VectorT::value_type(0);
    typename VectorT::value_type result = *(i++);
    for( ; i != end ; ++i )
      if ( *i > result )
        result = *i;
    return result;
  }

  /// Minimium of elements
  template <class VectorT>
  inline typename VectorT::value_type min( VectorBase<VectorT> const& v ) {
    typename VectorT::const_iterator i = v.impl().begin(), end = v.impl().end();
    if( i == end ) return typename VectorT::value_type(0);
    typename VectorT::value_type result = *(i++);
    for( ; i != end ; ++i )
      if ( *i < result )
        result = *i;
    return result;
  }

  /// Fill all elements with a given value
  template <class VectorT, class ValT>
  inline void fill( VectorBase<VectorT> &v, ValT const& val  ) {
    typename VectorT::iterator i = v.impl().begin(), end = v.impl().end();
    if( i == end ) return;
    for( ; i != end ; ++i ) *i = val;
  }


  // *******************************************************************
  // Assorted mathematical vector functions.
  // *******************************************************************

  /// Returns a normalized (unit) vector.
  template <class VectorT>
  VectorUnaryFunc<VectorT, ArgValQuotientFunctor<double> >
  inline normalize( VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ArgValQuotientFunctor<double> >( v.impl(), norm_2(v) );
  }

  /// Returns a homogenized vector.
  template <class VectorT>
  Vector<typename VectorT::value_type>
  inline hom( VectorBase<VectorT> const& v ) {
    Vector<typename VectorT::value_type> r;
    r.set_size( v.impl().size()+1 );
    subvector( r, 0, v.impl().size() ) = subvector( v, 0, v.impl().size() );
    r(v.impl().size()) = 1;
    return r;
  }

  /// Returns a dehomogenized (projected) vector.
  template <class VectorT>
  Vector<typename VectorT::value_type>
  inline dehom( VectorBase<VectorT> const& v ) {
    return subvector( VectorUnaryFunc<VectorT, ArgValQuotientFunctor<typename VectorT::value_type> >( v.impl(), v.impl()(v.impl().size()-1) ), 0, v.impl().size()-1 );
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
  inline operator*( VectorTranspose<Vector1T> const& v1, VectorBase<Vector2T> const& v2 ) {
    return dot_prod( v1.child(), v2 );
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

  /// Real part of of a (complex) vector.
  template <class VectorT>
  VectorUnaryFunc<VectorT, ArgRealFunctor>
  inline real( VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ArgRealFunctor>( v.impl() );
  }

  /// Imaginary part of of a (complex) vector.
  template <class VectorT>
  VectorUnaryFunc<VectorT, ArgImagFunctor>
  inline imag( VectorBase<VectorT> const& v ) {
    return VectorUnaryFunc<VectorT, ArgImagFunctor>( v.impl() );
  }

} // namespace math

  // Typedefs for commonly-used static vector types and using
  // directives for backwards compatability.
  using math::Vector;
  using math::VectorBase;
  using math::VectorProxy;
  typedef Vector<float64,2> Vector2;
  typedef Vector<float64,3> Vector3;
  typedef Vector<float64,4> Vector4;
  typedef Vector<float32,2> Vector2f;
  typedef Vector<float32,3> Vector3f;
  typedef Vector<float32,4> Vector4f;
  typedef Vector<int32,2> Vector2i;
  typedef Vector<int32,3> Vector3i;
  typedef Vector<int32,4> Vector4i;

} // namespace vw

#endif // __VW_MATH_VECTOR_H__
