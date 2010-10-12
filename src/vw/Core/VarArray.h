// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file VarArray.h
///
/// Defines a stripped-down variable-size array class.
///
/// This is essentially a stripped-down version of std::vector.  It
/// exists solely to work around the fact that on some systems the
/// std::vector::iterator performs bounds-checking when you advance
/// the iterator rather than when you dereference it.  This gets in
/// the way of clean implementations of certain classes like MatrixCol,
/// SubMatrix, SelectColView, and so on.
///
/// This class also has supports non-copying resize() semantics.
///
/// At some point we may want to extend this to deal with allocation
/// pools and so forth, if we want to start exposing that functionality
/// in matrices and images.
///
#ifndef __VW_CORE_VARARRAY_H__
#define __VW_CORE_VARARRAY_H__

#include <algorithm>
#include <boost/smart_ptr.hpp>

namespace vw {

  template <class T>
  class VarArray {
    boost::shared_array<T> m_data;
    size_t m_size;
  public:
    VarArray() : m_data(), m_size(0) {}

    VarArray( VarArray const& other ) : m_data(new T[other.size()]), m_size(other.size()) {
      std::copy( other.begin(), other.end(), begin() );
    }

    VarArray( size_t size ) : m_data(new T[size]), m_size(size) {
      std::fill(begin(),end(),T());
    }

    template <class IterT>
    VarArray( IterT b, IterT e ) : m_data(new T[e-b]), m_size(e-b) {
      std::copy(b,e,begin());
    }

    VarArray& operator=( VarArray const& other ) {
      boost::shared_array<T> new_data( new T[other.size()] );
      std::copy( other.begin(), other.end(), new_data.get() );
      m_data = new_data;
      m_size = other.size();
      return *this;
    }

    inline T& operator[]( size_t i ) {
      return m_data[i];
    }

    inline T const& operator[]( size_t i ) const {
      return m_data[i];
    }

    typedef T* iterator;
    typedef const T* const_iterator;

    iterator begin() { return m_data.get(); }
    iterator end() { return m_data.get() + m_size; }

    const_iterator begin() const { return m_data.get(); }
    const_iterator end() const { return m_data.get() + m_size; }

    size_t size() const { return m_size; }

    void resize( size_t new_size, bool preserve = true ) {
      if( new_size == m_size ) return;
      if( new_size == 0 ) {
        m_data.reset();
        m_size = 0;
        return;
      }
      boost::shared_array<T> new_data( new T[new_size] );
      if( preserve ) {
        if( new_size > m_size ) {
          std::copy( m_data.get(), m_data.get()+m_size, new_data.get() );
          std::fill( new_data.get()+m_size, new_data.get()+new_size, T() );
        }
        else std::copy( m_data.get(), m_data.get()+new_size, new_data.get() );
      }
      else std::fill( new_data.get(), new_data.get()+new_size, T() );
      m_data = new_data;
      m_size = new_size;
    }

    void swap( VarArray& other ) {
      std::swap( m_data, other.m_data );
      std::swap( m_size, other.m_size );
    }
  };

} // namespace vw

#endif // __VW_CORE_VARARRAY_H__
