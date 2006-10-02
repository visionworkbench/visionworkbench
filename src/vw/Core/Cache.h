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

/// \file Core/Cache.h
/// 
/// Types and functions to assist cacheing regeneratable data.
///
#ifndef __VW_CORE_CACHE_H__
#define __VW_CORE_CACHE_H__

#include <boost/smart_ptr.hpp>

#include <vw/Core/Exception.h>

namespace vw {

  // Cache contains a list of pointers to CacheLine
  // CacheLine is virtual and contains {generator,object,valid}
  // Handle contains a shared pointer to CacheLine

  // An LRU-based regeneratable-data cache
  class Cache {

    // CacheLineBase
    class CacheLineBase {
      Cache& m_cache;
      CacheLineBase *m_prev, *m_next;
      size_t m_size;
      friend class Cache;
    protected:
      inline void allocate( size_t size ) { m_cache.allocate(m_size=size); }
      inline void deallocate() { m_cache.deallocate(m_size); }
      inline void validate() { m_cache.validate(this); }
    public:
      CacheLineBase( Cache& cache ) : m_cache(cache), m_prev(0), m_next(0) { m_cache.invalidate( this ); }
      virtual ~CacheLineBase() { m_cache.remove( this ); }
      virtual inline void invalidate() { m_cache.invalidate(this); }
    };
    friend class CacheLineBase;

    // CacheLine<>
    template <class GeneratorT>
    class CacheLine : public CacheLineBase {
      GeneratorT m_generator;
      typename boost::shared_ptr<typename GeneratorT::value_type> m_value;
    public:
      CacheLine( Cache& cache, GeneratorT const& generator )
        : CacheLineBase(cache), m_generator(generator) {}
      
      virtual ~CacheLine() {}
      
      virtual void invalidate() {
        CacheLineBase::invalidate();
        m_value.reset();
        CacheLineBase::deallocate();
      }
      
      typename boost::shared_ptr<typename GeneratorT::value_type> value() {
        if( !m_value ) {
          CacheLineBase::allocate( m_generator.size() );
          m_value = m_generator.generate();
        }
        CacheLineBase::validate();
        return m_value;
      }
    };

    CacheLineBase *m_first_valid, *m_last_valid, *m_first_invalid;
    size_t m_size, m_max_size;

    void allocate( size_t size );
    void deallocate( size_t size );
    void validate( CacheLineBase *line );
    void invalidate( CacheLineBase *line );
    void remove( CacheLineBase *line );

  public:

    // Handle<>
    template <class GeneratorT>
    class Handle {
      boost::shared_ptr<CacheLine<GeneratorT> > m_line_ptr;
    public:
      Handle() {}
      Handle( boost::shared_ptr<CacheLine<GeneratorT> > line_ptr ) : m_line_ptr(line_ptr) {}
      boost::shared_ptr<typename GeneratorT::value_type> operator->() const { 
        VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
        return m_line_ptr->value();
      }
      
    };

    Cache( size_t max_size ) : m_first_valid(0), m_last_valid(0), m_first_invalid(0), m_size(0), m_max_size(max_size) {}

    template <class GeneratorT>
    Handle<GeneratorT> insert( GeneratorT const& generator ) {
      boost::shared_ptr<CacheLine<GeneratorT> > line( new CacheLine<GeneratorT>( *this, generator ) );
      VW_ASSERT( line, NullPtrErr() << "Error creating new cache line!" );
      return Handle<GeneratorT>( line );
    }

    void resize( size_t size );

    static Cache& system_cache();
  };

} // namespace vw

#endif  // __VW_CORE_CACHE_H__
