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

/// \file Core/Cache.h
/// 
/// Types and functions to assist cacheing regeneratable data.
///
#ifndef __VW_CORE_CACHE_H__
#define __VW_CORE_CACHE_H__

// Uncomment one of these to enable or disable cache debug messages
//#define VW_CACHE_DEBUG(x) x
#define VW_CACHE_DEBUG(x)

#include <boost/smart_ptr.hpp>

#include <vw/Core/Exception.h>
#include <iostream>
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
      inline void allocate() { m_cache.allocate(m_size); }
      inline void deallocate() { m_cache.deallocate(m_size); }
      inline void validate() { m_cache.validate(this); }
    public:
      CacheLineBase( Cache& cache, size_t size ) : m_cache(cache), m_prev(0), m_next(0), m_size(size) { m_cache.invalidate( this ); }
      virtual ~CacheLineBase() { m_cache.deallocate(m_size); m_cache.remove( this ); }
      virtual inline void invalidate() { m_cache.invalidate(this); }
      virtual size_t size() const { return m_size; }
      void deprioritize() { m_cache.deprioritize(this); }
    };
    friend class CacheLineBase;

    // CacheLine<>
    template <class GeneratorT>
    class CacheLine : public CacheLineBase {
      GeneratorT m_generator;
      typename boost::shared_ptr<typename GeneratorT::value_type> m_value;
    public:
      CacheLine( Cache& cache, GeneratorT const& generator )
        : CacheLineBase(cache,generator.size()), m_generator(generator)
      {
        VW_CACHE_DEBUG( vw_out(VerboseDebugMessage) << "Cache creating CacheLine " << this << std::endl; )
      }
      
      virtual ~CacheLine() {
        VW_CACHE_DEBUG( vw_out(VerboseDebugMessage) << "Cache destroying CacheLine " << this << std::endl; )
      }
      
      virtual void invalidate() {
        VW_CACHE_DEBUG( vw_out(VerboseDebugMessage) << "Cache invalidating CacheLine " << this << std::endl; )
        CacheLineBase::invalidate();
        m_value.reset();
        CacheLineBase::deallocate();
      }
      
      typename boost::shared_ptr<typename GeneratorT::value_type> value() {
        VW_CACHE_DEBUG( vw_out(VerboseDebugMessage) << "Cache accessing CacheLine " << this << std::endl; )
        if( !m_value ) {
          CacheLineBase::allocate();
          m_value = m_generator.generate();
        }
        CacheLineBase::validate();
        return m_value;
      }

      bool valid() const { return m_value; }

      void deprioritize() { if( valid() ) CacheLineBase::deprioritize(); }
    };

    CacheLineBase *m_first_valid, *m_last_valid, *m_first_invalid;
    size_t m_size, m_max_size;

    void allocate( size_t size );
    void deallocate( size_t size );
    void validate( CacheLineBase *line );
    void invalidate( CacheLineBase *line );
    void remove( CacheLineBase *line );
    void deprioritize( CacheLineBase *line );

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
      typename GeneratorT::value_type const& operator*() const {
        VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
        return *(m_line_ptr->value());
      }
      operator boost::shared_ptr<typename GeneratorT::value_type>() const {
        VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
        return m_line_ptr->value();
      }
      size_t size() const { return m_line_ptr->size(); }
      bool valid() const { return m_line_ptr->valid(); }
      void deprioritize() const { return m_line_ptr->deprioritize(); }
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
