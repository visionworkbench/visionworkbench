// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Core/Cache.h
/// 
/// The Vision Workbench provides a thread-safe system for caching
/// regeneratable data.  When the cache is full, the least recently
/// used object is "invalidated" to make room for new objects.
/// Invalidated objects have had the resource associated with them
/// (e.g. memory or other resources) deallocated or freed, however,
/// the object can be "regenerated" (that is, the resource is
/// regenerated automatically by the cache) when the object is next
/// accessed.

/// The vw::Cache object defined in src/vw/Core/Cache.h can be used to
/// store any resource.  For example, one common usage would be to
/// create a cache of image blocks in memory.  In this case, the cache
/// enforces a maximum memory footprint for image block storage, and
/// it regenerates the blocks (e.g. reloads them from a file on disk)
/// when necessary if a block is accessed.
///
/// Types and functions to assist cacheing regeneratable data.
///
/// The main public API is thread-safe:
///  Cache::insert(GeneratorT const&)
///  Cache::system_cache()
///  Cache::resize(size_t)
///  The entire Handle<GeneratorT> class
///
/// No other functions are guaranteed to be thread-safe.  There are
/// two levels of synchronization: one lock per cache to protect the
/// cache data structure itself, and one lock per cache line to
/// protect the m_value pointer and synchronize the (potentially very
/// expensive) generation operation.  However, the lock on the cache
/// line ends just before the generate() method is called on the
/// m_value object itself, so that object is responsible for its own
/// thread safety.
///
/// Note also that the valid() function is only useful as a heuristic:
/// there is no guarantee that the cache line won't be invalidated
/// between when the function checks the state and when you examine
/// the result.
///
#ifndef __VW_CORE_CACHE_H__
#define __VW_CORE_CACHE_H__

// Uncomment one of these to enable or disable cache debug messages
#define VW_CACHE_DEBUG(x) x
//#define VW_CACHE_DEBUG(x)

#include <typeinfo>
#include <boost/smart_ptr.hpp>

#include <vw/Core/Exception.h>
#include <vw/Core/Thread.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Core/Log.h>

#include <iostream>

namespace vw {

  // Cache contains a list of pointers to CacheLine CacheLine is
  // virtual and contains {generator,object,valid} Handle contains a
  // shared pointer to CacheLine

  // An LRU-based regeneratable-data cache
  class Cache {

    // The abstract base class for all cache line objects.
    class CacheLineBase {
      Cache& m_cache;
      CacheLineBase *m_prev, *m_next;
      const size_t m_size;
      friend class Cache;
    protected:
      Cache& cache() const { return m_cache; }
      inline void allocate() { m_cache.allocate(m_size); }
      inline void deallocate() { m_cache.deallocate(m_size); }
      inline void validate() { m_cache.validate(this); }
      inline void remove() { m_cache.remove( this ); }
      inline void deprioritize() { m_cache.deprioritize(this); }
    public:
      CacheLineBase( Cache& cache, size_t size ) : m_cache(cache), m_prev(0), m_next(0), m_size(size) {}
      virtual ~CacheLineBase() {}
      virtual inline void invalidate() { m_cache.invalidate(this); }
      virtual size_t size() const { return m_size; }
    };
    friend class CacheLineBase;

    // CacheLine<>
    template <class GeneratorT>
    class CacheLine : public CacheLineBase {
      GeneratorT m_generator;
      typename boost::shared_ptr<typename GeneratorT::value_type> m_value;
      Mutex m_mutex; // Mutex for m_value and generation of this cache line
      unsigned m_generation_count;

    public:
      CacheLine( Cache& cache, GeneratorT const& generator )
        : CacheLineBase(cache,generator.size()), m_generator(generator), m_generation_count(0)
      {
        VW_CACHE_DEBUG( vw_out(DebugMessage, "cache") << "Cache creating CacheLine " << info() << "\n"; )
        Mutex::Lock cache_lock(cache.m_mutex);
        CacheLineBase::invalidate();
      }
      
      virtual ~CacheLine() {
        Mutex::Lock cache_lock(cache().m_mutex);
        invalidate();
        VW_CACHE_DEBUG( vw_out(DebugMessage, "cache") << "Cache destroying CacheLine " << info() << "\n"; )
        remove();
      }
      
      virtual void invalidate() {
        Mutex::Lock line_lock(m_mutex);
        if( ! m_value ) return;
        VW_CACHE_DEBUG( vw_out(DebugMessage, "cache") << "Cache invalidating CacheLine " << info() << "\n"; );
        CacheLineBase::invalidate();
        CacheLineBase::deallocate();
        m_value.reset();
      }

      std::string info() {
        std::ostringstream oss;
        oss << typeid(this).name() << " " << this
            << " (size " << (int)size() << ", gen count " << m_generation_count << ")";
        return oss.str();
      }
      
      typename boost::shared_ptr<typename GeneratorT::value_type> const& value() {
        Mutex::Lock line_lock(m_mutex);
        if( !m_value ) {
          m_generation_count++;
          VW_CACHE_DEBUG( vw_out(DebugMessage, "cache") << "Cache generating CacheLine " << info() << "\n"; )
          {
            Mutex::Lock cache_lock(cache().m_mutex);
            CacheLineBase::allocate();
          }
          ScopedWatch sw((std::string("Cache ")
                          + (m_generation_count == 1 ? "generating " : "regenerating ")
                          + typeid(this).name()).c_str());
          m_value = m_generator.generate();
        }
        {
          Mutex::Lock cache_lock(cache().m_mutex);
          CacheLineBase::validate();
        }
        return m_value;
      }

      bool valid() {
        Mutex::Lock line_lock(m_mutex);
        return (bool)m_value;
      }

      void deprioritize() {
        Mutex::Lock line_lock(m_mutex);
        if( m_value ) {
          Mutex::Lock cache_lock(cache().m_mutex);
          CacheLineBase::deprioritize();
        }
      }
    };


    CacheLineBase *m_first_valid, *m_last_valid, *m_first_invalid;
    size_t m_size, m_max_size;
    Mutex m_mutex;

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
      bool valid() const { return m_line_ptr->valid(); }
      size_t size() const { return m_line_ptr->size(); }
      void reset() { m_line_ptr.reset(); }
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
    size_t max_size() { return m_max_size; }

  };

  /// Use this method to return a reference to the Vision Workbench
  /// system cache.  This cache is used by default for all new
  /// BlockImageView<>'s such as DiskImageView<>.
  Cache& vw_system_cache();

} // namespace vw

#endif  // __VW_CORE_CACHE_H__
