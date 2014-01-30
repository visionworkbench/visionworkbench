// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
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

#include <vw/Core/Exception.h>
#include <vw/Core/Thread.h>
#include <vw/Core/Log.h>
#include <vw/Core/FundamentalTypes.h>

#include <typeinfo>
#include <sstream>
#include <stddef.h>
#include <string>

#include <boost/smart_ptr/shared_ptr.hpp>

namespace vw {
namespace core {
namespace detail {

  template <typename T>
  T* pointerish(boost::shared_ptr<T> x) {return x.get();}

  template <typename T>
  T* pointerish(T& x) {return &x;}

  template <typename T>
  struct GenValue {
    typedef typename T::value_type type;
  };

  template <typename T>
  struct GenValue<boost::shared_ptr<T> > {
    typedef typename T::value_type type;
  };
}}} // namespace vw::core::detail

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
      inline void allocate() { m_cache.allocate(m_size, this); }
      inline void deallocate() { m_cache.deallocate(m_size, this); }
      inline void validate() { m_cache.validate(this); }
      inline void remove() { m_cache.remove(this); }
      inline void deprioritize() { m_cache.deprioritize(this); }
    public:
      CacheLineBase( Cache& cache, size_t size ) : m_cache(cache), m_prev(0), m_next(0), m_size(size) {}
      virtual ~CacheLineBase() {}
      virtual inline void invalidate() { m_cache.invalidate(this); }
      virtual inline bool try_invalidate() { m_cache.invalidate(this); return true; }
      virtual size_t size() const { return m_size; }
    };
    friend class CacheLineBase;

    // CacheLine<>
    //
    // Always follow the order of mutexs is:
    // ACQUIRE LINE FIRST
    // ACQUIRE CACHE's LINE MGMT SECOND
    template <class GeneratorT>
    class CacheLine : public CacheLineBase {
      GeneratorT m_generator;
      typedef typename boost::shared_ptr<typename core::detail::GenValue<GeneratorT>::type> value_type;
      value_type m_value;
      Mutex m_mutex; // Mutex for m_value and generation of this cache line
      uint64 m_generation_count;

    public:
      CacheLine( Cache& cache, GeneratorT const& generator )
        : CacheLineBase(cache,core::detail::pointerish(generator)->size()), m_generator(generator), m_generation_count(0)
      {
        VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache creating CacheLine " << info() << "\n"; )
        CacheLineBase::invalidate();
      }

      virtual ~CacheLine() {
        VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache destroying CacheLine " << info() << "\n"; )
        invalidate();
        remove();
      }

      virtual void invalidate() {
        VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache invalidating CacheLine " << info() << "\n"; );

        Mutex::WriteLock line_lock(m_mutex);
        if( !m_value ) return;

        CacheLineBase::deallocate(); // Calls invalidate internally
        m_value.reset();
      }

      virtual bool try_invalidate() {
        bool have_lock = m_mutex.try_lock();
        if ( !have_lock ) return false;
        if ( !m_value ) {
          m_mutex.unlock();
          return true;
        }

        VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache invalidating CacheLine " << info() << "\n"; );
        CacheLineBase::deallocate(); // Calls invalidate internally
        m_value.reset();

        m_mutex.unlock();
        return true;
      }

      std::string info() {
        Mutex::ReadLock line_lock(m_mutex);
        std::ostringstream oss;
        oss << typeid(this).name() << " " << this
            << " (size " << (int)size() << ", gen count " << m_generation_count << ")";
        return oss.str();
      }

      // This grabs a lock and never releases it! User must unlock
      // themselves because we are passing them a pointer to an object
      // that another thread could delete.
      value_type const& value() {
        m_mutex.lock_shared();
        bool hit = (bool)m_value;
        {
          { // This should be abstracted into a call
            Mutex::WriteLock cache_lock( cache().m_stats_mutex );
            if (hit)
              cache().m_hits++;
            else
              cache().m_misses++;
          }
        }
        if( !hit ) {
          VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache generating CacheLine " << info() << "\n"; );
          m_mutex.unlock_shared(); // Loose shared
          m_mutex.lock_upgrade();  // Get upgrade status
          m_mutex.unlock_upgrade_and_lock(); // Get exclusive access
          CacheLineBase::allocate(); // Call validate internally

          m_generation_count++;
          m_value = core::detail::pointerish(m_generator)->generate();
          m_mutex.unlock_and_lock_upgrade();
          m_mutex.unlock_upgrade_and_lock_shared();
        }
        return m_value;
      }

      void release() {
        m_mutex.unlock_shared();
      }

      bool valid() {
        Mutex::WriteLock line_lock(m_mutex);
        return m_value;
      }

      void deprioritize() {
        bool exists = valid();
        if ( exists ) {
          Mutex::WriteLock line_lock(m_mutex);
          CacheLineBase::deprioritize();
        }
      }
    };

    CacheLineBase *m_first_valid, *m_last_valid, *m_first_invalid;
    size_t m_size, m_max_size, m_num_warnings;
    RecursiveMutex m_line_mgmt_mutex;
    Mutex m_stats_mutex;
    volatile vw::uint64 m_hits, m_misses, m_evictions;

    void allocate( size_t size, CacheLineBase *line );
    void deallocate( size_t size, CacheLineBase *line );
    void validate( CacheLineBase *line );
    void invalidate( CacheLineBase *line );
    void remove( CacheLineBase *line );
    void deprioritize( CacheLineBase *line );

  public:

    // Handle<>
    template <class GeneratorT>
    class Handle {
      boost::shared_ptr<CacheLine<GeneratorT> > m_line_ptr;
      mutable bool m_is_locked;
    public:
      typedef typename core::detail::GenValue<GeneratorT>::type value_type;

      Handle() : m_is_locked(false) {}
      Handle( boost::shared_ptr<CacheLine<GeneratorT> > line_ptr ) : m_line_ptr(line_ptr), m_is_locked(false) {}
      ~Handle() {
        if (m_is_locked)
          m_line_ptr->release();
      }
      boost::shared_ptr<value_type> operator->() const {
        VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
        m_is_locked = true;
        return m_line_ptr->value();
      }
      value_type const& operator*() const {
        VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
        m_is_locked = true;
        return *(m_line_ptr->value());
      }
      operator boost::shared_ptr<value_type>() const {
        VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
        m_is_locked = true;
        return m_line_ptr->value();
      }
      void release() const {
        m_is_locked = false;
        m_line_ptr->release();
      }
      bool valid() const {
        VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
        return m_line_ptr->valid();
      }
      size_t size() const {
        VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
        return m_line_ptr->size();
      }
      void reset() { m_line_ptr.reset(); }
      void deprioritize() const {
        VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
        return m_line_ptr->deprioritize();
      }
      bool attached() const {
        return (bool)m_line_ptr;
      }
    };

    Cache( size_t max_size ) :
      m_first_valid(0), m_last_valid(0), m_first_invalid(0),
      m_size(0), m_max_size(max_size), m_num_warnings(0),
      m_hits(0), m_misses(0), m_evictions(0) {}

    template <class GeneratorT>
    Handle<GeneratorT> insert( GeneratorT const& generator ) {
      boost::shared_ptr<CacheLine<GeneratorT> > line( new CacheLine<GeneratorT>( *this, generator ) );
      VW_ASSERT( line, NullPtrErr() << "Error creating new cache line!" );
      return Handle<GeneratorT>( line );
    }

    void resize( size_t size );
    size_t max_size();

    // Statistics functions
    uint64 hits();
    uint64 misses();
    uint64 evictions();
    void clear_stats();
  };
} // namespace vw

#endif  // __VW_CORE_CACHE_H__
