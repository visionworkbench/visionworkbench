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
/// data.  When the cache is full, the least recently
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
/// Types and functions to assist caching data.
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
#include <stddef.h>
#include <string>

#include <boost/smart_ptr/shared_ptr.hpp>

namespace vw {
namespace core {
namespace detail {

  // These template tools are used to extract a desired output type for several
  //  possible cases of the input type.
  
  template <typename T>
  T* getPtr(boost::shared_ptr<T> x) {return x.get();}

  template <typename T>
  T* getPtr(T& x) {return &x;}

  template <typename T>
  struct GenValue {
    typedef typename T::value_type type;
  };

  template <typename T>
  struct GenValue<boost::shared_ptr<T> > {
    typedef typename T::value_type type;
  };
}}} // namespace vw::core::detail_changed

namespace vw {

  //TODO: This class design is a tangled mess!

  // Cache contains a list of pointers to CacheLine CacheLine is
  // virtual and contains {generator,object,valid} Handle contains a
  // shared pointer to CacheLine

  /// An LRU-based data cache
  /**
    - This class contains three pointers (*m_first_valid, *m_last_valid,
      *m_first_invalid) which keep track of two double-linked lists, the valid
      list and the invalid list.  
    - Each list is made up of CacheLine objects, each each CacheLine object has
      m_prev and m_next member variables which are used to maintain the lists.
    - The four private functions validate(), invalidate(), remove(),
      deprioritize() rearrange the position of CacheLine objects in the valid
      and invalid lists.

    - The Cache class itself does not directly allocate or free any memory.  It
      manages the two lists, monitors total reported memory usage, and calls
      functions on the CacheLine objects.  It also records cache hit and miss
      statistics.
    - The CacheLine class is where objects are created and destroyed (using
      smart pointers and the provided GeneratorT class))

    - TODO: The Cache class should support different allocation strategies
      besides just deleting the oldest resource every time!

    User interface:
    - Call insert() to add a new GeneratorT object (internally wrapped in a
      CacheLine object) to the Cache and you will get out a Handle object.

  */
  class Cache {
  // Do some forward declarations so we can put the public interface first
  private:
    class CacheLineBase;
    template <class GeneratorT> class CacheLine;
  public:
    template <class GeneratorT> class Handle;
    
    // Cache public functions

    /// Constructor
    inline Cache( size_t max_size );

    /// Wrap a GeneratorT in a CacheLine in a Handle object and return it.
    /// - By creating the CacheLine object it is automatically registered with the Cache object.
    ///   Retrieving the value of the CacheLine object will cause it to be added to the Cache.
    template <class GeneratorT>
    Handle<GeneratorT> insert( GeneratorT const& generator );

    void   resize( size_t size ); ///< Change the maximum size in bytes of the Cache.
    size_t max_size();            ///< Return the maximum permissible size in bytes.
 
    // Statistics functions to query and clear hit, miss, and eviction counts.
    uint64 hits       () { Mutex::ReadLock  cache_lock(m_stats_mutex);  return m_hits;      }
    uint64 misses     () { Mutex::ReadLock  cache_lock(m_stats_mutex);  return m_misses;    }
    uint64 evictions  () { Mutex::ReadLock  cache_lock(m_stats_mutex);  return m_evictions; }
    void   clear_stats() { Mutex::WriteLock cache_lock(m_stats_mutex);  m_hits=m_misses=m_evictions = 0; }
    
    /// Interface class for safe user access to CacheLine objects.
    template <class GeneratorT>
    class Handle {
      boost::shared_ptr<CacheLine<GeneratorT>> m_line_ptr;
      mutable bool m_is_locked;

    public:
      typedef typename core::detail::GenValue<GeneratorT>::type value_type;

      /// Default constructor
      Handle() : m_is_locked(false) {}
      
      /// Constructor with a CacheLine object
      Handle( boost::shared_ptr<CacheLine<GeneratorT> > line_ptr ) : m_line_ptr(line_ptr), m_is_locked(false) {}

      /// Destructor - release the Cacheline object
      ~Handle();
      
      boost::shared_ptr<value_type> operator->() const; ///< Get pointer   to the Cacheline's value
      value_type const&             operator* () const; ///< Get const ref to the Cacheline's value
      
      // TODO: Can we delete this?
      /// Define conversion to the underlying data type
      operator boost::shared_ptr<value_type>() const;
      
      // These functions just redirect to the underlying Cacheline object
      void   release     () const; ///< Release the shared mutex to the underlying data.
      bool   valid       () const; ///< Return true if the data is in memory.
      size_t size        () const; ///< Return the size in bytes of the underlying data.
      void   reset       ();       ///< Disconnect the handle from the underlying data.
      void   deprioritize() const; ///< Send the underlying data to the front of the "next to free" list.
      bool   attached    () const; ///< Return true if there is a wrapped Cacheline object.
    }; // End class Handle
    
  private:

    // Cache class private variables
    CacheLineBase      *m_first_valid, 
                       *m_last_valid, 
                       *m_first_invalid;
    size_t              m_size,     ///< Currently loaded size in bytes
                        m_max_size; ///< Maximum permissible size in bytes
    RecursiveMutex      m_line_mgmt_mutex; ///< Mutex for adjusting the CacheLineBase pointers above.
    Mutex               m_stats_mutex;     ///< Separate mutex for the statistics variables below.
    volatile vw::uint64 m_hits, m_misses, m_evictions, ///< Cache statistics
                        m_last_size; ///< Record the last size at which we printed a size warning to screen!

    // Cache class private functions
    
    /// Call validate() on the line, increment m_size, and then clear up old CacheLine objects
    /// if we went over the size limit.
    void allocate  ( size_t size, CacheLineBase *line );
    
    /// Call invalidate() on the line then decrement m_size.
    void deallocate( size_t size, CacheLineBase *line );
    
    void validate    ( CacheLineBase *line ); ///< Move the cache line to the top of the valid list.
    void invalidate  ( CacheLineBase *line ); ///< Move the cache line to the top of the invalid list.
    void remove      ( CacheLineBase *line ); ///< Remove the cache line from the cache lists.
    void deprioritize( CacheLineBase *line ); ///< Move the cache line to the bottom of the valid list.

    //TODO: Why does CacheLineBase exist?  Do we need a base class?  It is not used outside this file.
    
    /// The abstract base class for all cache line objects.
    /// - This just redirects everything to a referenced Cache object.
    class CacheLineBase {
    private:
      /// Reference to parent Cache object
      Cache& m_cache;
      /// These are used to form an ordered linked list of CacheLine objects
      CacheLineBase *m_prev, *m_next; 
      /// Size in bytes of the CacheLine data object.
      const size_t m_size;
      friend class Cache;
      
    protected:
      Cache& cache() const { return m_cache; }
      
      inline void allocate    () { m_cache.allocate  (m_size, this); }
      inline void deallocate  () { m_cache.deallocate(m_size, this); }
      inline void validate    () { m_cache.validate    (this); }
      inline void remove      () { m_cache.remove      (this); }
      inline void deprioritize() { m_cache.deprioritize(this); }
      
    public:
      CacheLineBase( Cache& cache, size_t size ) : m_cache(cache), 
                                                   m_prev(0), m_next(0), 
                                                   m_size(size) {}
      virtual ~CacheLineBase() {}
      
      virtual inline void   invalidate    ()       { m_cache.invalidate(this); }
      virtual inline bool   try_invalidate()       { m_cache.invalidate(this); return true; }
      virtual inline size_t size          () const { return m_size; }
    }; // End class CacheLineBase
    friend class CacheLineBase; // Make this a friend of the Cache class

    /// Private class to wrap a data generator object and keep a pointer to the generated data.
    // Always follow the order of mutexs is:
    // ACQUIRE LINE FIRST
    // ACQUIRE CACHE's LINE MGMT SECOND
    template <class GeneratorT>
    class CacheLine : public CacheLineBase {
    
      typedef typename boost::shared_ptr<typename core::detail::GenValue<GeneratorT>::type> value_type;
      GeneratorT m_generator;
      value_type m_value;
      Mutex      m_mutex; // Mutex for m_value and generation of this cache line
      uint64     m_generation_count;

    public:
      /// Constructor
      CacheLine( Cache& cache, GeneratorT const& generator );

      virtual ~CacheLine();

      /// If data is in memory deallocate it, otherwise do nothing.  Blocks.
      /// - Note that this does not call Base::invalidate()!
      /// - Maybe this function should have been called something else?
      virtual void invalidate();

      /// Non-blocking version of invalidate.
      virtual bool try_invalidate();

      /// Print some information about this object.
      std::string info();

      /// Load the data if needed and then return a pointer to it.
      /// - This grabs a (shared) lock and never releases it! User must unlock
      ///   themselves because we are passing them a pointer to an object
      ///   that another thread could delete.
      value_type const& value();
      
      /// Release all access to the data.
      void release();
      
      /// Check whether the data is currently loaded into memory.
      bool valid();

      /// Call deprioritize from the Cache class
      void deprioritize();
    }; // End class Cacheline

  }; // End class Cache

// Start class CacheLine


template <class GeneratorT>
Cache::CacheLine<GeneratorT>::CacheLine( Cache& cache, GeneratorT const& generator )
  : CacheLineBase(cache,core::detail::getPtr(generator)->size()), m_generator(generator), m_generation_count(0)
{
  VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache creating CacheLine " << info() << "\n"; )
  CacheLineBase::invalidate(); // Move to the start of the Cache class invalid list.
}

template <class GeneratorT>
Cache::CacheLine<GeneratorT>::~CacheLine() {
  VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache destroying CacheLine " << info() << "\n"; )
  invalidate(); // Clean up the allocated data.
  remove();
}

template <class GeneratorT>
void Cache::CacheLine<GeneratorT>::invalidate() {
  VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache invalidating CacheLine " << info() << "\n"; );

  Mutex::WriteLock line_lock(m_mutex); // Grab a lock until the function exits.
  if (m_value.get() == NULL) return; // Not in memory, don't need to do anything.

  CacheLineBase::deallocate(); // Calls invalidate internally which redirects to the parent Cache class
  m_value.reset(); // After the base class function is done, delete the last shared pointer to the data.
}

template <class GeneratorT>
bool Cache::CacheLine<GeneratorT>::try_invalidate() {
  bool have_lock = m_mutex.try_lock();
  if ( !have_lock ) 
    return false;
  if (m_value.get() == NULL) { // Data is not loaded, unlock and quit.
    m_mutex.unlock();
    return true;
  }

  VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache invalidating CacheLine " << info() << "\n"; );
  CacheLineBase::deallocate(); // Calls invalidate internally
  m_value.reset();

  m_mutex.unlock();
  return true;
}

template <class GeneratorT>
std::string Cache::CacheLine<GeneratorT>::info() {
  Mutex::ReadLock line_lock(m_mutex);
  std::ostringstream oss;
  oss << typeid(this).name() << " " << this
      << " (size " << (int)size() << ", gen count " << m_generation_count << ")";
  return oss.str();
}

template <class GeneratorT>
typename Cache::CacheLine<GeneratorT>::value_type const& Cache::CacheLine<GeneratorT>::value() {

  m_mutex.lock_shared(); // Grab a shared lock
  bool hit = (m_value.get() != NULL);
  { // Grab a temporary mutex to update our cache statistics
    { // TODO: This should be abstracted into a call
      Mutex::WriteLock cache_lock( cache().m_stats_mutex );
      if (hit)
        cache().m_hits++;
      else
        cache().m_misses++;
    }
  }
  if( !hit ) { // Then we need to load the data into memory.
    VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache generating CacheLine " << info() << "\n"; );
    m_mutex.unlock_shared(); // Release shared
    m_mutex.lock_upgrade();  // Get upgrade status
    m_mutex.unlock_upgrade_and_lock(); // Upgrade to exclusive access
    CacheLineBase::allocate(); // Call validate internally

    //TODO: Why allocate and then generate?
    m_generation_count++; // Update stats
    m_value = core::detail::getPtr(m_generator)->generate();
    // Downgrade from exclusive access down to shared access
    m_mutex.unlock_and_lock_upgrade();
    m_mutex.unlock_upgrade_and_lock_shared();
  }
  return m_value; // Now the data is in memory, return a pointer to it.
}

template <class GeneratorT>
void Cache::CacheLine<GeneratorT>::release() {
  m_mutex.unlock_shared();
}

template <class GeneratorT>
bool Cache::CacheLine<GeneratorT>::valid() {
  Mutex::WriteLock line_lock(m_mutex);
  return (m_value.get() != NULL);
}

template <class GeneratorT>
void Cache::CacheLine<GeneratorT>::deprioritize() {
  bool exists = valid();
  if ( exists ) {
    Mutex::WriteLock line_lock(m_mutex);
    CacheLineBase::deprioritize();
  }
}

template <class GeneratorT>
Cache::Handle<GeneratorT>::~Handle() {
  if (m_is_locked)
    m_line_ptr->release();
}

template <class GeneratorT>
boost::shared_ptr<typename Cache::Handle<GeneratorT>::value_type> Cache::Handle<GeneratorT>::operator->() const {
  VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
  m_is_locked = true;
  return m_line_ptr->value();
}

template <class GeneratorT>
typename Cache::Handle<GeneratorT>::value_type const& Cache::Handle<GeneratorT>::operator*() const {
  VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
  m_is_locked = true;
  return *(m_line_ptr->value());
}

// TODO: Can we delete this?
template <class GeneratorT>
Cache::Handle<GeneratorT>::operator boost::shared_ptr<value_type>() const {
  VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
  m_is_locked = true;
  return m_line_ptr->value();
}

template <class GeneratorT>
void Cache::Handle<GeneratorT>::release() const {
  m_is_locked = false;
  m_line_ptr->release();
}

template <class GeneratorT>
bool Cache::Handle<GeneratorT>::valid() const {
  VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
  return m_line_ptr->valid();
}

template <class GeneratorT>
size_t Cache::Handle<GeneratorT>::size() const {
  VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
  return m_line_ptr->size();
}

template <class GeneratorT>
void Cache::Handle<GeneratorT>::reset() { 
  m_line_ptr.reset(); 
}

template <class GeneratorT>
void Cache::Handle<GeneratorT>::deprioritize() const {
  VW_ASSERT( m_line_ptr, NullPtrErr() << "Invalid cache handle!" );
  return m_line_ptr->deprioritize();
}

template <class GeneratorT>
bool Cache::Handle<GeneratorT>::attached() const {
  return (bool)m_line_ptr;
}

Cache::Cache( size_t max_size ) :
  m_first_valid(0), m_last_valid(0), m_first_invalid(0),
  m_size(0), m_max_size(max_size),
  m_hits(0), m_misses(0), m_evictions(0), m_last_size(0) {
}

template <class GeneratorT>
Cache::Handle<GeneratorT> Cache::insert( GeneratorT const& generator ) {
  boost::shared_ptr<CacheLine<GeneratorT> > line( new CacheLine<GeneratorT>( *this, generator ) );
  VW_ASSERT( line, NullPtrErr() << "Error creating new cache line!" );
  return Handle<GeneratorT>( line );
}
  
} // namespace vw

#endif  // __VW_CORE_CACHE_H__
