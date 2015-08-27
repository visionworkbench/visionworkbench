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


/// \file Core/Cache.tcc


// =============== Start class CacheLine ========================================


template <class GeneratorT>
Cache::CacheLine<GeneratorT>::CacheLine( Cache& cache, GeneratorT const& generator )
  : CacheLineBase(cache,core::detail::pointerish(generator)->size()), m_generator(generator), m_generation_count(0)
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
  if( !m_value ) return; // Not in memory, don't need to do anything.

  CacheLineBase::deallocate(); // Calls invalidate internally which redirects to the parent Cache class
  m_value.reset(); // After the base class function is done, delete the last shared pointer to the data.
}

template <class GeneratorT>
bool Cache::CacheLine<GeneratorT>::try_invalidate() {
  bool have_lock = m_mutex.try_lock();
  if ( !have_lock ) 
    return false;
  if ( !m_value ) { // Data is not loaded, unlock and quit.
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
  bool hit = (bool)m_value;
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
    m_value = core::detail::pointerish(m_generator)->generate();
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
  return m_value;
}

template <class GeneratorT>
void Cache::CacheLine<GeneratorT>::deprioritize() {
  bool exists = valid();
  if ( exists ) {
    Mutex::WriteLock line_lock(m_mutex);
    CacheLineBase::deprioritize();
  }
}



// =============== Start class Handle ========================================


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


// ============= Start class Cache ========================================================


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









