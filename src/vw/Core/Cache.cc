// __BEGIN_LICENSE__
//  Copyright (c) 2006-2012, United States Government as represented by the
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


/// \file Core/Cache.cc
///
/// Types and functions to assist cacheing regeneratable data.
///
#include <vw/Core/Cache.h>
#include <vw/Core/Debugging.h>

void vw::Cache::allocate( size_t size, CacheLineBase* line ) {
  // WARNING! YOU CAN NOT HOLD THE CACHE MUTEX AND THEN CALL
  // INVALIDATE. That's a line -> cache -> line mutex hold. A
  // deadlock!

  // This step also implies the need to call validate.

  uint64 local_evictions = 0;
  size_t local_size, local_max_size;
  CacheLineBase* local_last_valid;
  size_t invalidate_loop_count = 0;
  { // Locally buffer variables that should be accessed behind mutex
    RecursiveMutex::Lock cache_lock( m_line_mgmt_mutex );
    validate( line ); // Call here to insure that last_valid is not us!
    m_size += size;
    local_size = m_size;
    local_max_size = m_max_size;
    local_last_valid = m_last_valid;
    VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache allocated " << size << " bytes (" << m_size << " / " << m_max_size << " used)" << "\n"; );
  }
  while ( local_size > local_max_size ) {
    if ( local_last_valid == line ) {
      // We're trying to deallocate self! Call validate again to get
      // our line back on top since some other threads have moved us.
      RecursiveMutex::Lock cache_lock( m_line_mgmt_mutex );
      validate( line );
      local_last_valid = m_last_valid;
    }
    if ( local_last_valid == line || !local_last_valid ) {
      VW_OUT(WarningMessage, "console") << "Warning: Cached object (" << size << ") larger than requested maximum cache size (" << m_max_size << "). Current Size = " << m_size << "\n";
      VW_OUT(WarningMessage, "cache") << "Warning: Cached object (" << size << ") larger than requested maximum cache size (" << m_max_size << "). Current Size = " << m_size << "\n";
      break;
    }
    bool invalidated = local_last_valid->try_invalidate();
    if (invalidated) {
      local_evictions++;
      { // Update local buffer by grabbing cache buffer
        RecursiveMutex::Lock cache_lock( m_line_mgmt_mutex );
        local_size = m_size;
        local_last_valid = m_last_valid;
      }
    } else {
      RecursiveMutex::Lock cache_lock( m_line_mgmt_mutex );
      local_last_valid = local_last_valid->m_prev;
      local_size = m_size;
      if (!local_last_valid) {
        invalidate_loop_count++;
        if ( invalidate_loop_count == 3 ) {
          // We tried hard to deallocate enough to do our
          // allocation. However there is a time to give up and just
          // move on. In practice, even with this cop out, we do stay
          // pretty close to our allocations amount.
          break;
        } {
          RecursiveMutex::Lock cache_lock( m_line_mgmt_mutex );
          local_last_valid = m_last_valid;
        }
      }
    }
  }
  {
    Mutex::WriteLock cache_lock( m_stats_mutex );
    m_evictions += local_evictions;
  }
}

void vw::Cache::resize( size_t size ) {
  // WARNING! YOU CAN NOT HOLD THE CACHE MUTEX AND THEN CALL
  // INVALIDATE. That's a line -> cache -> line mutex hold. A
  // deadlock!
  size_t local_size, local_max_size;
  CacheLineBase* local_last_valid;
  { // Locally buffer variables that require Cache Mutex
    RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);
    m_max_size = size;
    local_size = m_size;
    local_max_size = m_max_size;
    local_last_valid = m_last_valid;
  }
  while ( local_size > local_max_size ) {
    VW_ASSERT( local_last_valid, LogicErr() << "Cache is empty but has nonzero size!" );
    m_last_valid->invalidate(); // Problem ( probably grabs a line's mutex too )
    { // Update local buffer by grabbing cache buffer
      RecursiveMutex::Lock cache_lock( m_line_mgmt_mutex );
      local_size = m_size;
      local_last_valid = m_last_valid;
    }
  }
}

size_t vw::Cache::max_size() {
  RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);
  return m_max_size;
}

void vw::Cache::deallocate( size_t size, CacheLineBase *line ) {
  RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);

  // This call implies the need to call invalidate
  invalidate( line );

  m_size -= size;
  VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache deallocated " << size << " bytes (" << m_size << " / " << m_max_size << " used)" << "\n"; )
}

// Move the cache line to the top of the valid list.
void vw::Cache::validate( CacheLineBase *line ) {
  RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);
  if( line == m_first_valid ) return;
  if( line == m_last_valid ) m_last_valid = line->m_prev;
  if( line == m_first_invalid ) m_first_invalid = line->m_next;
  if( line->m_next ) line->m_next->m_prev = line->m_prev;
  if( line->m_prev ) line->m_prev->m_next = line->m_next;
  line->m_next = m_first_valid;
  line->m_prev = 0;
  if( m_first_valid ) m_first_valid->m_prev = line;
  m_first_valid = line;
  if( ! m_last_valid ) m_last_valid = line;
}

// Move the cache line to the top of the invalid list.
void vw::Cache::invalidate( CacheLineBase *line ) {
  RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);
  if( line == m_first_valid ) m_first_valid = line->m_next;
  if( line == m_last_valid ) m_last_valid = line->m_prev;
  if( line->m_next ) line->m_next->m_prev = line->m_prev;
  if( line->m_prev ) line->m_prev->m_next = line->m_next;
  line->m_next = m_first_invalid;
  line->m_prev = 0;
  if( m_first_invalid ) m_first_invalid->m_prev = line;
  m_first_invalid = line;
}

// Remove the cache line from the cache lists.
void vw::Cache::remove( CacheLineBase *line ) {
  RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);
  if( line == m_first_valid ) m_first_valid = line->m_next;
  if( line == m_last_valid ) m_last_valid = line->m_prev;
  if( line == m_first_invalid ) m_first_invalid = line->m_next;
  if( line->m_next ) line->m_next->m_prev = line->m_prev;
  if( line->m_prev ) line->m_prev->m_next = line->m_next;
  line->m_next = line->m_prev = 0;
}

// Move the cache line to the bottom of the valid list.
void vw::Cache::deprioritize( CacheLineBase *line ) {
  RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);
  if( line == m_last_valid ) return;
  if( line == m_first_valid ) m_first_valid = line->m_next;
  if( line->m_next ) line->m_next->m_prev = line->m_prev;
  if( line->m_prev ) line->m_prev->m_next = line->m_next;
  line->m_prev = m_last_valid;
  line->m_next = 0;
  m_last_valid->m_next = line;
  m_last_valid = line;
}

// Statistics request methods
vw::uint64 vw::Cache::hits() {
  Mutex::ReadLock cache_lock( m_stats_mutex );
  return m_hits;
}

vw::uint64 vw::Cache::misses() {
  Mutex::ReadLock cache_lock( m_stats_mutex );
  return m_misses;
}

vw::uint64 vw::Cache::evictions() {
  Mutex::ReadLock cache_lock( m_stats_mutex );
  return m_evictions;
}

void vw::Cache::clear_stats() {
  Mutex::WriteLock cache_lock( m_stats_mutex );
  m_hits = m_misses = m_evictions = 0;
}
