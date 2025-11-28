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

/// \file Core/Cache.cc
///
/// Types and functions to assist caching data that can be regenerated
/// on-the-fly.
///
#include <vw/Core/Cache.h>

// Note that this function does not actually load the data,
// it is up to the calling function to do that.
void vw::Cache::allocate( size_t size, CacheLineBase* line ) {

  // Put the current cache line at the top of the list (so the most
  // recently used). If the cache size is beyond the storage limit,
  // de-allocate the least recently used elements.

  // Note: Doing allocation implies the need to call validate.

  // WARNING! YOU CAN NOT HOLD THE CACHE MUTEX AND THEN CALL
  // INVALIDATE. That's a line -> cache -> line mutex hold. A deadlock!

  // The lock below is recursive, so if a resource is locked by a
  // thread, it can still be accessed by this thread, but not by others.
  RecursiveMutex::Lock cache_lock( m_line_mgmt_mutex );

  uint64 local_evictions = 0;

  // Call here to insure that last_valid is not us. This places the line at the
  // beginning of the valid list.
  validate( line ); 
                    
  m_size += size;   // Update the size after adding the new line
  VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache allocated " << size
                  << " bytes (" << m_size << " / " << m_max_size << " used)" << "\n"; );

  // Grab the oldest CacheLine object
  CacheLineBase* local_last_valid = m_last_valid;

  while ( m_size > m_max_size ) {

    if ( local_last_valid == line || !local_last_valid ) {
      // De-allocated all lines except the current one which are not
      // held up currently by other threads.
      break;
    }

    // Deallocate the oldest CacheLine object if nothing is using it.
    bool invalidated = local_last_valid->try_invalidate();
    if (invalidated) { // If we were able to clear it...
      local_evictions++;
      local_last_valid = m_last_valid;  // Update the local pointer to the new oldest CacheLine.
    } else {
      // If we can't deallocate current line,
      // switch to the one used a bit more recently.
      local_last_valid = local_last_valid->m_prev;
    }
  }

  {
    Mutex::WriteLock cache_lock( m_stats_mutex );
    m_evictions += local_evictions; // Update class evictions stat
    
    // Warn about exceeding the cache size. Note that the warning is
    // printed only if the size now is a multiple of the previous size
    // at which the warning was printed, so it will warn say when the
    // cache size is 1.5^n GB. This will limit the number of warnings
    // to a representative subset.
    double factor = 1.5;
    double MB = 1024.0 * 1024.0; // 1 MB in bytes
    if ((m_size > m_max_size) && (m_size > factor*m_last_size)) {
      VW_OUT(WarningMessage, "cache")
        << "Cache size (" << m_size / MB
        << " MB) is larger than the requested maximum cache size (" << m_max_size / MB
        << " MB). Consider increasing --cache-size-mb for this program.\n";
      m_last_size = m_size;
    }
    
  }
}

void vw::Cache::resize( size_t size ) {
  // WARNING! YOU CAN NOT HOLD THE CACHE MUTEX AND THEN CALL
  // INVALIDATE. That's a line -> cache -> line mutex hold. A deadlock!
  size_t local_size, local_max_size;
  CacheLineBase* local_last_valid;
  { // Locally buffer variables that require Cache Mutex
    RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);
    m_max_size       = size;
    local_size       = m_size;
    local_max_size   = m_max_size;
    local_last_valid = m_last_valid;
  }
  // Keep deallocating objects until we shrink under the new size limit
  while ( local_size > local_max_size ) {
    VW_ASSERT( local_last_valid, LogicErr() << "Cache is empty but has nonzero size" );
    // Deallocate the last invalid CacheLine object
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

// Note that this call does not actually deallocate the data from the CacheLine object.
// It is up to the originating call to do that.  This call only removes all reference in 
// the Cache class to the CacheLine object.
void vw::Cache::deallocate( size_t size, CacheLineBase *line ) {
  RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);

  // This call implies the need to call invalidate (move to top of invalid list)
  invalidate( line );

  m_size -= size; // Remove the given size contribution.
  VW_CACHE_DEBUG( VW_OUT(DebugMessage, "cache") << "Cache deallocated " << size << " bytes (" << m_size << " / " << m_max_size << " used)" << "\n"; )
}

// TODO: Could we use some sort of linked list class to handle this stuff?

void vw::Cache::validate( CacheLineBase *line ) {
  RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);
  // If the input line is already most valid, done!
  if( line == m_first_valid ) 
    return;
  // This is the last line, we need to retreat the last valid pointer by one.
  if( line == m_last_valid ) 
    m_last_valid = line->m_prev;
  // If this is the first in the invalid list, we need to advance the first invalid pointer by one.
  if( line == m_first_invalid ) 
    m_first_invalid = line->m_next;
  // Adjust the elements before and after the input element to restore the linked list
  //  with the current element removed. TODO: Make this a function?
  if( line->m_next ) 
    line->m_next->m_prev = line->m_prev;
  if( line->m_prev ) 
    line->m_prev->m_next = line->m_next;
  // Make whatever is now first valid come after the input line
  line->m_next = m_first_valid;
  line->m_prev = 0; // The new line is first, nothing before it!
  
  // Update first valid pointer to point to the new object
  if( m_first_valid ) 
    m_first_valid->m_prev = line;
  m_first_valid = line;
  
  // Handle case where this is the first valid element to be validated
  if( ! m_last_valid ) 
    m_last_valid = line;
}

void vw::Cache::invalidate( CacheLineBase *line ) {
  RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);
  // Update first and last pointers if they point to the line
  if( line == m_first_valid ) m_first_valid = line->m_next;
  if( line == m_last_valid  ) m_last_valid  = line->m_prev;
  // Extract the line from its current location in the linked list
  if( line->m_next ) line->m_next->m_prev = line->m_prev;
  if( line->m_prev ) line->m_prev->m_next = line->m_next;
  // Set the line to the first place in the list
  line->m_next = m_first_invalid;
  line->m_prev = 0;
  if( m_first_invalid ) m_first_invalid->m_prev = line;
  m_first_invalid = line;
}

void vw::Cache::remove( CacheLineBase *line ) {
  RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);
  // Update list pointers if they pointed to the line
  if( line == m_first_valid   ) m_first_valid   = line->m_next;
  if( line == m_last_valid    ) m_last_valid    = line->m_prev;
  if( line == m_first_invalid ) m_first_invalid = line->m_next;
  // Extract the line from its current location in the linked list
  if( line->m_next ) line->m_next->m_prev = line->m_prev;
  if( line->m_prev ) line->m_prev->m_next = line->m_next;
  line->m_next = line->m_prev = 0;
}

void vw::Cache::deprioritize( CacheLineBase *line ) {
  RecursiveMutex::Lock cache_lock(m_line_mgmt_mutex);
  // Already the last item, done!
  if( line == m_last_valid  ) return;
  // Update the first valid pointer if needed
  if( line == m_first_valid ) m_first_valid = line->m_next;
  // Extract the line from its current location in the linked list
  if( line->m_next ) line->m_next->m_prev = line->m_prev;
  if( line->m_prev ) line->m_prev->m_next = line->m_next;
  // Set the line to the last place in the list
  line->m_prev = m_last_valid;
  line->m_next = 0;
  m_last_valid->m_next = line;
  m_last_valid = line;
}
