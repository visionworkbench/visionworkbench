#include <vw/Plate/BlobManager.h>

// Vision Workbench
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>

// -------------------------------------------------------------------
//                            BLOB_MANAGER
// -------------------------------------------------------------------

void vw::platefile::BlobManager::BlobManager::next_blob_index() {
  // Move to the next blob
  m_blob_index++;
  if (m_blob_index >= m_blob_locks.size())
    m_blob_index = 0;
}

/// Create a new blob manager.  The max_blob_size is specified in
/// units of megabytes.
vw::platefile::BlobManager::BlobManager(int64 max_blob_size, int nblobs) : 
  m_max_blob_size(max_blob_size * 1024 * 1024), m_blob_index(0) { 
  
  m_blob_locks.resize(nblobs);
  
  for (int i=0; i < m_blob_locks.size(); ++i) {
    m_blob_locks[i] = false;
  }

}

/// Return the number of blobs currently in use.
int vw::platefile::BlobManager::num_blobs() {
  Mutex::Lock lock(m_mutex);
  return m_blob_locks.size();
}

vw::int64 vw::platefile::BlobManager::max_blob_size() { 
  Mutex::Lock lock(m_mutex);
  return m_max_blob_size; 
}

/// Request a blob to write to that has sufficient space to write at
/// least 'size' bytes.  Returns the blob index of a locked blob
/// that you have sole access to write to.
///
/// size is specified in bytes.
//
// TODO: This is pretty simple logic, and would not be very
// efficient because it blocks on write if it catches up to a blob
// that is still locked.  We should add real blob selection logic
// here at a later date.
int vw::platefile::BlobManager::request_lock(int64 size) {
  Mutex::Lock lock(m_mutex);

  // First, we check to see if the next blob is free.  If not, we
  // wait for a release event and recheck.
  while(m_blob_locks[m_blob_index] != false)
    m_blob_release_condition.wait(lock);
  
  // Then we lock it, increment the blob index, and return it.
  m_blob_locks[m_blob_index] = true;      
  int idx = m_blob_index;
  next_blob_index();
  return idx;
}

// Release the blob lock and update its write index (essentially
// "committing" the write to the blob when you are finished with it.).
int vw::platefile::BlobManager::release_lock(int blob_id) {
  Mutex::Lock lock(m_mutex);
  m_blob_locks[blob_id] = false;
  m_blob_release_condition.notify_all();
}

