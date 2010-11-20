// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_BLOB_MANAGER__
#define __VW_PLATE_BLOB_MANAGER__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------
  //                            BLOB_MANAGER
  // -------------------------------------------------------------------

  /// The BlobManager keeps track of how much data has been stored in
  /// each blob so far.  This allows it to pick blobs that have enough
  /// space to store a new block of a given size.  The blob manager also
  /// controls the locking/unlocking of blobs, and can load balance
  /// blobs writes by alternating which blob is offered up for writing
  /// data.
  ///
  /// The BlobManager is thread safe.
  class BlobManager {

    struct BlobCacheRecord {
      bool locked;
      uint64 current_blob_offset;
      //      time_t lock_time;

      // Note: with the end_of_file_ptr written at the beginning of
      // the blob file, the initial blob_offset should be
      // 3 * 8 bytes = 24 bytes.
      BlobCacheRecord() : locked(false), current_blob_offset(0) {}

      void lock() {
        this->locked = true;
        //        this->lock_time = time(0);
      }

      // void unlock_if_timeout() {
      //   time_t now = time(0);
      //   if (this->locked && now - lock_time > 60) { // 60 second timeout
      //     this->locked = false;
      //   }
      // }

      void unlock(uint64 current_blob_offset) {
        this->current_blob_offset = current_blob_offset;
        this->locked = false;
      }
    };

    vw::uint64 m_max_blob_size;
    unsigned int m_max_blobs;
    std::vector<BlobCacheRecord> m_blob_locks;
    int m_blob_index;
    vw::Mutex m_mutex;

    // A method to poll for an available blob.  Returns -1 if there
    // are no blobs available.
    int get_next_available_blob();

    // Helper function for incrementing blob ids, and wrapping around the end.
    void increment_blob_index(int &blob_index);

  public:

    /// Create a new blob manager.  The max_blob_size is specified in
    /// units of megabytes.
    BlobManager(vw::uint64 max_blob_size = 1792, int initial_nblobs = 1, unsigned max_blobs = 16384);

    /// Return the number of blobs currently in use.
    unsigned num_blobs();

    /// Return the max size of a blob (in megabytes)
    uint64 max_blob_size();

    /// Request a blob to write to that has sufficient space to write at
    /// least 'size' bytes.  Returns the blob index of a locked blob
    /// that you have sole access to write to.
    ///
    /// size is specified in bytes.
    //
    // TODO: This is pretty simple logic so far, and would not be very
    // efficient because it blocks on write if it catches up to a blob
    // that is still locked.  We should add real blob selection logic
    // here at a later date.
    int request_lock(vw::uint64 &size);

    // Release the blob lock and update its write index (essentially
    // "committing" the write to the blob when you are finished with
    // it.).  The blob_offset is used to roughly gauge how full the
    // blob is.  Once the blob_offset moves past max_blob_size defined
    // above, this blob is considered "full" and is no longer offered
    // when locks are requested.
    void release_lock(int blob_id, uint64 blob_offset);
  };


}} // namespace vw::platefile

#endif // __VW_PLATE_BLOB_MANAGER__
