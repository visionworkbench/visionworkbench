// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_BLOB_MANAGER__
#define __VW_PLATE_BLOB_MANAGER__

#include <string>
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
    vw::int64 m_max_blob_size;
    std::vector<bool> m_blob_locks;
    int m_blob_index;
    vw::Mutex m_mutex;
    vw::Condition m_blob_release_condition;

    void next_blob_index();

  public:

    /// Create a new blob manager.  The max_blob_size is specified in
    /// units of megabytes.
    BlobManager(vw::int64 max_blob_size = 2048, int nblobs = 2);

    /// Return the number of blobs currently in use.
    int num_blobs();

    vw::int64 max_blob_size();

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
    int request_lock(vw::int64 size);

    // Release the blob lock and update its write index (essentially
    // "committing" the write to the blob when you are finished with it.).
    int release_lock(int blob_id);
  };


}} // namespace vw::platefile

#endif // __VW_PLATE_BLOB_MANAGER__
