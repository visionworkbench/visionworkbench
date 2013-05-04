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


#ifndef __VW_PLATE_BLOB_MANAGER__
#define __VW_PLATE_BLOB_MANAGER__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/member.hpp>

namespace vw {
namespace platefile {

  // The BlobManager keeps track of how much data has been stored in each blob
  // so far.  This allows it to pick blobs that have enough space to store a
  // new block of a given size.  The blob manager also controls the
  // locking/unlocking of blobs, and can load balance blobs writes by
  // alternating which blob is offered up for writing data.
  //
  // The BlobManager is thread safe.
  class BlobManager {

    static const uint64 BLOB_MIN_WRITABLE_SCORE;
    struct BlobKey {
      uint64 size;
      uint32 id;
      bool locked;

      BlobKey(uint64 size, uint32 id, bool locked) : size(size), id(id), locked(locked) {}

      uint64 score() const;
      bool can_write() const;

      struct SetLock {
        bool lock;
        SetLock(bool lock) : lock(lock) {}
        void operator()(BlobKey& k) {k.locked = lock;}
      };

      struct SetUnlockSize {
        uint64 size;
        SetUnlockSize(uint64 size) : size(size) {}
        void operator()(BlobKey& k) {k.size = size; k.locked = false;}
      };
    };

    // largest unlocked blob (ordered_non_unique<score>)
    // largest blob id       (ordered_unique<id>)
    // blob id key           (ordered_unique<id>)
    typedef boost::multi_index_container<BlobKey,
              boost::multi_index::indexed_by<boost::multi_index::ordered_unique<boost::multi_index::member<BlobKey, uint32, &BlobKey::id> >,
                                             boost::multi_index::ordered_non_unique<boost::multi_index::const_mem_fun<BlobKey, uint64, &BlobKey::score>, std::greater<uint64> > > >
            blob_tracker_t;

    typedef blob_tracker_t::nth_index<0>::type blob_by_id_t;
    typedef blob_tracker_t::nth_index<1>::type blob_by_score_t;

    mutable vw::Mutex m_mutex;
    std::string m_directory;
    blob_tracker_t m_blobs;

    uint32 locked_add_blob();

  public:

    // Return the number of tracked blobs
    uint32 num_blobs() const;

    // Returns the size of a blob in bytes
    uint64 blob_size(uint32 blob_id) const;

    // Create a new blob manager.
    BlobManager(const std::string& directory);

    // Request a blob to write to that has sufficient space. Returns the blob
    // index of a locked blob that you have sole access to write to.
    uint32 request_lock();

    // Given a blob id, return the filename of the corresponding blob
    std::string name_from_id(uint32 blob_id) const;

    // Release the blob lock
    void release_lock(uint32 blob_id);
  };


}} // namespace vw::platefile

#endif // __VW_PLATE_BLOB_MANAGER__
