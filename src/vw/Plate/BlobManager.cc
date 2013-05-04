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


#include <vw/Plate/BlobManager.h>
#include <vw/Plate/Exception.h>
#include <vw/Core/Debugging.h>

// Vision Workbench
#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>

#include <boost/format.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

namespace fs = boost::filesystem;

namespace {
  static const vw::uint64 BLOB_MAX_SIZE = 1717986920; // 1.6 GB
  static const boost::format blob_tmpl("%s/plate_%u.blob");
}

#define WHEREAMI if(::vw::vw_log().is_enabled(VerboseDebugMessage, "platefile.blob")) ::vw::vw_out(VerboseDebugMessage, "platefile.blob")

namespace vw {
namespace platefile {

const uint64 BlobManager::BLOB_MIN_WRITABLE_SCORE = 2;

uint64 BlobManager::BlobKey::score() const {
  if (locked)          return 1;
  if (size > BLOB_MAX_SIZE) return 0;
  return size + BLOB_MIN_WRITABLE_SCORE;
}

bool BlobManager::BlobKey::can_write() const {
  return score() >= BLOB_MIN_WRITABLE_SCORE;
}

std::string BlobManager::name_from_id(uint32 blob_id) const {
  boost::format blob_name(blob_tmpl);
  return boost::str(blob_name % m_directory % blob_id);
}

uint32 BlobManager::num_blobs() const {
  Mutex::Lock lock(m_mutex);
  return m_blobs.size();
}

uint64 BlobManager::blob_size(uint32 blob_id) const {
  Mutex::Lock lock(m_mutex);
  const blob_by_id_t& lookup = m_blobs.get<0>();
  blob_by_id_t::const_iterator i = lookup.find(blob_id);
  VW_ASSERT(i != lookup.end(), ArgumentErr() << "No such blob id " << blob_id);
  return i->size;
}

uint32 BlobManager::request_lock() {
  WHEREAMI << std::endl;
  Mutex::Lock lock(m_mutex);
  blob_by_score_t& lookup = m_blobs.get<1>();

  blob_by_score_t::iterator i = lookup.begin();
  if (i == lookup.end() || !i->can_write())
    return locked_add_blob();

  VW_ASSERT(!i->locked, LogicErr() << "Tried to lock an already-locked blob");
  lookup.modify(i, BlobKey::SetLock(true));

  return i->id;
}

uint32 BlobManager::locked_add_blob() {
  WHEREAMI << std::endl;

  uint32 next_id = 0;
  const blob_by_id_t& lookup = m_blobs.get<0>();
  blob_by_id_t::const_reverse_iterator i = lookup.rbegin();

  if (i != lookup.rend())
    next_id = i->id + 1;

  std::pair<blob_tracker_t::iterator, bool> ret = m_blobs.insert(BlobKey(0, next_id, true));
  VW_ASSERT(ret.second, LogicErr() << "Failed to add blob " << next_id);
  return next_id;
}

void BlobManager::release_lock(uint32 blob_id) {
  WHEREAMI << "release " << blob_id << std::endl;
  Mutex::Lock lock(m_mutex);
  blob_by_id_t& lookup = m_blobs.get<0>();
  blob_by_id_t::iterator i = lookup.find(blob_id);
  VW_ASSERT(i != lookup.end(), ArgumentErr() << "No such blob id " << blob_id);
  VW_ASSERT(i->locked, LogicErr() << "Tried to unlock already-unlocked blob");

  std::string fn = name_from_id(blob_id);

  uint64 size = 0;
  if (fs::exists(fn))
    size = fs::file_size(fn);

  lookup.modify(i, BlobKey::SetUnlockSize(size));
}

BlobManager::BlobManager(const std::string& directory)
  : m_directory(directory)
{
  if (!fs::exists(directory))
    return;
  boost::regex re("plate_(\\d+)\\.blob");
  typedef fs::directory_iterator iter_t;

  BOOST_FOREACH(const fs::path& p, boost::make_iterator_range(iter_t(directory), iter_t())) {
    boost::cmatch matches;
    if (!boost::regex_match(p.filename().c_str(), matches, re))
      continue;

    std::string blob_id_str(matches[1].first, matches[1].second);
    uint32 blob_id = boost::lexical_cast<uint32>(blob_id_str);

    std::pair<blob_tracker_t::iterator, bool> ret = m_blobs.insert(BlobKey(fs::file_size(p), blob_id, false));
    VW_ASSERT(ret.second, LogicErr() << "Failed to add blob " << blob_id);
  }
}

}} // namespace vw::platefile
