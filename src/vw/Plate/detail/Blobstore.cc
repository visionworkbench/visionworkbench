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


#include <vw/Plate/detail/Blobstore.h>
#include <vw/Plate/detail/Index.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/Exception.h>
#include <boost/foreach.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>

namespace fs = boost::filesystem;
namespace io = boost::iostreams;

namespace {
  static const size_t DEFAULT_BLOB_CACHE_SIZE = 8;
  static const boost::format blob_tmpl("%s/plate_%u.blob");

  class BlobWriteState : public vw::platefile::WriteState {
    public:
      BlobWriteState(const vw::platefile::Transaction& id) : transaction(id) {}
      virtual std::string what() const {
        return std::string("BlobWriteState[id = " + vw::stringify(blob_id) + "]");
      }
      vw::platefile::Transaction transaction;
      boost::shared_ptr<vw::platefile::Blob> blob;
      vw::uint32 blob_id;
  };

}

namespace vw { namespace platefile { namespace detail {

class BlobOpener {
    std::string fn;
  public:
    typedef ReadBlob value_type;
    BlobOpener(const Index& idx, uint32 blob_id) {
      boost::format blob_name(blob_tmpl);
      fn = str(blob_name % idx.platefile_name() % blob_id);
    }
    size_t size() const {return 1;}
    boost::shared_ptr<value_type> generate() const {
      return boost::shared_ptr<value_type>(new ReadBlob(fn));
    }
};

void Blobstore::init() {
  const std::string& name = m_index->platefile_name();

  VW_ASSERT(fs::exists(name),       ArgumentErr() << "Plate directory " << name << " must exist. (This datastore does not support remote data)");
  VW_ASSERT(fs::is_directory(name), LogicErr() << "Plate " << name << " is not a directory.");
}

Blobstore::Blobstore(const Url& u)
  : m_index(Index::construct_open(u)) {init();}

Blobstore::Blobstore(const Url& u, const IndexHeader& d)
  : m_index(Index::construct_create(u, d)) {init();}

boost::shared_ptr<ReadBlob> Blobstore::open_read_blob(uint32 blob_id) {
  Mutex::Lock lock(m_mutex);

  write_cache_t::const_iterator i = m_write_cache.find(blob_id);
  if (i != m_write_cache.end())
    return i->second;

  read_cache_by_filename_t& fn_cache = m_read_cache.get<1>();

  boost::format blob_name(blob_tmpl);
  const std::string fn = boost::str(blob_name % m_index->platefile_name() % blob_id);

  read_cache_by_filename_t::const_iterator j = fn_cache.find(fn);
  if (j != fn_cache.end()) {
    read_cache_by_age_t::iterator k = m_read_cache.project<0>(j);
    m_read_cache.relocate(m_read_cache.begin(), k);
    return *k;
  } else {
    boost::shared_ptr<ReadBlob> blob(new ReadBlob(fn));
    std::pair<read_cache_t::iterator, bool> i = m_read_cache.push_front(blob);
    VW_ASSERT(i.second, LogicErr() << "We just checked for the key " << fn << " and it's gone now!");
    if (m_read_cache.size() > DEFAULT_BLOB_CACHE_SIZE)
      m_read_cache.pop_back();
    return *i.first;
  }
}

boost::shared_ptr<Blob> Blobstore::open_write_blob(uint32 blob_id) {
  Mutex::Lock lock(m_mutex);

  VW_ASSERT(m_write_cache.count(blob_id) == 0, LogicErr() << "Cannot open a blob for writing more than once [tried " << blob_id << "]");

  boost::format blob_name(blob_tmpl);
  const std::string fn = boost::str(blob_name % m_index->platefile_name() % blob_id);
  boost::shared_ptr<Blob>& blob = m_write_cache[blob_id];
  blob.reset(new Blob(fn));

  // Expire the blob from the read cache, since we just opened it for write. It
  // will be opened from the write cache next time.
  // TODO: We should repopulate the cache entry, not just remove it
  read_cache_by_filename_t& fn_cache = m_read_cache.get<1>();
  fn_cache.erase(fn);
  return blob;
}

Transaction Blobstore::transaction_begin(const std::string& description, TransactionOrNeg override) {
  return m_index->transaction_request(description, override);
}
void Blobstore::transaction_end(Transaction transaction_id, bool update_read_cursor) {
  m_index->transaction_complete(transaction_id, update_read_cursor);
}

Datastore::TileSearch& Blobstore::head(TileSearch& buf, uint32 level, uint32 row, uint32 col, TransactionRange range, uint32 limit) {
  vw_out(VerboseDebugMessage, "datastore") << "Blobstore::head("
    << "level=" << level
    << ", row=" << row
    << ", col=" << col
    << ", range=" << range
    << ", limit=" << limit << ")" << std::endl;

  buf.clear();
  {
    std::list<TileHeader> hdrs = m_index->search_by_location(col, row, level, range.first(), range.last());
    size_t len = hdrs.size();
    if (limit > 0 && len > limit) {
      hdrs.resize(limit);
      len = limit;
    }

    buf.reserve(len);
    std::copy(hdrs.begin(), hdrs.end(), std::back_inserter(buf));
  }
  return buf;
}

Datastore::TileSearch& Blobstore::head(TileSearch& buf, uint32 level,   const BBox2u& region, TransactionRange range, uint32 limit) {
  vw_out(VerboseDebugMessage, "datastore") << "Blobstore::head("
    << "level=" << level
    << ", region=" << region
    << ", range=" << range
    << ", limit=" << limit << ")" << std::endl;

  buf.clear();
  {
    std::list<TileHeader> hdrs = m_index->search_by_region(level, region, range.first(), range.last());
    size_t len = hdrs.size();
    if (limit > 0 && len > limit) {
      hdrs.resize(limit);
      len = limit;
    }

    buf.reserve(len);
    std::copy(hdrs.begin(), hdrs.end(), std::back_inserter(buf));
  }
  return buf;
}

struct SortByPage {
  const Index& idx;
  SortByPage(const Index& idx) : idx(idx) {}
  bool operator()(const Tile& a, const Tile& b) {
    return idx.page_id(a.hdr.col(), a.hdr.row(), a.hdr.level()) < idx.page_id(b.hdr.col(), b.hdr.row(), b.hdr.level());
  }
};

struct SortByIndexRecord {
  bool operator()(const IndexRecord& a, const IndexRecord& b)
  {
    if (a.blob_id() == b.blob_id())
      return a.blob_offset() < b.blob_offset();
    return a.blob_id() < b.blob_id();
  }
};

bool HasData(const Tile& t) {
  return t.data;
}

Datastore::TileSearch& Blobstore::populate(TileSearch& hdrs) {
  // first, sort by page to keep page accesses together
  std::sort(hdrs.begin(), hdrs.end(), SortByPage(*m_index));

  typedef std::map<IndexRecord, Tile*, SortByIndexRecord> map_t;
  map_t recs;

  BOOST_FOREACH(Tile& tile, hdrs) {
    try {
      IndexRecord rec = m_index->read_request(tile.hdr.col(), tile.hdr.row(), tile.hdr.level(), tile.hdr.transaction_id(), true);
      recs[rec] = &tile;
    } catch (const TileNotFoundErr& e) {
      // I don't think this is possible if the hdrs are coming from get(). If
      // they're not coming from get(), the user created their own TileHeader
      // array incorrectly. In either case, this is a canary of something
      // really bad being wrong.
      vw_throw(LogicErr() << "Blobstore::populate(): cannot populate nonexistent tile: " << tile.hdr);
    }
  }
  VW_ASSERT(recs.size() == hdrs.size(), LogicErr() << "TileHeaders and IndexRecords should be 1:1");

  // keys in a std::map are sorted in ascending order according to the
  // comparison function.  SortByIndexRecord sorts by by blob and then by
  // offset, so we keep blob reads together and in order. We need to keep the
  // headers sorted in the same order, too

  bool prune = false;
  BOOST_FOREACH(map_t::value_type& tile, recs) {
    const IndexRecord& rec = tile.first;
    const TileHeader&  hdr = tile.second->hdr;
    Tile& t = *tile.second;

    try {
      boost::shared_ptr<ReadBlob> blob = open_read_blob(rec.blob_id());
      BlobTileRecord tile_rec = blob->read_record(rec.blob_offset());
      if (tile_rec.hdr.col() != hdr.col()
          || tile_rec.hdr.row() != hdr.row()
          || tile_rec.hdr.level() != hdr.level()
          || tile_rec.hdr.transaction_id() != hdr.transaction_id()
          || (hdr.has_filetype() && hdr.filetype().size() && tile_rec.hdr.filetype() != hdr.filetype())) {
        vw_out(ErrorMessage) << "output TileHeader doesn't match IndexRecord. skipping. [" << tile_rec.hdr << "] vs [" << hdr << "]\n";
        continue;
      }
      // Must copy the tilerec one, because we won't necessarily have a filetype in the search
      t.hdr  = tile_rec.hdr;
      t.data = tile_rec.data;
    } catch (const BlobIoErr& e) {
      // These are bad, and might indicate corruption, but probably shouldn't kill everything.
      error_log()() << "BlobIoErr while reading tile " << hdr << ": " << e.what() << std::endl;
      prune = true;
    } catch (const IOErr& e) {
      // These are bad, and might indicate corruption, but probably shouldn't kill everything.
      error_log()() << "IOErr while reading tile " << hdr << ": " << e.what() << std::endl;
      prune = true;
    }
  }

  if (prune) {
    TileSearch::iterator i = std::partition(hdrs.begin(), hdrs.end(), HasData);
    hdrs.erase(i, hdrs.end());
  }

  return hdrs;
}

WriteState* Blobstore::write_request(const Transaction& id) {
  std::auto_ptr<BlobWriteState> state(new BlobWriteState(id));

  state->blob_id = m_index->write_request();
  state->blob    = open_write_blob(state->blob_id);

  vw_out(DebugMessage, "blob") << "Opened blob " << state->blob_id << " ( size = " << state->blob->size() << " )\n";
  return state.release();
}

void Blobstore::write_update(WriteState& state_, uint32 level, uint32 row, uint32 col, const std::string& filetype, const uint8* data, uint64 size) {
  if (filetype == "auto")
    vw_throw(NoImplErr() << "write_update() does not support filetype 'auto'");

  BlobWriteState* state = dynamic_cast<BlobWriteState*>(&state_);
  VW_ASSERT(state, LogicErr() << "Cannot pass write states between different implementations!");

  vw_out(DebugMessage, "blob") << "Attempting to write " << size << " to blob "
                               << state->blob_id << " ( size = " << state->blob->size() << " )\n";

  // Check to see if we are exceeding our 32 bit pointer. All of the
  // Index Code and Blob Store code actually passes around a 64bit
  // pointer. It's only the Blob I/O that uses 32 bit pointers on
  // seek. We could make that 64bit, but instead lets write another
  // blob file. This way we have precision to know when we are going
  // to overflow instead of relying on break down.
  if ( size + state->blob->size() >= (uint64(0x1) << 32) ) {
    write_complete( state_ );
    // We don't call write_request as that would require us modifing a
    // pointer (WriteState) which is owned by a shared pointer in
    // PlateFile.
    state->blob_id = m_index->write_request();
    state->blob    = open_write_blob(state->blob_id);
    vw_out(DebugMessage, "blob") << "Switched over to blob " << state->blob_id << " ( size = " << state->blob->size() << " )\n";
  }

  TileHeader header;
  header.set_level(level);
  header.set_row(row);
  header.set_col(col);
  header.set_transaction_id(state->transaction);
  header.set_filetype(filetype);

  // 1. Write the data into the blob
  uint64 blob_offset = state->blob->write(header, data, size);

  // 2. Update the index
  IndexRecord write_record;
  write_record.set_blob_id(state->blob_id);
  write_record.set_blob_offset(blob_offset);
  write_record.set_filetype(header.filetype());

  m_index->write_update(header, write_record);
}

void Blobstore::write_complete(WriteState& state_) {
  BlobWriteState* state = dynamic_cast<BlobWriteState*>(&state_);
  VW_ASSERT(state, LogicErr() << "Cannot pass write states between different implementations!");

  // Fetch the size from the blob.
  uint64 new_blob_size = state->blob->size();

  {
    Mutex::Lock lock(m_mutex);

    // The blob might still technically be open for writing (in the read
    // cache), so flush it and drop our reference to it.
    state->blob->flush();
    state->blob.reset();
  }

  // For debugging:
  vw_out(DebugMessage, "blob") << "Closed blob " << state->blob_id << " ( size = " << new_blob_size << " )\n";

  // Release the blob lock.
  m_index->write_complete(state->blob_id);

  // Release cache
  m_write_cache.erase(state->blob_id);

  state->blob_id = 0;
}

void Blobstore::flush() {
  m_index->sync();
}

IndexHeader Blobstore::index_header() const     { return m_index->index_header();  }
uint32 Blobstore::num_levels() const            { return m_index->num_levels();    }
uint32 Blobstore::id() const                    { return m_index->platefile_id();  }
uint32 Blobstore::tile_size() const             { return m_index->tile_size();     }
std::string Blobstore::tile_filetype() const    { return m_index->tile_filetype(); }
PixelFormatEnum Blobstore::pixel_format() const { return m_index->pixel_format();  }
ChannelTypeEnum Blobstore::channel_type() const { return m_index->channel_type();  }

// LOGGING
Datastore::Logger Blobstore::audit_log() const {
  return boost::bind(&Index::log, boost::ref(m_index));
}

namespace {
  struct LogHelper {
    Index& i;
    typedef io::tee_device<std::ostream, std::ostream> tee_t;
    typedef io::stream<tee_t> log_t;
    boost::shared_ptr<tee_t> tee;
    boost::shared_ptr<log_t> log;
    LogHelper(Index& i) : i(i) {}
    std::ostream& operator()() {
      tee.reset(new tee_t(i.log(), vw_out(ErrorMessage, "console")));
      log.reset(new log_t(*tee));
      return *log;
    }
  };
}

Datastore::Logger Blobstore::error_log() const {
  return LogHelper(*m_index);
}


}}} // vw::platefile::detail

