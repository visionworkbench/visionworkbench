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
  : m_index(Index::construct_open(u)), m_read_cache(DEFAULT_BLOB_CACHE_SIZE) {init();}

Blobstore::Blobstore(const Url& u, const IndexHeader& d)
  : m_index(Index::construct_create(u, d)), m_read_cache(DEFAULT_BLOB_CACHE_SIZE) {init();}

boost::shared_ptr<ReadBlob> Blobstore::open_read_blob(uint32 blob_id) {
  Mutex::Lock lock(m_handle_lock);
  if (m_handles.count(blob_id) == 0)
    m_handles.insert(std::make_pair(blob_id, m_read_cache.insert(BlobOpener(*m_index, blob_id))));
  return m_handles[blob_id];
}

boost::shared_ptr<Blob> Blobstore::open_write_blob(uint32 blob_id) {
  Mutex::Lock lock(m_handle_lock);
  boost::format blob_name(blob_tmpl);
  // TODO: Should I invalidate a read blob if you write to it? I don't think it's necessary...
  return boost::shared_ptr<Blob>(new Blob(str(blob_name % m_index->platefile_name() % blob_id)));
}

Transaction Blobstore::transaction_begin(const std::string& description, TransactionOrNeg override) {
  return m_index->transaction_request(description, override);
}
void Blobstore::transaction_end(Transaction transaction_id, bool update_read_cursor) {
  m_index->transaction_complete(transaction_id, update_read_cursor);
}

Datastore::meta_range Blobstore::head(uint32 level, uint32 row, uint32 col, TransactionRange range, uint32 limit) {
  typedef std::vector<TileHeader> vec_t;
  boost::shared_ptr<vec_t> tiles;

  vw_out(VerboseDebugMessage, "datastore") << "Blobstore::head("
    << "level=" << level
    << ", row=" << row
    << ", col=" << col
    << ", range=" << range
    << ", limit=" << limit << ")" << std::endl;

  {
    std::list<TileHeader> hdrs = m_index->search_by_location(col, row, level, range.first(), range.last());
    size_t len = hdrs.size();
    if (limit > 0 && len > limit) {
      hdrs.resize(limit);
      len = limit;
    }

    tiles.reset(new vec_t(len));
    std::copy(hdrs.begin(), hdrs.end(), tiles->begin());
  }
  return make_const_shared_range(tiles);
}

Datastore::meta_range Blobstore::head(uint32 level,   const BBox2u& region, TransactionRange range, uint32 limit) {
  typedef std::vector<TileHeader> vec_t;
  boost::shared_ptr<vec_t> tiles;

  vw_out(VerboseDebugMessage, "datastore") << "Blobstore::head("
    << "level=" << level
    << ", region=" << region
    << ", range=" << range
    << ", limit=" << limit << ")" << std::endl;

  {
    std::list<TileHeader> hdrs = m_index->search_by_region(level, region, range.first(), range.last(), 0);
    size_t len = hdrs.size();
    if (limit > 0 && len > limit) {
      hdrs.resize(limit);
      len = limit;
    }

    tiles.reset(new vec_t(len));
    std::copy(hdrs.begin(), hdrs.end(), tiles->begin());
  }
  return make_const_shared_range(tiles);
}

struct SortByPage {
  const Index& idx;
  SortByPage(const Index& idx) : idx(idx) {}
  bool operator()(const TileHeader& a, const TileHeader& b) {
    return idx.page_id(a.col(), a.row(), a.level()) < idx.page_id(b.col(), b.row(), b.level());
  }
};

struct SortByIndexRecord {
  typedef std::pair<IndexRecord, TileHeader> pair_t;
  bool operator()(const pair_t& a, const pair_t& b)
  {
    if (a.first.blob_id() == b.first.blob_id())
      return a.first.blob_offset() < b.first.blob_offset();
    return a.first.blob_id() < b.first.blob_id();
  }
};

Datastore::tile_range Blobstore::populate(TileHeader* hdrs_, size_t len) {
  // first, sort by page to keep page accesses together
  std::sort(hdrs_, hdrs_+len, SortByPage(*m_index));

  std::vector<std::pair<IndexRecord, TileHeader> > recs;
  recs.reserve(len);

  //if (hdr.filetype() != rec.filetype())
  //  vw_out(WarningMessage) << "input TileHeader doesn't match IndexRecord [filetype] [" << hdr.filetype() << " vs " << rec.filetype() << "]\n";

  BOOST_FOREACH(const TileHeader& hdr, boost::make_iterator_range(hdrs_, hdrs_+len)) {
    try {
      IndexRecord rec = m_index->read_request(hdr.col(), hdr.row(), hdr.level(), hdr.transaction_id(), true);
      recs.push_back(std::make_pair(rec, hdr));
    } catch (const TileNotFoundErr& e) {
      // I don't think this is possible if the hdrs are coming from get(). If
      // they're not coming from get(), the user created their own TileHeader
      // array incorrectly. In either case, this is a canary of something
      // really bad being wrong.
      vw_throw(LogicErr() << "Blobstore::populate(): cannot populate nonexistent tile: " << hdr);
    }
  }

  // Now we have index records... sort them by blob and then by offset, so we
  // keep blob reads together and in order. We need to keep the headers sorted
  // in the same order, too
  std::sort(recs.begin(), recs.end(), SortByIndexRecord());

  typedef std::vector<Tile> vec_t;
  boost::shared_ptr<vec_t> tiles(new vec_t());
  tiles->reserve(len);

  for (size_t i = 0; i < len; ++i) {
    const IndexRecord& rec = recs[i].first;
    const TileHeader&  hdr = recs[i].second;

    try {
      boost::shared_ptr<ReadBlob> blob = open_read_blob(rec.blob_id());
      BlobTileRecord tile_rec = blob->read_record(rec.blob_offset());
      if (tile_rec.hdr != hdr) {
        vw_out(ErrorMessage) << "output TileHeader doesn't match IndexRecord. skipping.\n";
        continue;
      }
      Tile tile;
      tile.hdr  = tile_rec.hdr;
      tile.data = tile_rec.data;
      tiles->push_back(tile);
    } catch (const IOErr& e) {
      // These are bad, and might indicate corruption, but probably shouldn't kill everything.
      error_log()() << "IOErr while reading tile " << hdr << ": " << e.what() << std::endl;
    }
  }

  return make_const_shared_range(tiles);
}

WriteState* Blobstore::write_request(const Transaction& id) {
  uint64 last_size;
  std::auto_ptr<BlobWriteState> state(new BlobWriteState(id));

  state->blob_id = m_index->write_request(last_size);
  state->blob = open_write_blob(state->blob_id);

  if (last_size != 0 && last_size != state->blob->size()) {
    error_log()() << "last close size did not match current size when opening "
         << state->blob->filename()
         << "  ( " << last_size << " != " << state->blob->size() << " )\n";
  }

  vw_out(DebugMessage, "blob") << "Opened blob " << state->blob_id << " ( size = " << state->blob->size() << " )\n";
  return state.release();
}

void Blobstore::write_update(WriteState& state_, uint32 level, uint32 row, uint32 col, const std::string& filetype, const uint8* data, uint64 size) {
  if (filetype == "auto")
    vw_throw(NoImplErr() << "write_update() does not support filetype 'auto'");

  BlobWriteState* state = dynamic_cast<BlobWriteState*>(&state_);
  VW_ASSERT(state, LogicErr() << "Cannot pass write states between different implementations!");

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

  // Close the blob for writing.
  state->blob.reset();

  // For debugging:
  vw_out(DebugMessage, "blob") << "Closed blob " << state->blob_id << " ( size = " << new_blob_size << " )\n";

  // Release the blob lock.
  m_index->write_complete(state->blob_id, new_blob_size);

  state->blob_id = 0;
}

void Blobstore::flush() {
  m_index->sync();
}

IndexHeader Blobstore::index_header() const {
  return m_index->index_header();
}

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

