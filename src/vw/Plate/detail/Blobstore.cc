#include <vw/Plate/detail/Blobstore.h>
#include <vw/Plate/detail/Index.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/Exception.h>
#include <boost/foreach.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/format.hpp>

namespace fs = boost::filesystem;

namespace {
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

void Blobstore::init() {
  const std::string& name = m_index->platefile_name();

  VW_ASSERT(fs::exists(name),       ArgumentErr() << "Plate directory " << name << " must exist. (This datastore does not support remote data)");
  VW_ASSERT(fs::is_directory(name), LogicErr() << "Plate " << name << " is not a directory.");

  m_error_log.add(vw_out(ErrorMessage, "console"));
  m_error_log.add(audit_log());
  m_read_blob_id = std::numeric_limits<uint32>::max();
}

Blobstore::Blobstore(const Url& u)
  : m_index(Index::construct_open(u)) {init();}

Blobstore::Blobstore(const Url& u, const IndexHeader& d)
  : m_index(Index::construct_create(u, d)) {init();}

boost::shared_ptr<Blob> Blobstore::open_blob(uint32 blob_id, bool readonly) {
  boost::format blob_name(blob_tmpl);

  if (readonly) {
    if (blob_id != m_read_blob_id) {
      m_read_blob.reset(new Blob(str(blob_name % m_index->platefile_name() % blob_id), readonly));
      m_read_blob_id = blob_id;
    }
    return m_read_blob;
  } else {
    return boost::shared_ptr<Blob>(new Blob(str(blob_name % m_index->platefile_name() % blob_id), readonly));
  }
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
    std::list<TileHeader> hdrs = m_index->search_by_location(col, row, level, range.first(), range.last(), false);
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
    std::list<TileHeader> hdrs = m_index->search_by_region(level, region, range.first(), range.last(), 0, false);
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

Datastore::tile_range Blobstore::populate(const TileHeader* hdrs, size_t len) {
  typedef std::vector<Tile> vec_t;
  boost::shared_ptr<vec_t> tiles(new vec_t());

  BOOST_FOREACH(const TileHeader& hdr, boost::make_iterator_range(hdrs, hdrs+len)) {
    Tile tile;

    try {
      IndexRecord rec = m_index->read_request(hdr.col(), hdr.row(), hdr.level(), hdr.transaction_id(), true);

      if (hdr.filetype() != rec.filetype())
        vw_out(WarningMessage) << "input TileHeader doesn't match IndexRecord [filetype] [" << hdr.filetype() << " vs " << rec.filetype() << "]\n";

      boost::shared_ptr<Blob> blob = open_blob(rec.blob_id(), true);

      // This read_header is unnecessary (we should alredy have it from hdr
      // but this check makes sure the plate is consistent
      tile.hdr  = blob->read_header(rec.blob_offset());

      if (tile.hdr != hdr)
        vw_out(ErrorMessage) << "output TileHeader doesn't match IndexRecord\n";

      tile.data = blob->read_data(rec.blob_offset());
      tiles->push_back(tile);
    } catch (const TileNotFoundErr& e) {
      // I don't think this is possible if the hdrs are coming from get(). If
      // they're not coming from get(), the user created their own TileHeader
      // array incorrectly. In either case, this is a canary of something
      // really bad being wrong.
      vw_throw(LogicErr() << "Blobstore::populate(): cannot populate nonexistent tile: " << hdr);
    } catch (const IOErr& e) {
      // These are bad, and might indicate corruption, but probably shouldn't kill everything.
      error_log() << "IOErr while reading tile " << hdr << ": " << e.what() << std::endl;
    }
  }
  return make_const_shared_range(tiles);
}

WriteState* Blobstore::write_request(const Transaction& id) {
  uint64 last_size;
  std::auto_ptr<BlobWriteState> state(new BlobWriteState(id));

  state->blob_id = m_index->write_request(last_size);
  state->blob = open_blob(state->blob_id, false);

  if (last_size != 0 && last_size != state->blob->size()) {
    error_log() << "last close size did not match current size when opening "
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
std::ostream& Blobstore::audit_log() {
  return m_index->log();
}
std::ostream& Blobstore::error_log() {
  return m_error_log;
}


}}} // vw::platefile::detail

