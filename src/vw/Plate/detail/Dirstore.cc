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


#include <vw/Plate/detail/Dirstore.h>
#include <vw/Plate/Exception.h>
#include <vw/FileIO/TemporaryFile.h>
#include <vw/Plate/detail/Seed.h>
#include <boost/foreach.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/format.hpp>
#include <boost/interprocess/sync/file_lock.hpp>


#include <sys/types.h>
#include <dirent.h>
#include <limits.h>
#include <stddef.h>
#include <unistd.h>
#include <set>

namespace fs  = boost::filesystem;
namespace ipc = boost::interprocess;
namespace sys = boost::system;

namespace {
  using namespace vw;
  using namespace vw::platefile;
  typedef fs::directory_iterator iter_t;

  class DirWriteState : public WriteState {
    public:
      DirWriteState(const Transaction& id) : id(id) {}
      virtual std::string what() const {
        return std::string("DirWriteState[tid = " + stringify(id) + "]");
      }
      Transaction id;
  };

  const uint32 BUCKET_SIZE = 32;
  uint64 calc_bucket(uint32 row, uint32 col) {
    return
      uint32(row / BUCKET_SIZE) * BUCKET_SIZE +
      uint32(col / BUCKET_SIZE);
  }

  void noop() {}
  void cleanup(const fs::path& path) {
    fs::remove(path);
  }
  struct raii {
    typedef boost::function<void (void)> FuncT;
    FuncT m_leave;
    raii(FuncT enter, FuncT leave, bool skip = false) : m_leave(skip ? &noop : leave) { if (!skip) enter(); }
    ~raii() VW_NOTHROW {
      try {
        m_leave();
      } catch (const std::exception& e) {
        vw_out(ErrorMessage) << "Raii destructor for " << this << " tried to throw exception: " << e.what() << std::endl;
      } catch (...) {
        vw_out(ErrorMessage) << "Raii destructor for " << this << " tried to throw unknown exception." << std::endl;
      }
    }
    void disarm() VW_NOTHROW {
      m_leave = &noop;
    }
  };

  TileHeader make_hdr(uint32 level, uint32 row, uint32 col, Transaction id, const std::string& filetype) {
    TileHeader hdr;
    hdr.set_col(col);
    hdr.set_row(row);
    hdr.set_level(level);
    hdr.set_transaction_id(id);
    hdr.set_filetype(filetype);
    return hdr;
  }

  bool HasData(const Tile& t) {
    return t.data;
  }

# define FMT_TID    "%010u"
# define FMT_LEVEL  "%03u"
# define FMT_ROWCOL "%08u"
# define FMT_BUCKET "%07u"

  std::string format_tid(const Transaction& id) {
    static const boost::format fmt_(FMT_TID);
    boost::format fmt(fmt_);
    return boost::str(fmt % id);
  }
  std::string format_rowcol(const uint32& rowcol) {
    static const boost::format fmt_(FMT_ROWCOL);
    boost::format fmt(fmt_);
    return boost::str(fmt % rowcol);
  }

  void sort_and_limit_size(Datastore::TileSearch& tiles, uint32 limit) {
    Datastore::TileSearch::iterator end = tiles.end();
    if (limit > 0 && tiles.size() > limit)
      end = tiles.begin() + limit;

    std::partial_sort(tiles.begin(), end, tiles.end(), OrderTileByTidDesc());
    if (end != tiles.end())
      tiles.resize(limit);
  }

  bool check_dir(const fs::path& p) {
    boost::system::error_code err;
    iter_t i(p, err);
    if (err == boost::system::errc::no_such_file_or_directory)
      return false;
    else if (err)
      vw_throw(IOErr() << "Dirstore: failed to open directory " << p << ": " << err.message());
    return true;
  }
#if 0
  struct remove_t {
    typedef bool result_type;
    bool operator()(const fs::path& a) {
      std::cerr << "remove " << a << std::endl;
      return fs::remove(a);
    }
  };

  struct link_t {
    typedef void result_type;
    void operator()(const fs::path& a, const fs::path& b) {
      std::cerr << "link " << b << " -> " << a << std::endl;
      fs::create_hard_link(a, b);
    }
  };

  struct rename_t {
    typedef void result_type;
    void operator()(const fs::path& a, const fs::path& b) {
      std::cerr << "rename " << a << " -> " << b << std::endl;
      fs::rename(a, b);
    }
  };
#endif

}

#if 0
  BUCKET_SIZE = 32 (32x32 slot buckets)
  bucket = int(row / BUCKET_SIZE) * BUCKET_SIZE + int(col / BUCKET_SIZE)

  4 possible cases
    1) Single-ID, Single-Location
      -> map(/by-tid/tid/level/bucket/row/col/+) reduce()
    2) Multiple-ID, Single-Location
      -> map(/by-loc/level/bucket/row/col/+) reduce(match-ids)
    3) Single-ID, Multiple-Location
      -> map(/by-tid/tid/level/bucket/+) reduce(match-loc)
    4) Multiple-ID, Multiple-Loc
      -> map(/by-loc/level/bucket/+) reduce(match-ids)

Tile add procedure
  tid_dir = /by-tid/<tid>/<level>/<bucket>/<row>/<col>/
  loc_dir = /by-loc/<level>/<bucket>/<row>/<col>/<tid>/
  bak_dir   = /tmpdir/XXXXXX

  open,write,close bak_dir/new

  mkdir $tid_dir, ignore EEXIST
  mkdir $loc_dir, ignore EEXIST

  boost::file_lock(tid_dir/lock)

  old_filetype = get_filetype($tid_dir)

  rename $tid_dir/old_filetype -> $bak_dir/tid [sentinel: reverse rename]
  rename $loc_dir/old_filetype -> $bak_dir/loc [sentinel: reverse rename]

  hardlink($bak_dir/new, $tid_dir/$new_filetype); [sentinel: rm $tid_dir/$new_filetype]
  hardlink($bak_dir/new, $loc_dir/$new_filetype); [sentinel: rm $loc_dir/$new_filetype]

  disarm_sentinels

  ~file_lock()
  rm -rf bak_dir
#endif

namespace vw { namespace platefile { namespace detail {


std::string Dirstore::path_by_tid(const platefile::Transaction& id, uint32 level, uint32 row, uint32 col, const std::string& filetype) const {
  // path, tid, level, bucket, row, col, filetype
  static const boost::format fmt_("%s/by-tid/" FMT_TID "/" FMT_LEVEL "/" FMT_BUCKET "/" FMT_ROWCOL "/" FMT_ROWCOL "/%s");
  boost::format fmt(fmt_);
  return boost::str(fmt % m_plate_path % id % level % calc_bucket(row, col) % row % col % filetype);
}

std::string Dirstore::path_by_loc(const platefile::Transaction& id, uint32 level, uint32 row, uint32 col, const std::string& filetype) const {
  // path, level, bucket, row, col, tid, filetype
  static const boost::format fmt_("%s/by-loc/" FMT_LEVEL "/" FMT_BUCKET "/" FMT_ROWCOL "/" FMT_ROWCOL "/" FMT_TID "/%s");
  boost::format fmt(fmt_);
  return boost::str(fmt % m_plate_path % level % calc_bucket(row, col) % row % col % id % filetype);
}

std::string Dirstore::path_by_loc_no_tid(uint32 level, uint32 row, uint32 col) const {
  // path, level, bucket, row, col
  static const boost::format fmt_("%s/by-loc/" FMT_LEVEL "/" FMT_BUCKET "/" FMT_ROWCOL "/" FMT_ROWCOL);
  boost::format fmt(fmt_);
  return boost::str(fmt % m_plate_path % level % calc_bucket(row, col) % row % col);
}

std::string Dirstore::path_by_bucket(uint32 level, uint64 bucket) const {
  // path, level, bucket, row, col
  static const boost::format fmt_("%s/by-loc/" FMT_LEVEL "/" FMT_BUCKET);
  boost::format fmt(fmt_);
  return boost::str(fmt % m_plate_path % level % bucket);
}

std::string Dirstore::tid_dir(const Transaction& id) const {
  static const boost::format fmt_("%s/by-tid/" FMT_TID);
  boost::format fmt(fmt_);
  return boost::str(fmt % m_plate_path % id);
}

# undef FMT_TID
# undef FMT_LEVEL
# undef FMT_ROWCOL
# undef FMT_BUCKET

void Dirstore::save_index_file() const {
  fs::path header_path = m_plate_path + "/header";
  std::ofstream f( header_path.string().c_str() );
  VW_ASSERT(f.is_open(), IOErr() << "Dirstore::save_index_file() : Could not create index file for writing.");
  bool worked = m_hdr.SerializeToOstream(&f);
  VW_ASSERT(worked, IOErr() << "Dirstore::save_index_file(): Failed to write index header");
  f.close();
}

void Dirstore::init() {
  fs::path header_path = m_plate_path + "/header";

  VW_ASSERT(fs::is_directory(m_plate_path), ArgumentErr() << "Plate " << m_plate_path << " is not a directory.");
  VW_ASSERT(fs::exists(header_path),        ArgumentErr() << "Index file " << header_path << " must exist! (This datastore does not support remote data)");

  std::ifstream f( header_path.string().c_str() );
  VW_ASSERT(f.is_open(), IOErr() << "Dirstore: Could not load index file for reading.");
  bool worked = m_hdr.ParseFromIstream(&f);
  VW_ASSERT(worked, IOErr() << "Dirstore: Could not parse index header " << header_path);
  f.close();
}

Dirstore::Dirstore(const Url& u)
  : m_plate_path(u.path())
{
  VW_ASSERT(u.scheme() == "dir", ArgumentErr() << "Dirstore must be local");
  init();
}

Dirstore::Dirstore(const Url& u, const IndexHeader& d)
  : m_plate_path(u.path())
{
  VW_ASSERT(u.scheme() == "dir", ArgumentErr() << "Dirstore must be local");

  if (fs::exists(m_plate_path)) {
    init();
#define WARN_IF_DIFFERENT(field) do {if (d.field() != m_hdr.field()) vw_out(WarningMessage, "plate") << "Refusing to change default for field " << #field << " from " << m_hdr.field() << " to " << d.field() << std::endl;} while(0)
    // Ignore attempts to change the header, but warn if anything is different.
    WARN_IF_DIFFERENT(tile_size);
    WARN_IF_DIFFERENT(tile_filetype);
    WARN_IF_DIFFERENT(pixel_format);
    WARN_IF_DIFFERENT(channel_type);
    WARN_IF_DIFFERENT(type);
    WARN_IF_DIFFERENT(description);
#undef WARN_IF_DIFFERENT
    return;
  }

  fs::create_directories(m_plate_path + "/tmp");

  m_hdr = d;

#define WARN_IF_SET(field) do {if (m_hdr.has_ ## field()) vw_out(ErrorMessage, "plate") << #field << " is a private field. Ignoring your value " << m_hdr.field() << std::endl;} while(0)
  WARN_IF_SET(platefile_id);
  WARN_IF_SET(version);
  WARN_IF_SET(transaction_read_cursor);
  WARN_IF_SET(transaction_write_cursor);
  WARN_IF_SET(num_levels);
#undef WARN_IF_SET

  plate_seed_random();
  m_hdr.set_platefile_id(vw::uint32(random()));

  // Set up the IndexHeader and write it to disk.
  m_hdr.set_version(0);
  m_hdr.set_transaction_read_cursor(0);  // Transaction 0 is the empty mosaic
  m_hdr.set_transaction_write_cursor(1); // Transaction 1 is the first valid transaction to write into
  m_hdr.set_num_levels(0);               // Index initially contains zero levels
  this->save_index_file();
  init();
}

Transaction Dirstore::transaction_begin(const std::string& description, TransactionOrNeg override) {
  boost::optional<Transaction> id;
  if (override == -1) {
    id = m_hdr.transaction_write_cursor();
    m_hdr.set_transaction_write_cursor(id.get() + 1);
    this->save_index_file();
  } else {
    id = override.promote();
    if (id.get() >= m_hdr.transaction_write_cursor()) {
      m_hdr.set_transaction_write_cursor(id.get()+1);
      this->save_index_file();
    }
  }
  fs::create_directories(tid_dir(id.get()));

  this->audit_log()() << "Transaction " << id.get() << " started: " << description << "\n";
  return id.get();
}

void Dirstore::transaction_end(Transaction id, bool update_read_cursor) {
   if ( update_read_cursor ) {
     if (id > m_hdr.transaction_read_cursor()) {
       m_hdr.set_transaction_read_cursor(id);
       this->save_index_file();
      }
   }

   if (!fs::is_directory(tid_dir(id)))
     this->error_log()() << "Transcation " << id << " was never started!" << std::endl;

   this->audit_log()() << "Transaction " << id << " complete.  "
                       << "[ read_cursor = " << m_hdr.transaction_read_cursor() << " ]\n";
}

// returns "" if the col directory doesn't exist
std::string Dirstore::get_filetype(uint32 level, uint32 row, uint32 col, const Transaction& id) const {
  const std::string dir_name(path_by_tid(id, level, row, col, ""));

  if (!check_dir(dir_name))
    return "";

  std::string filetype;
  BOOST_FOREACH(const fs::path& p, std::make_pair(iter_t(dir_name), iter_t())) {
    const std::string fn = p.filename().string();
    // if we've set filetype already, we already saw a file in this dir
    VW_ASSERT(filetype.empty(), LogicErr() << "Directory " << dir_name << " has multiple files (at least " << filetype << " and " << fn << ")");
    filetype = fn;
  }
  return filetype;
}

struct PathCompareByLeaf {
  bool operator()(const fs::path& a, const fs::path& b) {
    return a.filename() < b.filename();
  }
};

void Dirstore::add_one_id(Datastore::TileSearch& tiles, uint32 level, uint32 row, uint32 col, const Transaction& id) const {
  std::string type(get_filetype(level, row, col, id));
  if (!type.empty())
    tiles.push_back(make_hdr(level, row, col, id, type));
}

void Dirstore::add_top_id(Datastore::TileSearch& tiles, uint32 level, uint32 row, uint32 col) const {
  const std::string dir_name(path_by_loc_no_tid(level, row, col));
  const std::string high = format_tid(m_hdr.transaction_read_cursor());

  if (!check_dir(dir_name))
    return;

  fs::path top;
  BOOST_FOREACH(const fs::path& p, std::make_pair(iter_t(dir_name), iter_t())) {
    if (p > top && p <= high)
      top = p;
  }

  if (top.empty())
    return;

  uint32 tid = boost::lexical_cast<uint32>(top.filename());
  std::string type = get_filetype(level, row, col, tid);
  VW_ASSERT(!type.empty(), LogicErr() << "Dirstore: matching tile exists in " << dir_name << "but not in by-tid in " << path_by_tid(tid, level, row, col, ""));
  tiles.push_back(make_hdr(level, row, col, tid, type));
}

void Dirstore::add_ids_in_range(Datastore::TileSearch& tiles, uint32 level, uint32 row, uint32 col, const TransactionRange& r, uint32 limit) const {
  const std::string dir_name(path_by_loc_no_tid(level, row, col));
  const std::string low  = format_tid(r.first().promote());
  const std::string high = format_tid(r.last().newest() ? Transaction::MAX_POSSIBLE() : uint32(r.last().promote()));

  if (!check_dir(dir_name))
    return;

  Datastore::TileSearch new_tiles;
  BOOST_FOREACH(const fs::path& p, std::make_pair(iter_t(dir_name), iter_t())) {
    const std::string tid_s = p.filename().string();
    if (tid_s < low || tid_s > high)
      continue;
    uint32 tid = boost::lexical_cast<uint32>(tid_s);

    std::string type = get_filetype(level, row, col, tid);
    VW_ASSERT(!type.empty(), LogicErr() << "Dirstore: matching tile exists in " << dir_name << "but not in by-tid in " << path_by_tid(tid, level, row, col, ""));
    new_tiles.push_back(make_hdr(level, row, col, tid, type));
  }
  sort_and_limit_size(new_tiles, limit);
  tiles.insert(tiles.end(), new_tiles.begin(), new_tiles.end());
}


Datastore::TileSearch& Dirstore::one_id_one_loc(Datastore::TileSearch& tiles, uint32 level, uint32 row, uint32 col, const Transaction& id) const {
  add_one_id(tiles, level, row, col, id);
  return tiles;
}
Datastore::TileSearch& Dirstore::top_id_one_loc(Datastore::TileSearch& tiles, uint32 level, uint32 row, uint32 col) const {
  add_top_id(tiles, level, row, col);
  return tiles;
}

Datastore::TileSearch& Dirstore::many_id_one_loc(Datastore::TileSearch& tiles, uint32 level, uint32 row, uint32 col, const TransactionRange& r, uint32 limit) const {
  add_ids_in_range(tiles, level, row, col, r, limit);
  return tiles;
}

namespace {
  std::vector<uint64> bucketize_region(const BBox2u& region) {
    uint32 col = (region.min().x() / BUCKET_SIZE) * BUCKET_SIZE;
    uint32 row = (region.min().y() / BUCKET_SIZE) * BUCKET_SIZE;
    const uint32& end_col   = region.max().x(),
                  end_row   = region.max().y();

    std::vector<uint64> buckets;
    for (; col < end_col; col += BUCKET_SIZE)
      for (; row < end_row; row += BUCKET_SIZE)
        buckets.push_back(calc_bucket(row, col));
    return buckets;
  }
}

void Dirstore::bucket_iterate(uint32 level, const BBox2u& region, add_func_t add_func) const {

  const std::string row_min = format_rowcol(region.min().y()),
                    row_max = format_rowcol(region.max().y()),
                    col_min = format_rowcol(region.min().x()),
                    col_max = format_rowcol(region.max().x());

  std::vector<uint64> buckets = bucketize_region(region);
  BOOST_FOREACH(const uint64& bucket, buckets) {
    std::string bucket_path(path_by_bucket(level, bucket));
    if (!check_dir(bucket_path))
      continue;
    BOOST_FOREACH(const fs::path& p, std::make_pair(iter_t(bucket_path), iter_t())) {
      const std::string row_s = p.filename().string();
      if (row_s < row_min || row_s >= row_max)
        continue;
      if (!check_dir(p))
        continue;
      BOOST_FOREACH(const fs::path& p, std::make_pair(iter_t(p), iter_t())) {
        const std::string col_s = p.filename().string();
        if (col_s < col_min || col_s >= col_max)
          continue;
        uint32 row = boost::lexical_cast<uint32>(row_s);
        uint32 col = boost::lexical_cast<uint32>(col_s);
        add_func(row, col);
      }
    }
  }
}

Datastore::TileSearch& Dirstore::one_id_many_loc(Datastore::TileSearch& tiles, uint32 level, const BBox2u& region, const Transaction& id) const {
  bucket_iterate(level, region, boost::bind(&Dirstore::add_one_id, this, boost::ref(tiles), level, _1, _2, id));
  return tiles;
}

Datastore::TileSearch& Dirstore::top_id_many_loc(Datastore::TileSearch& tiles, uint32 level, const BBox2u& region) {
  bucket_iterate(level, region, boost::bind(&Dirstore::add_top_id, this, boost::ref(tiles), level, _1, _2));
  return tiles;
}

Datastore::TileSearch& Dirstore::many_id_many_loc(Datastore::TileSearch& tiles, uint32 level, const BBox2u& region, const TransactionRange& r, uint32 limit) const {
  bucket_iterate(level, region, boost::bind(&Dirstore::add_ids_in_range, this, boost::ref(tiles), level, _1, _2, r, limit));
  return tiles;
}

Datastore::TileSearch& Dirstore::head(Datastore::TileSearch& tiles, uint32 level, uint32 row, uint32 col, TransactionRange range, uint32 limit) {
  bool one_id  = (range.first() == range.last());
  bool top_id  = (one_id && range.first().newest());

  tiles.clear();

  if (top_id)
    return top_id_one_loc(tiles, level, row, col);
  else if (one_id)
    return one_id_one_loc(tiles, level, row, col, range.first().promote());
  else
    return many_id_one_loc(tiles, level, row, col, range, limit);
}

Datastore::TileSearch& Dirstore::head(Datastore::TileSearch& tiles, uint32 level, const BBox2u& region, TransactionRange range, uint32 limit) {
  bool one_id  = (range.first() == range.last());
  bool top_id  = (one_id && range.first().newest());
  bool one_loc = (region.width() == 1 && region.height() == 1);

  const Vector<uint32, 2>& origin = region.min();
  const uint32& row = origin.y();
  const uint32& col = origin.x();

  tiles.clear();

  if (top_id) {
    if (one_loc)
      return top_id_one_loc(tiles, level, row, col);
    else
      return top_id_many_loc(tiles, level, region);
  } else if (one_id) {
    if (one_loc)
      return one_id_one_loc(tiles, level, row, col, range.first().promote());
    else
      return one_id_many_loc(tiles, level, region, range.first().promote());
  } else {
    if (one_loc)
      return many_id_one_loc(tiles, level, row, col, range, limit);
    else
      return many_id_many_loc(tiles, level, region, range, limit);
  }
}

namespace {
  using namespace vw;
  using namespace vw::platefile;

  TileData slurp(const std::string& filename) {
    TileData data;

    std::ifstream f(filename.c_str(), std::ios::binary);
    VW_ASSERT(f.is_open(), IOErr() << "Failed to open tile file: " << filename);

    f.seekg(0, std::ios::end);
    size_t size = f.tellg();
    f.seekg(0, std::ios::beg);
    VW_ASSERT(!f.fail(), IOErr() << "Failed to find file size: " << filename);
    VW_ASSERT(size > 0,  IOErr() << "Zero-length file");

    data.reset(new std::vector<uint8>(size));
    f.read(reinterpret_cast<char*>(&data->operator[](0)), size);
    VW_ASSERT(!f.fail(), IOErr() << "Failed to read tile file: " << filename);
    return data;
  }

  void spit(const std::string& filename, const uint8* data, size_t size) {
    std::ofstream f(filename.c_str(), std::ios::binary|std::ios::out|std::ios::trunc);
    VW_ASSERT(f.is_open(), IOErr() << "Failed to open tile file: " << filename);
    f.write(reinterpret_cast<const char*>(data), size);
    if (f.fail()) {
      f.close();
      fs::remove(filename);
      vw_throw(IOErr() << "Failed to write tile file: " << filename);
    }
  }
}

Datastore::TileSearch& Dirstore::populate(TileSearch& tiles) {
  bool prune = false;
  BOOST_FOREACH(Tile& t, tiles) {
    const std::string filename = path_by_tid(t.hdr.transaction_id(), t.hdr.level(), t.hdr.row(), t.hdr.col(), t.hdr.filetype());
    t.data = slurp(filename.c_str());
    if (!t.data)
      prune = true;
  }
  if (prune) {
    TileSearch::iterator i = std::partition(tiles.begin(), tiles.end(), HasData);
    tiles.erase(i, tiles.end());
  }
  return tiles;
}

WriteState* Dirstore::write_request(const Transaction& id) {
  return new DirWriteState(id);
}

std::string Dirstore::get_lockfile(const Transaction& id) const {
  std::string lockfile = tid_dir(id) + "/lock";
  if (!fs::exists(lockfile))
    std::ofstream f(lockfile.c_str());
  return lockfile;
}

void Dirstore::write_update(WriteState& state_, uint32 level, uint32 row, uint32 col, const std::string& filetype, const uint8* data, uint64 size) {
  VW_ASSERT(filetype != "auto", ArgumentErr() << "write_update(): unsupported filetype 'auto'");
  VW_ASSERT(!filetype.empty(),  ArgumentErr() << "write_update(): unsupported empty filetype");

  DirWriteState* state = dynamic_cast<DirWriteState*>(&state_);
  VW_ASSERT(state, LogicErr() << "Cannot pass write states between different implementations!");

  const fs::path tid_dir  = path_by_tid(state->id, level, row, col, "");
  const fs::path loc_dir  = path_by_loc(state->id, level, row, col, "");

  TemporaryDir tmp_dir(m_plate_path + "/tmp", true, "update");
  const std::string tmp_tile = tmp_dir.filename() + "/new";

  spit(tmp_tile, data, size);
  fs::create_directories(tid_dir);
  fs::create_directories(loc_dir);

  const fs::path old_filetype = get_filetype(level, row, col, state->id);
  const fs::path old_tile_tid = tid_dir / old_filetype;
  const fs::path old_tile_loc = loc_dir / old_filetype;
  const fs::path old_tile_bak = fs::path(tmp_dir.filename()) / "old";

  const std::string lockfile = get_lockfile(state->id);
  ipc::file_lock lock(lockfile.c_str());

  bool (*remove)(const fs::path&)                  = &fs::remove;
  void   (*link)(const fs::path&, const fs::path&) = &fs::create_hard_link;
  void (*rename)(const fs::path&, const fs::path&) = &fs::rename;

  raii loc_remove(boost::bind(remove, old_tile_loc),
                  boost::bind(  link, old_tile_tid, old_tile_loc),
                  old_filetype.empty());

  raii tid_save(boost::bind(rename, old_tile_tid, old_tile_bak),
                boost::bind(rename, old_tile_bak, old_tile_tid),
                old_filetype.empty());

  const fs::path new_tile_tid = tid_dir / filetype,
                 new_tile_loc = loc_dir / filetype;

  raii tid_create(boost::bind(  link, tmp_tile, new_tile_tid),
                  boost::bind(remove, new_tile_tid));
  raii loc_create(boost::bind(  link, tmp_tile, new_tile_loc),
                  boost::bind(remove, new_tile_loc));

  if (level >= num_levels()) {
    m_hdr.set_num_levels(level+1);
    save_index_file();
  }

  loc_create.disarm();
  tid_create.disarm();
  tid_save.disarm();
  loc_remove.disarm();
}

void Dirstore::write_complete(WriteState& state_) {
  DirWriteState* state = dynamic_cast<DirWriteState*>(&state_);
  VW_ASSERT(state, LogicErr() << "Cannot pass write states between different implementations!");
  state->id = 0;
}

void Dirstore::flush() { }

IndexHeader Dirstore::index_header() const {
  return m_hdr;
}

vw::uint32 Dirstore::num_levels() const            { return m_hdr.num_levels();   }
vw::uint32 Dirstore::id() const                    { return m_hdr.platefile_id(); }
vw::uint32 Dirstore::tile_size() const             { return m_hdr.tile_size();    }
std::string Dirstore::tile_filetype() const        { return m_hdr.tile_filetype();}
vw::PixelFormatEnum Dirstore::pixel_format() const { return static_cast<PixelFormatEnum>(m_hdr.pixel_format()); }
vw::ChannelTypeEnum Dirstore::channel_type() const { return static_cast<ChannelTypeEnum>(m_hdr.channel_type()); }

Datastore::Logger Dirstore::audit_log() const {
  return boost::bind(&vw_out, InfoMessage, "console");
}

Datastore::Logger Dirstore::error_log() const {
  return boost::bind(&vw_out, ErrorMessage, "console");
}


}}} // vw::platefile::detail
