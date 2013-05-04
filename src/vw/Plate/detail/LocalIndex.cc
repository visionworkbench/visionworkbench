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


#include <vw/Core/ProgressCallback.h>
#include <vw/Core/Log.h>
#include <vw/FileIO/TemporaryFile.h>
#include <vw/Plate/detail/LocalIndex.h>
#include <vw/Plate/detail/RemoteIndex.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/BlobManager.h>
#include <vw/Plate/detail/Seed.h>

using namespace vw;
using namespace vw::platefile;
using namespace vw::platefile::detail;

#include <fstream>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <algorithm>
namespace fs = boost::filesystem;

LocalIndexPage::LocalIndexPage(std::string filename, uint32 level, uint32 base_col,
                               uint32 base_row, uint32 page_width, uint32 page_height)
  : IndexPage(level, base_col, base_row, page_width, page_height),
    m_filename(filename), m_needs_saving(false)
{
  if (fs::exists(filename)) {
    this->deserialize();
  }
}

LocalIndexPage::~LocalIndexPage() {
  this->sync();
}

// Hijack this method momentartarily to mark the page as "dirty" by
// setting m_needs_saving to true.
void LocalIndexPage::set(TileHeader const& header, IndexRecord const& record) {

  // First call up to the parent class and let the original code run.
  IndexPage::set(header, record);

  // Then, mark this page is 'dirty' so that it gets saved to disk
  // when destroyed.
  m_needs_saving = true;
}

void LocalIndexPage::sync() {
  if (m_needs_saving) {
    this->serialize();
    m_needs_saving = false;
  }
}

void LocalIndexPage::serialize() {

  fs::path pagename(m_filename);
  fs::path parent(pagename.parent_path());

  try {
    fs::create_directories(parent);
  } catch ( const std::exception& e ) {
    vw_throw(IOErr() << "Could not create IndexPage directory: " << e.what());
  }

  std::string tmpname;
  {
    TemporaryFile tmp(parent.string(), false);
    IndexPage::serialize(tmp);
    tmpname = tmp.filename();
  }

  if (::rename(tmpname.c_str(), m_filename.c_str()) == -1)
    vw_throw(IOErr() << "LocalIndexPage::serialize(): failed to rename temp page to " << m_filename << ": " << ::strerror(errno));

  m_needs_saving = false;
}

void LocalIndexPage::deserialize() {

  std::ifstream istr(m_filename.c_str());
  if (!istr.good())
    vw_throw(IOErr() << "IndexPage::deserialize() failed.  Could not open " << m_filename << " for reading.");

  // Call up to superclass to finish deserializing.
  try {
    IndexPage::deserialize(istr);
  } catch (const vw::IOErr &e) {
    // Add more useful error reporting.
    vw_throw(IOErr() << "Failed to load page \"" << m_filename << "\": " << e.what());
  }

  m_needs_saving = false;
}

 // ----------------------------------------------------------------------
 //                    LOCAL INDEX PAGE GENERATOR
 // ----------------------------------------------------------------------

 LocalPageGenerator::LocalPageGenerator( std::string filename, uint32 level, uint32 base_col,
                                         uint32 base_row, uint32 page_width, uint32 page_height)
 : m_filename( filename ), m_level(level), m_base_col(base_col), m_base_row(base_row),
   m_page_width(page_width), m_page_height(page_height) {}

boost::shared_ptr<IndexPage>
LocalPageGenerator::generate() const {
 return boost::shared_ptr<IndexPage>(new LocalIndexPage(m_filename, m_level,
       m_base_col, m_base_row,
       m_page_width, m_page_height) );
}

boost::shared_ptr<PageGeneratorBase>
LocalPageGeneratorFactory::create(uint32 level, uint32 base_col, uint32 base_row, uint32 page_width, uint32 page_height)
{
  // Generate a filename
  std::ostringstream filename;
  filename << m_plate_filename
    << "/index"
    << "/" << level
    << "/" << base_row
    << "/" << base_col;

  // Create the proper type of page generator.
  boost::shared_ptr<PageGeneratorBase> page_gen(
      new LocalPageGenerator(filename.str(), level, base_col, base_row,
        page_width, page_height) );

  return page_gen;
}

std::string LocalPageGeneratorFactory::who() const {
  return fs::path(m_plate_filename).stem().string();
}

// -------------------------------------------------------------------
//                            LOCAL INDEX
// -------------------------------------------------------------------

// ------------------------ Private Methods --------------------------

std::string LocalIndex::index_filename() const {
  return m_plate_filename + "/plate.index";
}

std::string LocalIndex::log_filename() const {
  return m_plate_filename + "/plate.log";
}

std::vector<std::string> LocalIndex::blob_filenames() const {

  std::vector<std::string> result;

  if ( !fs::exists( m_plate_filename ) )
    vw_throw(IOErr() << "Could not list blob filenames for plate (" << m_plate_filename << "): no such directory");

  boost::regex re("plate_(\\d+)\\.blob");
  typedef fs::directory_iterator iter_t;

  BOOST_FOREACH(const fs::path& elt, std::make_pair(iter_t(m_plate_filename), iter_t())) {
    boost::cmatch matches;
    if (!boost::regex_match(elt.filename().c_str(), matches, re)) {
      vw_out(DebugMessage, "plate") << "Skipping non-blob-filename " << elt.filename() << std::endl;
      continue;
    }

    std::string blob_id_str(matches[1].first, matches[1].second);
    uint32 blob_id = boost::lexical_cast<uint32>(blob_id_str);
    if (result.size() < blob_id+1)
      result.resize(blob_id+1);
    result[blob_id] = elt.string();
  }

  return result;
}

void LocalIndex::save_index_file() const {

  std::ofstream ofstr( this->index_filename().c_str() );
  if (!ofstr.is_open())
    vw_throw(IOErr() << "LocalIndex::save_index_file() -- Could not create index file for writing.");

  // Serialize the header data.
  m_header.SerializeToOstream(&ofstr);
  ofstr.close();

}

 // ------------------------ Public Methods --------------------------

 /// Create a new index.
 LocalIndex::LocalIndex( std::string plate_filename, IndexHeader new_index_info)
  : PagedIndex(boost::shared_ptr<PageGeneratorFactory>( new LocalPageGeneratorFactory(plate_filename) )),
   m_plate_filename(plate_filename), m_blob_manager(boost::shared_ptr<BlobManager>( new BlobManager(m_plate_filename) ))
{
   // Create subdirectory for storing hard copies of index pages.
   std::string base_index_path = plate_filename + "/index";

   if (fs::exists(base_index_path)) {
      open_impl();
#define WARN_IF_DIFFERENT(field) do {if (new_index_info.field() != m_header.field()) vw_out(WarningMessage, "plate") << "Refusing to change default for field " << #field << " from " << m_header.field() << " to " << new_index_info.field() << std::endl;} while(0)
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

   fs::create_directories(base_index_path);

   // Start with the new_index_info, which provides the tile size, file
   // type, pixel and channel types, etc.  Then we augment it with a
   // few things to help track the new platefile.
   m_header = new_index_info;

   // Check to make sure we've selected a sane file type
   if (m_header.channel_type() == VW_CHANNEL_FLOAT32) {
     if (m_header.tile_filetype() == "png" || m_header.tile_filetype() == "jpg")
       vw_out(WarningMessage, "plate")
         << "Constructing 32-bit floating point platefile with a non-32-bit file-type ("
         << m_header.tile_filetype() << ").\n";
   }

#define WARN_IF_SET(field) do {if (new_index_info.has_ ## field()) vw_out(ErrorMessage, "plate") << #field << " is a private field. Ignoring your value " << m_header.field() << std::endl;} while(0)
   WARN_IF_SET(platefile_id);
   WARN_IF_SET(version);
   WARN_IF_SET(transaction_read_cursor);
   WARN_IF_SET(transaction_write_cursor);
   WARN_IF_SET(num_levels);
#undef WARN_IF_SET

   // Create a unique (random) platefile_id using the random()
   // function.  It would probably be better to use some sort of fancy
   // MD5 hash of the current time and date, or something like that,
   // but little 'ol random() will probably work just fine for our
   // humble purposes.
   plate_seed_random();
   m_header.set_platefile_id(vw::uint32(random()));

   // Set up the IndexHeader and write it to disk.
   m_header.set_version(VW_PLATE_INDEX_VERSION);
   m_header.set_transaction_read_cursor(0);  // Transaction 0 is the empty mosaic
   m_header.set_transaction_write_cursor(1); // Transaction 1 is the first valid transaction to write into
   m_header.set_num_levels(0);               // Index initially contains zero levels
   this->save_index_file();

   // Create the logging facility
   m_log = boost::shared_ptr<LogInstance>( new LogInstance(this->log_filename()) );
   this->log() << "Created new index: \"" << this->index_filename() << "\n"
               <<  m_header.DebugString() << "\n";
 }

  void LocalIndex::open_impl() {
    // First, check to make sure the platefile directory exists.
    if ( !fs::exists( m_plate_filename ) )
      vw_throw(ArgumentErr() << "Could not open local platefile (" << m_plate_filename << "): no such directory");

    // The load the index file & parse the IndexHeader protobuffer
    std::ifstream ifstr( this->index_filename().c_str() );
    if (!ifstr.is_open())
      vw_throw(IOErr() << "LocalIndex::LocalIndex() -- Could not load index file for reading.");

    bool worked = m_header.ParseFromIstream(&ifstr);
    if (!worked)
      vw_throw(IOErr() << "LocalIndex::LocalIndex() -- Could not parse index information from "
          << this->index_filename());

    ifstr.close();

    // Load Index Levels for PagedIndex
    for (uint32 level = 0; level < this->num_levels(); ++level) {
      boost::shared_ptr<IndexLevel> new_level( new IndexLevel(m_page_gen_factory, level,
            m_page_width, m_page_height,
            m_default_cache_size) );
      m_levels.push_back(new_level);
    }

    m_log = boost::shared_ptr<LogInstance>( new LogInstance(this->log_filename()) );
    this->log() << "Reopened index \"" << this->index_filename() << "\n" << m_header.DebugString() << "\n";
  }

/// Open an existing index from a file on disk.
LocalIndex::LocalIndex(std::string plate_filename) :
  PagedIndex(boost::shared_ptr<PageGeneratorFactory>( new LocalPageGeneratorFactory(plate_filename) ) ),  // superclass constructor
  m_plate_filename(plate_filename), m_blob_manager(boost::shared_ptr<BlobManager>( new BlobManager(m_plate_filename) ))
{
  open_impl();
}

 /// Use this to send data to the index's logfile like this:
 ///
 ///   index_instance.log() << "some text for the log...\n";
 ///
 std::ostream& LocalIndex::log() {
   return (*m_log)(InfoMessage, "console");
 }

// Clients are expected to make a transaction request whenever
// they start a self-contained chunk of mosaicking work.  .
Transaction LocalIndex::transaction_request(std::string transaction_description,
                                            TransactionOrNeg transaction_id_override_neg) {

   if (transaction_id_override_neg != -1) {
     Transaction transaction_id_override = transaction_id_override_neg.promote();
     // If the user has chosen to override the transaction ID's, then we
     // use the transaction ID they specify.  We also increment the
     // transaction_write_cursor if necessary so that we dole out a
     // reasonable transaction ID if we are ever asked again without an
     // override.
     Transaction max_trans_id = std::max(m_header.transaction_write_cursor(), transaction_id_override+1);
     m_header.set_transaction_write_cursor(max_trans_id);
     this->save_index_file();
     this->log() << "Transaction " << transaction_id_override << " started: " << transaction_description << "\n";
     return transaction_id_override;

   } else {

     // Pick the next transaction ID, increment the cursor, and then save
     // the cursor to a file.  (Saving to a file gurantees that we won't
     // accidentally assign two transactions the same ID if the index
     // server crashes and has to be restarted.
     Transaction transaction_id = m_header.transaction_write_cursor();
     m_header.set_transaction_write_cursor(transaction_id + 1);
     this->save_index_file();
     this->log() << "Transaction " << transaction_id << " started: " << transaction_description << "\n";
     return transaction_id;
   }
 }

 // Once a chunk of work is complete, clients can "commit" their
 // work to the mosaic by issuing a transaction_complete method.
 void LocalIndex::transaction_complete(Transaction transaction_id, bool update_read_cursor) {

   // First we save (sync) the index pages to disk
   //this->sync();

   if ( update_read_cursor ) {
     Transaction max_trans_id = std::max(m_header.transaction_read_cursor(), transaction_id+0);
     m_header.set_transaction_read_cursor(max_trans_id);
     this->save_index_file();
   }

   this->log() << "Transaction " << transaction_id << " complete.  "
               << "[ read_cursor = " << m_header.transaction_read_cursor() << " ]\n";
 }

 // Once a chunk of work is complete, clients can "commit" their
 // work to the mosaic by issuing a transaction_complete method.
 void LocalIndex::transaction_failed(Transaction transaction_id) {

   // Log the failure.
   this->log() << "Transaction " << transaction_id << " FAILED.\n";

 }

 // Return the current location of the transaction cursor.  This
 // will be the last transaction id that refers to a coherent
 // version of the mosaic.
 Transaction LocalIndex::transaction_cursor() {
   return m_header.transaction_read_cursor();
 }

 // Load index entries by iterating through TileHeaders saved in the
 // blob file.  This function essentially rebuilds an index in memory
 // using entries that had been previously saved to disk.
 void LocalIndex::rebuild_index() {

   vw_out(InfoMessage) << "Rebuilding index: " << m_plate_filename <<"\n";

   // index in blob_filenames is the blobfile id
   std::vector<std::string> blob_files = this->blob_filenames();
   for (uint32 blob_id = 0; blob_id < blob_files.size(); ++blob_id) {
     const std::string& name = blob_files[blob_id];
     if (name.empty())
       continue;

     TerminalProgressCallback tpc("plate", "\t --> Rebuild from blob" + vw::stringify(blob_id) + " : ");
     tpc.report_progress(0);

     ReadBlob blob(name);

     IndexRecord rec;
     typedef Blob::iterator iter_t;
     for (iter_t tile = blob.begin(), end = blob.end(); tile != end; ++tile) {
       const BlobTileRecord& hdr = *tile;
       rec.set_blob_id(blob_id);
       rec.set_blob_offset(tile.current_base_offset());
       rec.set_filetype(hdr.hdr.filetype());
       this->write_update(hdr.hdr, rec);
       tpc.report_progress(float(tile.current_base_offset()) / float(blob.size()));
    }
    tpc.report_finished();
  }
}

// -----------------------    I/O      ----------------------

/// Writing, pt. 1: Reserve a blob lock
uint32 LocalIndex::write_request() {
  return m_blob_manager->request_lock();
}

/// Writing, pt. 1: Reserve a blob lock
void LocalIndex::write_update(TileHeader const& header, IndexRecord const& record) {

  // Store the number of tiles that are contained in the mosaic.
  size_t starting_size = m_levels.size();

  // Write the update to the PagedIndex superclass.
  PagedIndex::write_update(header, record);

  // If adding the record resulted in more levels, we save that
  // information to the index header.
  if (m_levels.size() != starting_size) {
    m_header.set_num_levels(boost::numeric_cast<uint32>(m_levels.size()));
    this->save_index_file();
  }
}

/// Writing, pt. 3: Signal the completion
void LocalIndex::write_complete(uint32 blob_id) {
  m_blob_manager->release_lock(blob_id);
}

