// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Core/ProgressCallback.h>
#include <vw/Core/Log.h>
#include <vw/Plate/LocalIndex.h>
#include <vw/Plate/RemoteIndex.h>
using namespace vw;
using namespace vw::platefile;

#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
namespace fs = boost::filesystem;

// -------------------------------------------------------------------
//                            LOCAL INDEX
// -------------------------------------------------------------------


// ------------------------ Private Methods --------------------------

std::string vw::platefile::LocalIndex::index_filename() const {
  return m_plate_filename + "/plate.index";
}

std::string vw::platefile::LocalIndex::log_filename() const {
  return m_plate_filename + "/plate.log";
}

std::vector<std::string> vw::platefile::LocalIndex::blob_filenames() const {

  std::vector<std::string> result;

  if ( !fs::exists( m_plate_filename ) ) 
    vw_throw(IOErr() << "LocalIndex::blob_filenames() -- could not open platefile directory.");

  // Create a regular expression for matching the pattern 'plate_<number>.blob"
  boost::regex re;
  re.assign("plate_\\d+\\.blob", boost::regex_constants::icase);

  fs::directory_iterator end_itr; // default construction yields past-the-end

  // Iterate through the files in the platefile directory and return
  // any that match the regex above.
  for ( fs::directory_iterator itr( m_plate_filename ); itr != end_itr; ++itr ) {
    if (boost::regex_match(itr->leaf(), re))
      result.push_back(itr->leaf());
  }

  return result;
}

// Load index entries by iterating through TileHeaders saved in the
// blob file.  This function essentially rebuilds an index in memory
// using entries that had been previously saved to disk.
void vw::platefile::LocalIndex::load_index(std::string plate_filename, 
                                           std::vector<std::string> const& blob_files) {

  std::cout << "\tLoading index: " << plate_filename <<"\n";

  for (unsigned int i = 0; i < blob_files.size(); ++i) {
    // this->log() << "Loading index entries from blob file: " 
    //             << m_plate_filename << "/" << blob_files[i] << "\n";
    
    
    TerminalProgressCallback tpc(InfoMessage, "\t    " + blob_files[i] + " : ");
    tpc.report_progress(0);

    // Extract the current blob id as an integer.
    boost::regex re;
    re.assign("(plate_)(\\d+)(\\.blob)", boost::regex_constants::icase);
    boost::cmatch matches;
    boost::regex_match(blob_files[i].c_str(), matches, re);
    if (matches.size() != 4)
      vw_throw(IOErr() << "LocalIndex::load_index() -- could not parse blob number from blob filename.");
    std::string blob_id_str(matches[2].first, matches[2].second);
    int current_blob_id = atoi(blob_id_str.c_str());
    
    Blob blob(m_plate_filename + "/" + blob_files[i], true);
    Blob::iterator iter = blob.begin();
    while (iter != blob.end()) {
      TileHeader hdr = *iter;
      IndexRecord rec;
      rec.set_blob_id(current_blob_id);
      rec.set_blob_offset(iter.current_base_offset());
      rec.set_status(INDEX_RECORD_VALID);
      //      rec.set_tile_size(iter.current_data_size());  // ?? do we still need this?
      m_root->insert(rec, hdr.col(), hdr.row(), hdr.depth(), hdr.transaction_id());
      tpc.report_progress(float(iter.current_base_offset()) / blob.size());
      ++iter;
    }
    tpc.report_finished();
  }
}

void vw::platefile::LocalIndex::save_index_file() const {

  // First, check to make sure the platefile directory exists.
  if ( !fs::exists( m_plate_filename ) )
    fs::create_directory(m_plate_filename);
  std::ofstream ofstr( this->index_filename().c_str() );
  if (!ofstr.is_open())
    vw_throw(IOErr() << "LocalIndex::save_index_file() -- Could not create index file for writing.");

  // Serialize the header data.
  m_header.SerializeToOstream(&ofstr);
  ofstr.close();

}

// ------------------------ Public Methods --------------------------

/// Create a new index.  User supplies a pre-configure blob manager.
vw::platefile::LocalIndex::LocalIndex( std::string plate_filename, IndexHeader new_index_info) :
  m_plate_filename(plate_filename),
  m_blob_manager(boost::shared_ptr<BlobManager>( new BlobManager() )), 
  m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {

  // Start with the new_index_info, which provides the tile size, file
  // type, pixel and channel types, etc.  Then we augment it with a
  // few things to help track the new platefile.
  m_header = new_index_info;

  // Create a unique (random) platefile_id using the random()
  // function.  It would probably be better to use some sort of fancy
  // MD5 hash of the current time and date, or something like that,
  // but little 'ol random() will probably work just fine for our
  // humble purposes.
  srandom(time(0));
  m_header.set_platefile_id(vw::int32(random()));
  
  // Set up the IndexHeader and write it to disk.
  m_header.set_version(VW_PLATE_INDEX_VERSION);
  m_header.set_transaction_read_cursor(0);   // Transaction 0 is the empty mosaic
  m_header.set_transaction_write_cursor(1);  // Transaction 1 is the first valid transaction

  this->save_index_file();

  // Create the logging facility
  m_log = boost::shared_ptr<LogInstance>( new LogInstance(this->log_filename()) );
  this->log() << "Created new index: \"" << this->index_filename() << "\n"
              <<  m_header.DebugString() << "\n";
}

/// Open an existing index from a file on disk.
vw::platefile::LocalIndex::LocalIndex(std::string plate_filename) : 
  m_plate_filename(plate_filename),                                    
  m_blob_manager(boost::shared_ptr<BlobManager>( new BlobManager() )),
  m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {

  // First, check to make sure the platefile directory exists.
  if ( !fs::exists( plate_filename ) ) 
    vw_throw(IOErr() << "Could not open platefile: \"" << plate_filename << "\" for reading.");

  // The load the index file & parse the IndexHeader protobuffer
  std::ifstream ifstr( this->index_filename().c_str() );
  if (!ifstr.is_open()) 
    vw_throw(IOErr() << "LocalIndex::LocalIndex() -- Could not load index file for reading.");

  bool worked = m_header.ParseFromIstream(&ifstr);
  if (!worked)
    vw_throw(IOErr() << "LocalIndex::LocalIndex() -- Could not parse index information from " 
             << this->index_filename()); 

  ifstr.close();

  // Create the logging facility
  try {
    m_log = boost::shared_ptr<LogInstance>( new LogInstance(this->log_filename()) );
    this->log() << "Reopened index \"" << this->index_filename() << "\n"
                << m_header.DebugString() << "\n";
  } catch (IOErr &e) {
    // If we fail to open the log, then it's ok for now.  However, the
    // program _will_ crash if we try to do something with the index
    // that generates a log message.
    vw_out(WarningMessage, "plate") << "WARNING: could not open index log file. "
                                    << "Proceed with caution.";
  }

  // Load the actual index data
  std::vector<std::string> blob_files = this->blob_filenames();
  this->load_index(plate_filename, blob_files);
}

/// Use this to send data to the index's logfile like this:
///
///   index_instance.log() << "some text for the log...\n";
///
std::ostream& vw::platefile::LocalIndex::log () {
  if (!m_log)
    vw_throw(LogicErr() << "LocalIndex::log() - attempted to write to log when the log wasn\'t open.");
  return (*m_log)(InfoMessage, "console");
}


/// Attempt to access a tile in the index.  Throws an
/// TileNotFoundErr if the tile cannot be found.
IndexRecord vw::platefile::LocalIndex::read_request(int col, int row, int depth, int transaction_id) {
  Mutex::Lock lock(m_mutex);
  IndexRecord rec;
  try {
    rec = m_root->search(col, row, depth, transaction_id);
  } catch (IndexErr &e) {
    vw_throw(TileNotFoundErr() << "Invalid index.  Tiles do not exist at the given location.\n");
  }
  return rec;
}
  
// Writing, pt. 1: Locks a blob and returns the blob id that can
// be used to write a tile.
int vw::platefile::LocalIndex::write_request(int size) {  
  return m_blob_manager->request_lock(size); 
}

// Writing, pt. 2: Supply information to update the index and
// unlock the blob id.
void vw::platefile::LocalIndex::write_complete(TileHeader const& header, IndexRecord const& record) {
  m_blob_manager->release_lock(record.blob_id()); 

  Mutex::Lock lock(m_mutex);
  m_root->insert(record, header.col(), header.row(), header.depth(), header.transaction_id());
  m_root->invalidate_records(header.col(), header.row(), header.depth());
}


// Clients are expected to make a transaction request whenever
// they start a self-contained chunk of mosaicking work.  .
int32 vw::platefile::LocalIndex::transaction_request(std::string transaction_description) {
  Mutex::Lock lock(m_mutex);

  // Pick the next transaction ID, increment the cursor, and then save
  // the cursor to a file.  (Saving to a file gurantees that we won't
  // accidentally assign two transactions the same ID if the index
  // server crashes and has to be restarted.
  int32 transaction_id = m_header.transaction_write_cursor();
  m_header.set_transaction_write_cursor(m_header.transaction_write_cursor() + 1);
  this->save_index_file();
  this->log() << "Transaction " << transaction_id << " started: " << transaction_description << "\n";
  return transaction_id;
}

// Once a chunk of work is complete, clients can "commit" their
// work to the mosaic by issueing a transaction_complete method.
void vw::platefile::LocalIndex::transaction_complete(int32 transaction_id) {
  Mutex::Lock lock(m_mutex);

  this->log() << "Transaction " << transaction_id << " finished.\n";
  
  // If this ID is the next one after the current write ID, then we
  // increment the read ID.
  if (transaction_id == m_header.transaction_read_cursor() + 1) {
    
    // Update the transaction read cursor
    m_header.set_transaction_read_cursor(transaction_id);
    this->save_index_file();

    // Replay the remaining transaction IDs to see if the next ids are outstanding.
    bool match;
    do {
      match = false;
      google::protobuf::RepeatedField<const int32>::iterator iter = m_header.complete_transaction_ids().begin();
      while (iter != m_header.complete_transaction_ids().end()) {
        if (*iter == m_header.transaction_read_cursor() + 1) {
          m_header.set_transaction_read_cursor(*iter);
          this->save_index_file();
          match = true;
        }
        ++iter;
      }
    } while (match);

  } else {
    // This wasn't the next ID.  We store it in the queue and save it for later.
    m_header.mutable_complete_transaction_ids()->Add(transaction_id);
    this->save_index_file();

    const int MAX_OUTSTANDING_TRANSACTION_IDS = 10;
    if (m_header.complete_transaction_ids().size() > MAX_OUTSTANDING_TRANSACTION_IDS) 
      vw_out(0) << "WARNING: Detected potential stall in the transaction pipeline.  "
                << "There are " << m_header.complete_transaction_ids().size() 
                << " completed IDs waiting in the queue.\n";
  }
}

// Return the current location of the transaction cursor.  This
// will be the last transaction id that refers to a coherent
// version of the mosaic.
int32 vw::platefile::LocalIndex::transaction_cursor() {
  Mutex::Lock lock(m_mutex);
  return m_header.transaction_read_cursor();
}

