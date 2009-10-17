// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Core/Log.h>
#include <vw/Plate/Index.h>
using namespace vw;
using namespace vw::platefile;

#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
namespace fs = boost::filesystem;

#define VW_CURRENT_PLATEFILE_VERSION 1

// -------------------------------------------------------------------
//                            PLATE INDEX
// -------------------------------------------------------------------



// ------------------------ Private Methods --------------------------

std::string vw::platefile::Index::index_filename() const {
  return m_plate_filename + "/plate.index";
}

std::string vw::platefile::Index::log_filename() const {
  return m_plate_filename + "/plate.log";
}

std::vector<std::string> vw::platefile::Index::blob_filenames() const {

  std::vector<std::string> result;

  if ( !fs::exists( m_plate_filename ) ) 
    vw_throw(IOErr() << "Index::blob_filenames() -- could not open platefile directory.");

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
void vw::platefile::Index::load_index(std::vector<std::string> const& blob_files) {

  for (unsigned int i = 0; i < blob_files.size(); ++i) {
    this->log() << "\t--> Loading index entries from blob file: " 
                << m_plate_filename << "/" << blob_files[i] << "\n";
    
    // Extract the current blob id as an integer.
    boost::regex re;
    re.assign("(plate_)(\\d+)(\\.blob)", boost::regex_constants::icase);
    boost::cmatch matches;
    boost::regex_match(blob_files[i].c_str(), matches, re);
    if (matches.size() != 4)
      vw_throw(IOErr() << "Index::load_index() -- could not parse blob number from blob filename.");
    std::string blob_id_str(matches[2].first, matches[2].second);
    int current_blob_id = atoi(blob_id_str.c_str());
    
    Blob blob(m_plate_filename + "/" + blob_files[i]);
    Blob::iterator iter = blob.begin();
    while (iter != blob.end()) {
      TileHeader hdr = *iter;
      IndexRecord rec;
      rec.set_blob_id(current_blob_id);
      rec.set_blob_offset(iter.current_base_offset());
      rec.set_status(INDEX_RECORD_VALID);
      rec.set_tile_size(iter.current_data_size());
      m_root->insert(rec, hdr.col(), hdr.row(), hdr.depth(), hdr.epoch());
      ++iter;
    }
  }
}

void vw::platefile::Index::save_index_file() const {

  // First, check to make sure the platefile directory exists.
  if ( !fs::exists( m_plate_filename ) )
    fs::create_directory(m_plate_filename);
  std::ofstream ofstr( this->index_filename().c_str() );
  if (!ofstr.is_open())
    vw_throw(IOErr() << "Index::Index() -- Could not create index file for writing.");

  // Serialize the header data.
  m_header.SerializeToOstream(&ofstr);
  ofstr.close();

}

// ------------------------ Public Methods --------------------------

/// Create a new index.  User supplies a pre-configure blob manager.
vw::platefile::Index::Index( std::string plate_filename,
                             int default_tile_size, std::string default_file_type) :
  m_plate_filename(plate_filename),
  m_blob_manager(boost::shared_ptr<BlobManager>( new BlobManager() )), 
  m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {
  
  // Set up the IndexHeader and write it to disk.
  m_header.set_platefile_version(VW_CURRENT_PLATEFILE_VERSION);
  m_header.set_default_tile_size(default_tile_size);
  m_header.set_default_file_type(default_file_type);
  m_header.set_transaction_read_cursor(0);
  m_header.set_transaction_write_cursor(1);

  this->save_index_file();

  // Create the logging facility
  m_log = boost::shared_ptr<LogInstance>( new LogInstance(this->log_filename()) );

  this->log() << "Created new index \"" << this->index_filename() 
              << "\"  (version " << VW_CURRENT_PLATEFILE_VERSION <<"\n";
  this->log() << "\t--> Tile size: " << m_header.default_tile_size()
              << " and file type: " << m_header.default_file_type() << "\n";
  this->log() << "\t--> Read cursor: " << m_header.transaction_read_cursor() 
              << "   Write cursor: " << m_header.transaction_write_cursor() << "\n";
}

/// Open an existing index from a file on disk.
vw::platefile::Index::Index(std::string plate_filename) : 
  m_plate_filename(plate_filename),                                    
  m_blob_manager(boost::shared_ptr<BlobManager>( new BlobManager() )),
  m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {

  // First, check to make sure the platefile directory exists.
  if ( !fs::exists( plate_filename ) ) 
    vw_throw(IOErr() << "Could not open platefile: \"" << plate_filename << "\" for reading.");

  // The load the index file & parse the IndexHeader protobuffer
  std::ifstream ifstr( this->index_filename().c_str() );
  if (!ifstr.is_open()) 
    vw_throw(IOErr() << "Index::Index() -- Could not load index file for reading.");

  bool worked = m_header.ParseFromIstream(&ifstr);
  if (!worked)
    vw_throw(IOErr() << "Index::Index() -- Could not parse index information from " 
             << this->index_filename()); 

  ifstr.close();

  // Create the logging facility
  m_log = boost::shared_ptr<LogInstance>( new LogInstance(this->log_filename()) );

  this->log() << "Reopened index \"" << this->index_filename() 
              << "\"  (version " << m_header.platefile_version() << "\n";
  this->log() << "\t--> Tile size: " << m_header.default_tile_size()
              << " and file type: " << m_header.default_file_type() << "\n";
  this->log() << "\t--> Read cursor: " << m_header.transaction_read_cursor() 
              << "   Write cursor: " << m_header.transaction_write_cursor() << "\n";
    

  // Load the actual index data
  std::vector<std::string> blob_files = this->blob_filenames();
  this->load_index(blob_files);
}

/// Use this to send data to the index's logfile like this:
///
///   index_instance.log() << "some text for the log...\n";
///
std::ostream& vw::platefile::Index::log () {
  return (*m_log)(InfoMessage, "console");
}


/// Attempt to access a tile in the index.  Throws an
/// TileNotFoundErr if the tile cannot be found.
IndexRecord vw::platefile::Index::read_request(int col, int row, int depth, int epoch) {
  Mutex::Lock lock(m_mutex);
  return m_root->search(col, row, depth, epoch);
}
  
// Writing, pt. 1: Locks a blob and returns the blob id that can
// be used to write a tile.
int vw::platefile::Index::write_request(int size) {  
  return m_blob_manager->request_lock(size); 
}

// Writing, pt. 2: Supply information to update the index and
// unlock the blob id.
void vw::platefile::Index::write_complete(TileHeader const& header, IndexRecord const& record) {
  m_blob_manager->release_lock(record.blob_id()); 

  Mutex::Lock lock(m_mutex);
  m_root->insert(record, header.col(), header.row(), header.depth(), header.epoch());
  m_root->invalidate_records(header.col(), header.row(), header.depth());
}


// Clients are expected to make a transaction request whenever
// they start a self-contained chunk of mosaicking work.  .
int32 vw::platefile::Index::transaction_request(std::string transaction_description) {
  Mutex::Lock lock(m_mutex);

  // Pick the next transaction ID, increment the cursor, and then save
  // the cursor to a file.  (Saving to a file gurantees that we won't
  // accidentally assign two transactions the same ID if the index
  // server crashes and has to be restarted.
  int32 transaction_id = m_header.transaction_write_cursor();
  m_header.set_transaction_write_cursor(m_header.transaction_write_cursor() + 1);
  this->save_index_file();
  this->log() << "Transaction " << m_header.transaction_write_cursor() << " started: " 
          << transaction_description << "\n";
  return transaction_id;
}

// Once a chunk of work is complete, clients can "commit" their
// work to the mosaic by issuding a transaction_complete method.
int32 vw::platefile::Index::transaction_complete(int32 transaction_id) {
  Mutex::Lock lock(m_mutex);

  this->log() << "Transaction " << m_header.transaction_write_cursor() << " finished.\n";
  
}
