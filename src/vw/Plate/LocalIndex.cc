// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Core/ProgressCallback.h>
#include <vw/Core/Log.h>
#include <vw/Plate/LocalIndex.h>
#include <vw/Plate/RemoteIndex.h>
#include <vw/Plate/Blob.h>
using namespace vw;
using namespace vw::platefile;

#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <algorithm>
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
  m_blob_manager(boost::shared_ptr<BlobManager>( new BlobManager() )) {

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
  m_header.set_num_levels(0);                // Index initially contains zero levels

  this->save_index_file();

  // Create the logging facility
  m_log = boost::shared_ptr<LogInstance>( new LogInstance(this->log_filename()) );
  this->log() << "Created new index: \"" << this->index_filename() << "\n"
              <<  m_header.DebugString() << "\n";
}

/// Open an existing index from a file on disk.
vw::platefile::LocalIndex::LocalIndex(std::string plate_filename) :
  m_plate_filename(plate_filename),
  m_blob_manager(boost::shared_ptr<BlobManager>( new BlobManager() )) {

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
    vw_out(WarningMessage, "plate") << "\nWARNING: could not open index log file. "
                                    << "Proceed with caution.";
  }
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

// Clients are expected to make a transaction request whenever
// they start a self-contained chunk of mosaicking work.  .
int32 vw::platefile::LocalIndex::transaction_request(std::string transaction_description,
                                                     int transaction_id_override) {

  if (transaction_id_override != -1) {

    // If the user has chosen to override the transaction ID's, then we
    // use the transaction ID they specify.  We also increment the
    // transaction_write_cursor if necessary so that we dole out a
    // reasonable transaction ID if we are ever asked again without an
    // override.
    int max_trans_id = std::max(m_header.transaction_write_cursor(), transaction_id_override+1);
    m_header.set_transaction_write_cursor(max_trans_id);
    this->save_index_file();
    this->log() << "Transaction " << transaction_id_override 
                << " started: " << transaction_description << "\n";
    return transaction_id_override;

  } else {
  
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
}

// Once a chunk of work is complete, clients can "commit" their
// work to the mosaic by issuing a transaction_complete method.
void vw::platefile::LocalIndex::transaction_complete(int32 transaction_id, bool update_read_cursor) {

  if ( update_read_cursor ) {
    int max_trans_id = std::max(m_header.transaction_read_cursor(), transaction_id);
    m_header.set_transaction_read_cursor(max_trans_id);
    this->save_index_file();
  }

  this->log() << "Transaction " << transaction_id << " complete.  "
              << "[ read_cursor = " << m_header.transaction_read_cursor() << " ]\n";
}

// Once a chunk of work is complete, clients can "commit" their
// work to the mosaic by issuing a transaction_complete method.
void vw::platefile::LocalIndex::transaction_failed(int32 transaction_id) {

  // Log the failure.
  this->log() << "Transaction " << transaction_id << " FAILED.\n";

}

// Return the current location of the transaction cursor.  This
// will be the last transaction id that refers to a coherent
// version of the mosaic.
int32 vw::platefile::LocalIndex::transaction_cursor() {
  return m_header.transaction_read_cursor();
}

