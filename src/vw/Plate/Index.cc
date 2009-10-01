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
    vw_out(InfoMessage, "plate") << "Loading index entries from blob file: " 
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
      rec.set_valid(true);
      m_root->insert(rec, hdr.col(), hdr.row(), hdr.depth(), hdr.epoch());
      ++iter;
    }
  }
}


// ------------------------ Public Methods --------------------------

/// Create a new index.  User supplies a pre-configure blob manager.
vw::platefile::Index::Index( std::string plate_filename,
                             int default_tile_size, std::string default_file_type) :
  m_plate_filename(plate_filename),
  m_blob_manager(boost::shared_ptr<BlobManager>( new BlobManager() )), 
  m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {

  // First, check to make sure the platefile directory exists.
  if ( !fs::exists( plate_filename ) )
    fs::create_directory(plate_filename);
  std::ofstream ofstr( this->index_filename().c_str() );
  
  if (!ofstr.is_open())
    vw_throw(IOErr() << "Index::Index() -- Could not create index file for writing.");

  // Set up the IndexHeader and write it to disk.
  m_header.set_platefile_version(VW_CURRENT_PLATEFILE_VERSION);
  m_header.set_default_tile_size(default_tile_size);
  m_header.set_default_file_type(default_file_type);
  m_header.SerializeToOstream(&ofstr);

  ofstr.close();
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

  // For testing
  std::vector<std::string> blob_files = this->blob_filenames();
  this->load_index(blob_files);


  //----

  // std::ifstream istr(plate_filename.c_str(), std::ios::binary);
  // if (!istr.good())
  //   vw_throw(IOErr() << "Could not open index file \"" << plate_filename << "\" for reading.");

  // // Read the index version and supporting info
  // istr.read( (char*)&m_index_version, sizeof(m_index_version) );
  // if (m_index_version != VW_PLATE_INDEX_VERSION) 
  //   vw_throw(IOErr() << "Could not open plate index.  " 
  //            << "Version " << m_index_version 
  //            << " is not compatible the current version (" 
  //            << VW_PLATE_INDEX_VERSION << ")");

  // // Read the index metadata
  // istr.read( (char*)&m_max_depth, sizeof(m_max_depth) );
  // istr.read( (char*)&m_default_tile_size, sizeof(m_default_tile_size) );
  // // for (unsigned i = 0; i < 4; ++i)
  // //   istr.read( (char*)(m_default_file_type+i), sizeof(*m_default_file_type) );

  // // Read blob manager info and create a blob manager.
  // int num_blobs;
  // int64 max_blob_size;
  // istr.read( (char*)&num_blobs, sizeof(num_blobs) );
  // istr.read( (char*)&max_blob_size, sizeof(max_blob_size) );
  // m_blob_manager = boost::shared_ptr<BlobManager>( new BlobManager(max_blob_size, 
  //                                                                  num_blobs));

  // // Deserialize the index tree
  // m_root->deserialize(istr);

  // istr.close();
}

// /// Save an index out to a file on disk.  This serializes the
// /// tree.
// void vw::platefile::Index::save(std::string const& filename) {
//   std::ofstream ostr(filename.c_str(), std::ios::binary);
//   if (!ostr.good())
//     vw_throw(IOErr() << "Could not open index file \"" << filename << "\" for writing.");

//   // Save basic index information & version.
//   ostr.write( (char*)&m_index_version, sizeof(m_index_version) );
//   ostr.write( (char*)&m_max_depth, sizeof(m_max_depth) );
//   ostr.write( (char*)&m_default_tile_size, sizeof(m_default_tile_size) );
//   for (unsigned i = 0; i < 4; ++i)
//     ostr.write( (char*)(m_default_file_type+i), sizeof(*m_default_file_type) );

//   // Save blob manager information
//   int num_blobs = m_blob_manager->num_blobs();
//   int64 max_blob_size = m_blob_manager->max_blob_size();
//   ostr.write( (char*)&num_blobs, sizeof(num_blobs) );
//   ostr.write( (char*)&max_blob_size, sizeof(max_blob_size) );

//   // Serialize the index tree
//   m_root->serialize(ostr);
//   ostr.close();
// }

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
}

