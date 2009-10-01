// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/Index.h>

using namespace vw::platefile;

// -------------------------------------------------------------------
//                            PLATE INDEX
// -------------------------------------------------------------------

/// Create a new index.  Uses default blob manager.
vw::platefile::Index::Index(int default_block_size, std::string default_file_type) :
  m_index_version(VW_PLATE_INDEX_VERSION), m_max_depth(0), 
  m_default_block_size(default_block_size),
  m_blob_manager( new BlobManager() ), 
  m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {
  strncpy(m_default_file_type, default_file_type.c_str(), 4);
}

/// Create a new index.  User supplies a pre-configure blob manager.
vw::platefile::Index::Index( boost::shared_ptr<BlobManager> blob_manager, 
                             int default_block_size, std::string default_file_type ) :
  m_index_version(VW_PLATE_INDEX_VERSION), m_max_depth(0), 
  m_default_block_size(default_block_size),
  m_blob_manager(blob_manager), 
  m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {
  strncpy(m_default_file_type, default_file_type.c_str(), 4);
}

/// Open an existing index from a file on disk.
vw::platefile::Index::Index(std::string index_filename) : 
  m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {
  std::ifstream istr(index_filename.c_str(), std::ios::binary);
  if (!istr.good())
    vw_throw(IOErr() << "Could not open index file \"" << index_filename << "\" for reading.");

  // Read the index version and supporting info
  istr.read( (char*)&m_index_version, sizeof(m_index_version) );
  if (m_index_version != VW_PLATE_INDEX_VERSION) 
    vw_throw(IOErr() << "Could not open plate index.  " 
             << "Version " << m_index_version 
             << " is not compatible the current version (" 
             << VW_PLATE_INDEX_VERSION << ")");

  // Read the index metadata
  istr.read( (char*)&m_max_depth, sizeof(m_max_depth) );
  istr.read( (char*)&m_default_block_size, sizeof(m_default_block_size) );
  for (unsigned i = 0; i < 4; ++i)
    istr.read( (char*)(m_default_file_type+i), sizeof(*m_default_file_type) );

  // Read blob manager info and create a blob manager.
  int num_blobs;
  int64 max_blob_size;
  istr.read( (char*)&num_blobs, sizeof(num_blobs) );
  istr.read( (char*)&max_blob_size, sizeof(max_blob_size) );
  m_blob_manager = boost::shared_ptr<BlobManager>( new BlobManager(max_blob_size, 
                                                                   num_blobs));

  // Deserialize the index tree
  m_root->deserialize(istr);

  istr.close();
}

/// Save an index out to a file on disk.  This serializes the
/// tree.
void vw::platefile::Index::save(std::string const& filename) {
  std::ofstream ostr(filename.c_str(), std::ios::binary);
  if (!ostr.good())
    vw_throw(IOErr() << "Could not open index file \"" << filename << "\" for writing.");

  // Save basic index information & version.
  ostr.write( (char*)&m_index_version, sizeof(m_index_version) );
  ostr.write( (char*)&m_max_depth, sizeof(m_max_depth) );
  ostr.write( (char*)&m_default_block_size, sizeof(m_default_block_size) );
  for (unsigned i = 0; i < 4; ++i)
    ostr.write( (char*)(m_default_file_type+i), sizeof(*m_default_file_type) );

  // Save blob manager information
  int num_blobs = m_blob_manager->num_blobs();
  int64 max_blob_size = m_blob_manager->max_blob_size();
  ostr.write( (char*)&num_blobs, sizeof(num_blobs) );
  ostr.write( (char*)&max_blob_size, sizeof(max_blob_size) );

  // Serialize the index tree
  m_root->serialize(ostr);
  ostr.close();
}

/// Attempt to access a tile in the index.  Throws an
/// TileNotFoundErr if the tile cannot be found.
IndexRecord vw::platefile::Index::read_request(int col, int row, int depth) {
  Mutex::Lock lock(m_mutex);
  return m_root->search(col, row, depth);
}
  
// Writing, pt. 1: Locks a blob and returns the blob id that can
// be used to write a block.
int vw::platefile::Index::write_request(int size) {  
  return m_blob_manager->request_lock(size); 
}

// Writing, pt. 2: Supply information to update the index and
// unlock the blob id.
void vw::platefile::Index::write_complete(int col, int row, int depth, IndexRecord record) {
  m_blob_manager->release_lock(record.blob_id()); 
  {
    Mutex::Lock lock(m_mutex);
    m_root->insert(record, col, row, depth);
    if (depth > m_max_depth) 
      m_max_depth = depth;
  }
}

