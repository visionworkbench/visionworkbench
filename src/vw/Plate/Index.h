#ifndef __VW_PLATEFILE_INDEX_H__
#define __VW_PLATEFILE_INDEX_H__

#include <vector>
#include <boost/smart_ptr.hpp>

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Thread.h>
#include <vw/Plate/Tree.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/IndexRecord.h>

#define VW_PLATE_INDEXRECORD_FILETYPE_SIZE 5
#define VW_PLATE_INDEX_VERSION 2

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------
  //                          INDEX BASE CLASS
  // -------------------------------------------------------------------
  class IndexBase {
  public:
    // Destructor
    virtual ~IndexBase() {}

    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    virtual IndexRecord read_request(int col, int row, int depth) = 0;
  
    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a block.
    virtual int write_request(int size) = 0;

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_complete(int col, int row, int depth, IndexRecord record) = 0;
  };



  // -------------------------------------------------------------------
  //                            PLATE INDEX
  // -------------------------------------------------------------------

  class Index : public IndexBase { 

    int m_index_version;
    int m_max_depth;
    int m_default_block_size;
    char m_default_file_type[4];
    boost::shared_ptr<BlobManager> m_blob_manager;
    boost::shared_ptr<TreeNode<IndexRecord> > m_root;
    Mutex m_mutex;

  public:

    /// Create a new index.  Uses default blob manager.
    Index(int default_block_size, std::string default_file_type) :
      m_index_version(VW_PLATE_INDEX_VERSION), m_max_depth(0), 
      m_default_block_size(default_block_size),
      m_blob_manager( new BlobManager() ), 
      m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {
      strncpy(m_default_file_type, default_file_type.c_str(), 4);
    }

    /// Create a new index.  User supplies a pre-configure blob manager.
    Index( boost::shared_ptr<BlobManager> blob_manager, 
           int default_block_size, std::string default_file_type ) :
      m_index_version(VW_PLATE_INDEX_VERSION), m_max_depth(0), 
      m_default_block_size(default_block_size),
      m_blob_manager(blob_manager), 
      m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {
      strncpy(m_default_file_type, default_file_type.c_str(), 4);
    }

    /// Open an existing index from a file on disk.
    Index(std::string index_filename) : 
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

    /// Destructor
    virtual ~Index() {}

    /// Save an index out to a file on disk.  This serializes the
    /// tree.
    void save(std::string const& filename) {
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

    int version() const { return m_index_version; }
    int max_depth() const { return m_max_depth; }
    int default_block_size() const { return m_default_block_size; }
    std::string default_file_type() const { return m_default_file_type; }

    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    virtual IndexRecord read_request(int col, int row, int depth) {
      Mutex::Lock lock(m_mutex);
      return m_root->search(col, row, depth);
    }
  
    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a block.
    virtual int write_request(int size) {  
      return m_blob_manager->request_lock(size); 
    }

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_complete(int col, int row, int depth, IndexRecord record) {
      m_blob_manager->release_lock(record.blob_id()); 
      {
        Mutex::Lock lock(m_mutex);
        m_root->insert(record, col, row, depth);
        if (depth > m_max_depth) 
          m_max_depth = depth;
      }
    }

    /// Use only for debugging small trees.
    void print() {
      m_root->print();
    }

  };


}} // namespace vw::plate

#endif // __VW_PLATEFILE_INDEX_H__
