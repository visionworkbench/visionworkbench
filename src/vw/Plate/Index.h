#ifndef __VW_PLATEFILE_INDEX_H__
#define __VW_PLATEFILE_INDEX_H__

#include <vector>
#include <boost/smart_ptr.hpp>

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Thread.h>
#include <vw/Plate/Tree.h>
#include <vw/Plate/Blob.h>

namespace vw {
namespace platefile {
  
  // -------------------------------------------------------------------
  //                          INDEX_RECORD
  // 
  // The IndexRecord stores all of the metadata needed to recover an
  // image from the blob.  It will be marked as valid if it contains
  // good data, or as invalid if it does not.  Records may be marked
  // as invalid if they represent low-res tiles that need to be
  // re-renedered from higher-res tiles that have been updated.
  // -------------------------------------------------------------------
  
  class IndexRecord {

    int32 m_blob_id;
    size_t m_blob_offset;
    size_t m_block_size;
    char m_block_filetype[5];
    bool m_valid;

  public:

    IndexRecord() : m_valid(false) {}

    IndexRecord(int blob_id, size_t blob_offset, 
                size_t block_size, std::string block_filetype) :
      m_blob_id(blob_id), m_blob_offset(blob_offset), 
      m_block_size(block_size), m_valid(true) {

      if (block_filetype.size() > 4) 
        vw_throw(ArgumentErr() << "IndexRecord: filetype argument must be 4 characters or fewer.");
      strncpy(m_block_filetype, block_filetype.c_str(), 5);

    }

    int32 blob_id() const { return m_blob_id; }
    void set_blob_id(int32 blob_id) { m_blob_id = blob_id; }
    
    size_t blob_offset() const { return m_blob_offset; }
    void set_blob_offset(size_t blob_offset) { m_blob_offset = blob_offset; }

    size_t block_size() const { return m_block_size; }
    void set_block_size(size_t block_size) { m_block_size = block_size; }

    std::string block_filetype() const { return m_block_filetype; }
    void set_block_filetype(std::string block_filetype) {
      if (block_filetype.size() > 4) 
        vw_throw(ArgumentErr() << "IndexRecord: filetype argument must be 4 characters or fewer.");
      strncpy(m_block_filetype, block_filetype.c_str(), 5);
    }

    bool valid() const { return m_valid; }
    void set_valid(bool valid) { m_valid = valid; }
  };


  // -------------------------------------------------------------------
  //                            PLATE_FILE
  // -------------------------------------------------------------------

  class Index { 

    boost::shared_ptr<BlobManager> m_blob_manager;
    boost::shared_ptr<TreeNode<IndexRecord> > m_root;
    Mutex m_mutex;

  public:

    /// Create a new index.  Uses default blob manager.
    Index() :
      m_blob_manager( new BlobManager() ), 
      m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {}

    /// Create a new index.  User supplies a pre-configure blob manager.
    Index( boost::shared_ptr<BlobManager> blob_manager) :
      m_blob_manager(blob_manager), 
      m_root(boost::shared_ptr<TreeNode<IndexRecord> >( new TreeNode<IndexRecord>() )) {}

    /// Open an existing index. 
    Index(std::string index_filename) {
      std::cout << "OPENING A PLATE INDEX IS NOT YET IMPLEMENTED!!!\n";
      exit(0);
    }

    void save(std::string filename) {
      std::cout << "SAVING A PLATE INDEX IS NOT YET IMPLEMENTED!!!\n";      
      exit(0);
    }

    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    IndexRecord read_request(int col, int row, int depth) {
      Mutex::Lock lock(m_mutex);
      return m_root->search(col, row, depth);
    }
  
    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a block.
    int write_request(int size) {  
      return m_blob_manager->request_lock(size); 
    }

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    void write_complete(int col, int row, int depth, IndexRecord record) {
      m_blob_manager->release_lock(record.blob_id()); 
      {
        Mutex::Lock lock(m_mutex);
        m_root->insert(record, col, row, depth);
      }
    }
  };

}} // namespace vw::plate

#endif // __VW_PLATEFILE_INDEX_H__
