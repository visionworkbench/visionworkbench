#ifndef __VW_PLATEFILE_INDEX_H__
#define __VW_PLATEFILE_INDEX_H__

#include <vector>
#include <boost/smart_ptr.hpp>

namespace vw {
namespace plate {


  /// TileNotFound exception
  VW_DEFINE_EXCEPTION(TileNotFoundErr, Exception);

  
  // -------------------------------------------------------------------
  //                    IndexRecord and IndexNode
  // -------------------------------------------------------------------
  
  class IndexRecord {

    int32 m_blob_id;
    size_t m_blob_offset;
    size_t m_block_size;
    char m_block_filetype[5];

  public:
    IndexRecord(int blob_id, size_t blob_offset, size_t block_size, std::string block_filetype) :
      m_blob_id(blob_id), m_blob_offset(blob_offset), m_block_size(block_size) {
      if (block_filetype.size() > 4) 
        vw_throw(ArgumentErr() << "IndexRecord: filetype argument must be 4 characters or fewer.");
      strncpy(m_block_filetype, block_filetype.c_str(), 5);
    }

    int32 blob_id() const { return m_blob_id; }
    size_t blob_offset() const { return m_blob_offset; }
    size_t block_size() const { return m_block_size; }
    const char* block_filetype() const { return m_block_filetype; }
  };

    
  class IndexNode {

    boost::shared_ptr<IndexNode> m_parent;
    std::vector<boost::shared_ptr<IndexNode> > m_children;
    IndexRecord m_record;

    // Search for a node at a given col, row, and level.
    IndexRecord search_impl(int col, int row, int level, int current_level) {

      // If we have reached the requested depth, our search terminates
      // here.
      if (current_level == level) {
        return m_record;
        
      // Otherwise, we go recurse deeper into the tree.
      } else {
        int tile_x = col / pow(2,level-current_level);
        int tile_y = row / pow(2,level-current_level);
        int child_id;
        if (tile_x == 0 && tile_y == 0) 
          child_id = 0;
        else if (tile_x == 1 && tile_y == 0) 
          child_id = 1;
        else if (tile_x == 0 && tile_y == 1) 
          child_id = 2;
        else 
          child_id = 3;
          
        if (m_children[child_id]) 
          // If the tile is found, we dive deeper.
          return search_impl(col, row, level, current_level+1);
        else
          // If not, we throw an exception
          vw_throw(TileNotFoundErr() << "Tile search [" col << " " << row << " " << depth 
                   << "] failed at depth " << current_level << "\n");
      }
    }

  public: 
    IndexNode(const boost::shared_ptr<IndexNode>& parent, IndexRecord const& record) :
      m_parent(parent), m_record(record) {
      m_children.resize(4);
    }

    // Return the child of this node with the 'id' according to the
    // following index scheme:
    //
    //    |---|---|
    //    | 0 | 1 |
    //    |---+---|
    //    | 2 | 3 |
    //    |---|---|
    //
    boost::shared_ptr<QTreeNode> child(int id) const { return m_children(id); }

    /// Return the parent of this node.
    boost::shared_ptr<QTreeNode> parent() const { return m_parent; }

    // Sets the child of this node with the 'id' according to the above
    // index scheme.
    //
    void set_child(int id, boost::shared_ptr<QTreeNode> node) {
      m_children(id) = node;
    }

    // Insert the child of this node, but preserves the previous child's descendents.
    //
    void insert_child(int id, boost::shared_ptr<QTreeNode> node) {
      
      // First we save the old child and replace it with the new one.
      boost::shared_ptr<QTreeNode> old_child = m_children[id];
      m_children[id] = node;

      // Then, if the old child existed, we transfer the old child's
      // children (our grandchildren) to the new child node.
      if (old_child)
        for (int i = 0; i < 4 ; ++i) 
          m_children[id]->set_child(i, old_child->child(i));
    }

    IndexRecord search(int col, int row, int level) {
      return search_impl(col, row, level, -1);
    }
         
  };

  // -------------------------------------------------------------------
  //                            BLOB_MANAGER
  // -------------------------------------------------------------------

  /// The BlobManager keeps track of how much data has been stored in
  /// each blob so far.  This allows it to pick blobs that have enough
  /// space to store a new block of a given size.  The blob manager also
  /// controls the locking/unlocking of blobs, and can load balance
  /// blobs writes by alternating which blob is offered up for writing
  /// data.
  ///
  /// The BlobManager is thread safe.
  class BlobManager {
    size_t m_max_blob_size;
    m_vector<int64> m_blob_write_indices;
    m_vector<bool> m_blob_locks;
    int m_blob_index = 0
      Mutex m_mutex;
    Condition m_blob_release_condition;

    void next_blob_index() {
      // Move to the next blob
      m_blob_index++;
      if (m_blob_index >= m_blob_locks.size())
        m_blob_index = 0;
    }

  public:

    BlobManager(size_t max_blob_size, int nblobs = 2) : 
      m_max_blob_size(max_blob_size) { 

      m_blob_write_indices.resize(nblobs);
      m_blob_locks.resize(nblobs);

      for (int i=0; i < m_blob_locks.size(); ++i) {
        m_blob_locks[i] = false;
        m_blob_write_indices[i] = 0;
      }

    }

    /// Return the number of blobs currently in use.
    int num_blobs() const {
      Mutex::Lock m_lock(m_mutex);
      return m_blob_locks.size();
    }

    /// Request a blob to write to that has sufficient space to write at
    /// least 'size' bytes.  Returns the blob index of a locked blob
    /// that you have sole access to write to.
    //
    // TODO: This is pretty simple logic, and would not be very
    // efficient because it blocks on write if it catches up to a blob
    // that is still locked.  We should add real blob selection logic
    // here at a later date.
    int request_lock(size_t size) {
      Mutex::Lock m_lock(m_mutex);

      // First, we check to see if the next blob is free.  If not, we
      // wait for a release event and recheck.
      while(m_blob_locks[m_blob_index] != false)
        m_blob_release_condition.wait(m_lock);

      // Then we lock it, increment the blob index, and return it.
      m_blob_locks[m_blob_index] = true;      
      int idx = m_blob_index;
      next_blob_index();
      return idx;
    }

    // Release the blob lock and update its write index (essentially
    // "committing" the write to the blob when you are finished with it.).
    int release(int blob_id, int size) {
      Mutex::Lock m_lock(m_mutex);
      m_blob_locks[blob_id] = false;
      m_blob_write_indices[blob_id] += size;
      m_blob_release_condition.notify_all();
    }

  };

  // -------------------------------------------------------------------
  //                            PLATE_FILE
  // -------------------------------------------------------------------


  class Index { 

    int32 m_depth;

    BlobManager m_blob_manager;
    boost::shared_ptr<QTreeNode> m_root;
    Mutex m_mutex;

    void insert( IndexRecord record, int col, int row, int depth) {
      // Implement me!!
    }

  public:

    Index(size_t blob_size) :
      m_blob_manager(blob_size) {
      // Initialize the tree
    }

    Index(std::string index_filename) {
      std::cout << "OPENING A PLATE INDEX IS NOT YET IMPLEMENTED!!!\n";
    }

    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    IndexRecord read_request(int col, int row, int depth) {
      Mutex::Lock lock(m_mutex);
      return m_root.search(col, row, depth, -1);
    }
  
    // Writing, pt. 1: Locks a blob and returns the blob id that can be
    // used to write a block.
    int write_request(int size) {  
      m_blob_manager.request_lock(size); 
    }

    // Writing, pt. 2: Supply information to update the index and unlock
    // the blob id.
    void write_complete(int col, int row, int depth, 
                        int blob_id, int blob_offset,
                        int block_size, std::string block_filetype) {  
      m_blob_manager.unlock(blob_id); 

      Mutex::Lock lock(m_mutex);
      IndexRecord record(blob_id, blob_offset, block_size, block_filetype);
      this->insert(record, col, row, depth);
    }


  
  };

}} // namespace vw::plate

#endif // __VW_PLATEFILE_INDEX_H__
