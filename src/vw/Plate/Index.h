#ifndef __VW_PLATEFILE_INDEX_H__
#define __VW_PLATEFILE_INDEX_H__

#include <vector>
#include <boost/smart_ptr.hpp>

#include <vw/Core/Exception.h>


namespace vw {
namespace platefile {

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
    bool m_valid;

  public:
    IndexRecord(int blob_id, size_t blob_offset, size_t block_size, std::string block_filetype) :
      m_blob_id(blob_id), m_blob_offset(blob_offset), m_block_size(block_size), m_valid(true) {
      if (block_filetype.size() > 4) 
        vw_throw(ArgumentErr() << "IndexRecord: filetype argument must be 4 characters or fewer.");
      strncpy(m_block_filetype, block_filetype.c_str(), 5);
    }

    IndexRecord() : m_valid(false) {}

    int32 blob_id() const { return m_blob_id; }
    size_t blob_offset() const { return m_blob_offset; }
    size_t block_size() const { return m_block_size; }
    const char* block_filetype() const { return m_block_filetype; }
    bool valid() const { return m_valid; }
  };

    
  class IndexNode {

    IndexNode *m_parent;
    std::vector<boost::shared_ptr<IndexNode> > m_children;
    int m_depth;
    IndexRecord m_record;

    // ------------------------ Private Methods -----------------------------

    // Do some quick arithmetic here to determine which child to follow
    // in the quad tree for a given col, row, level request.
    int compute_child_id(int col, int row, int level) {

      int tile_x = col / powl(2,level-m_depth);
      int tile_y = row / powl(2,level-m_depth);
      int child_id;
      if (tile_x == 0 && tile_y == 0) 
        return 0;
      else if (tile_x == 1 && tile_y == 0) 
        return 1;
      else if (tile_x == 0 && tile_y == 1) 
        return 2;
      else 
        return 3;
    }

  public: 

    // ------------------------ Public Methods -----------------------------

    /// Use this contructor for the root of the tree.  This sets up
    /// the correct depth (and no parent node).
    IndexNode(): m_parent(NULL), m_depth(-1) {
      m_children.resize(4);
    }

    /// Use this constructor if you intend to add actual record data later on.
    IndexNode(IndexNode *parent) :
      m_parent(parent), m_depth(parent->depth() + 1) {
      m_children.resize(4);
    }

    /// Use this contstructor to add record data immediately.
    IndexNode(IndexNode *parent, IndexRecord const& record) :
      m_parent(parent), m_depth(parent->depth() + 1), m_record(record) {
      m_children.resize(4);
    }
    
    int depth() const { return m_depth; }

    // Return the child of this node with the 'id' according to the
    // following index scheme:
    //
    //    |---|---|
    //    | 0 | 1 |
    //    |---+---|
    //    | 2 | 3 |
    //    |---|---|
    //
    boost::shared_ptr<IndexNode> child(int id) const { return m_children[id]; }

    // Sets the child of this node with the 'id' according to the above
    // index scheme.
    //
    void set_child(int id, boost::shared_ptr<IndexNode> node) {
      m_children[id] = node;
    }

    // Insert the child of this node, but preserves the previous child's descendents.
    //
    void insert_child(int id, boost::shared_ptr<IndexNode> node) {
      
      // First we save the old child and replace it with the new one.
      boost::shared_ptr<IndexNode> old_child = m_children[id];
      m_children[id] = node;

      // Then, if the old child existed, we transfer the old child's
      // children (our grandchildren) to the new child node.
      if (old_child)
        for (int i = 0; i < 4 ; ++i) 
          m_children[id]->set_child(i, old_child->child(i));
    }

    // Search for a node at a given col, row, and level.
    IndexRecord search(int col, int row, int level) {

      // If we have reached the requested depth, then we must be at
      // the node we want!  Return the IndexRecord!
      if (m_depth == level) {
        return m_record;
        
      // Otherwise, we go recurse deeper into the tree....
      } else {
        
        int child_id = compute_child_id(col, row, level);
        
        if (m_children[child_id]) {

          // If the tile is found, we dive deeper.
          return m_children[child_id]->search(col, row, level);

        } else {
          TileNotFoundErr err;

          // If not, we throw an exception.
          vw_throw(TileNotFoundErr() << "Tile search [" << col << " " << row << " " << m_depth 
                   << "] failed at depth " << m_depth << "\n");

        }
      }
    }

    // Insert an IndexRecord at a given position.  Intermediate nodes
    // in the tree are created (with empty IndexRecords) in the tree
    // along the way, as needed.
    void insert(IndexRecord const& record, int col, int row, int level) {

      // If we have reached the requested depth, then we must be at
      // the node we want!  Return the IndexRecord!
      if (m_depth == level) {

        // TODO: This is where we could keep the history of
        // IndexRecords for this node.  Let's add that feature
        // someday!!
        m_record = record;

      // Otherwise, we need to recurse further into the tree....
      } else {
        int child_id = this->compute_child_id(col, row, level);
        
        // If the child we need is not yet created, we create it, add
        //it as our child, and recurse down it.
        if (!m_children[child_id]) {
          boost::shared_ptr<IndexNode> node(new IndexNode(this));
          this->set_child(child_id, node);
        }
         
        m_children[child_id]->insert(record, col, row, level);
      }

      
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
    std::vector<int64> m_blob_write_indices;
    std::vector<bool> m_blob_locks;
    int m_blob_index;
    Mutex m_mutex;
    Condition m_blob_release_condition;

    void next_blob_index() {
      // Move to the next blob
      m_blob_index++;
      if (m_blob_index >= m_blob_locks.size())
        m_blob_index = 0;
    }

  public:

    /// Create a new blob manager.  The max_blob_size is specified in
    /// units of megabytes.
    BlobManager(size_t max_blob_size, int nblobs = 2) : 
      m_max_blob_size(max_blob_size), m_blob_index(0) { 

      m_blob_write_indices.resize(nblobs);
      m_blob_locks.resize(nblobs);

      for (int i=0; i < m_blob_locks.size(); ++i) {
        m_blob_locks[i] = false;
        m_blob_write_indices[i] = 0;
      }

    }

    /// Return the number of blobs currently in use.
    int num_blobs() {
      Mutex::Lock lock(m_mutex);
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
      Mutex::Lock lock(m_mutex);

      // First, we check to see if the next blob is free.  If not, we
      // wait for a release event and recheck.
      while(m_blob_locks[m_blob_index] != false)
        m_blob_release_condition.wait(lock);

      // Then we lock it, increment the blob index, and return it.
      m_blob_locks[m_blob_index] = true;      
      int idx = m_blob_index;
      next_blob_index();
      return idx;
    }

    // Release the blob lock and update its write index (essentially
    // "committing" the write to the blob when you are finished with it.).
    int release(int blob_id, int size) {
      Mutex::Lock lock(m_mutex);
      m_blob_locks[blob_id] = false;
      m_blob_write_indices[blob_id] += size;
      m_blob_release_condition.notify_all();
    }

  };

  // -------------------------------------------------------------------
  //                            PLATE_FILE
  // -------------------------------------------------------------------


  class Index { 

    BlobManager m_blob_manager;
    boost::shared_ptr<IndexNode> m_root;
    Mutex m_mutex;

  public:

    /// Create a new index.  Blob size is specified in Megabytes.
    Index(size_t blob_size = 2048) :
      m_blob_manager(blob_size), m_root(boost::shared_ptr<IndexNode>(new IndexNode())) {}

    Index(std::string index_filename) : m_blob_manager(2048) {
      std::cout << "OPENING A PLATE INDEX IS NOT YET IMPLEMENTED!!!\n";
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
      m_blob_manager.request_lock(size); 
    }

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    void write_complete(int col, int row, int depth, 
                        int blob_id, int blob_offset,
                        int block_size, std::string block_filetype) {  
      m_blob_manager.release(blob_id, block_size); 

      Mutex::Lock lock(m_mutex);
      IndexRecord record(blob_id, blob_offset, block_size, block_filetype);
      m_root->insert(record, col, row, depth);
    }
  };

}} // namespace vw::plate

#endif // __VW_PLATEFILE_INDEX_H__
