// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATEFILE_INDEX_H__
#define __VW_PLATEFILE_INDEX_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Thread.h>

#include <vw/Plate/Tree.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/BlobManager.h>

// Protocol Buffer
#include <vw/Plate/IndexRecord.pb.h>

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


    virtual void save(std::string const& filename) = 0;


    virtual int32 version() const = 0;
    virtual int32 default_block_size() const = 0;
    virtual std::string default_block_filetype() const = 0;
    virtual int32 max_depth() const = 0;
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
    Index(int default_block_size, std::string default_file_type);

    /// Create a new index.  User supplies a pre-configure blob manager.
    Index( boost::shared_ptr<BlobManager> blob_manager, 
           int default_block_size, std::string default_file_type );

    /// Open an existing index from a file on disk.
    Index(std::string index_filename);

    /// Destructor
    virtual ~Index() {}

    /// Save an index out to a file on disk.  This serializes the
    /// tree.
    virtual void save(std::string const& filename);

    virtual int version() const { return m_index_version; }
    virtual int max_depth() const { return m_max_depth; }
    virtual int32 default_block_size() const { return m_default_block_size; }
    virtual std::string default_block_filetype() const { return m_default_file_type; }

    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    virtual IndexRecord read_request(int col, int row, int depth);
  
    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a block.
    virtual int write_request(int size);

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_complete(int col, int row, int depth, IndexRecord record);

    /// Use only for debugging small trees.
    void print() {
      m_root->print();
    }

  };


}} // namespace vw::plate

#endif // __VW_PLATEFILE_INDEX_H__
