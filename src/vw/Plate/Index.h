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
#include <vw/Core/Log.h>

#include <vw/Image/PixelTypeInfo.h>

#include <vw/Plate/Tree.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/BlobManager.h>

// Protocol Buffer
#include <vw/Plate/ProtoBuffers.pb.h>

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

    // -------------------------- I/O ---------------------------

    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    virtual IndexRecord read_request(int col, int row, int depth, int transaction_id = -1) = 0;
  
    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(int size) = 0;

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_complete(TileHeader const& header, IndexRecord const& record) = 0;


    // ----------------------- PROPERTIES  ----------------------

    virtual int32 version() const = 0;
    virtual int32 max_depth() const = 0;

    virtual std::string platefile_name() const = 0;

    virtual int32 default_tile_size() const = 0;
    virtual std::string default_tile_filetype() const = 0;

    virtual PixelFormatEnum pixel_format() const = 0;
    virtual ChannelTypeEnum channel_type() const = 0;


    // --------------------- TRANSACTIONS ------------------------

    // Clients are expected to make a transaction request whenever
    // they start a self-contained chunk of mosaicking work.  .
    virtual int32 transaction_request(std::string transaction_description) = 0;

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id) = 0;

    virtual int32 transaction_cursor() = 0;


    // --------------------- UTILITIES ------------------------

    /// Iterate over all nodes in a tree, calling func for each
    /// location.  Note: this will only be implemented for local
    /// indexes.  This function will throw an error if called on a
    /// remote index.
    virtual void map(boost::shared_ptr<TreeMapFunc> func) { 
      vw_throw(NoImplErr() << "IndexBase::map() not implemented for this index type.");
    }
  };

  // -------------------------------------------------------------------
  //                            PLATE INDEX
  // -------------------------------------------------------------------

  class Index : public IndexBase { 
    
    std::string m_plate_filename;
    IndexHeader m_header;
    boost::shared_ptr<BlobManager> m_blob_manager;
    boost::shared_ptr<TreeNode<IndexRecord> > m_root;
    boost::shared_ptr<vw::LogInstance> m_log;
    Mutex m_mutex;

    void save_index_file() const;
    std::string index_filename() const;
    std::string log_filename() const;
    std::vector<std::string> blob_filenames() const;
    void load_index(std::vector<std::string> const& blob_files);

  public:

    /// Create a new, empty index.
    Index( std::string plate_filename, 
           int default_tile_size, std::string default_file_type,
           PixelFormatEnum default_pixel_format,
           ChannelTypeEnum default_channel_type);

    /// Open an existing index from a file on disk.
    Index(std::string plate_filename);

    /// Destructor
    virtual ~Index() {}

    /// Use this to send data to the index's logfile like this:
    ///
    ///   index_instance.log() << "some text for the log...\n";
    ///
    std::ostream& log ();

    // /// Save an index out to a file on disk.  This serializes the
    // /// tree.
    // virtual void save(std::string const& filename);

    virtual int version() const { return m_header.platefile_version(); }
    virtual int32 max_depth() const { return m_root->max_depth(); }
    
    virtual std::string platefile_name() const { return m_plate_filename; }

    virtual int32 default_tile_size() const { return m_header.default_tile_size(); }
    virtual std::string default_tile_filetype() const { return m_header.default_file_type(); }

    virtual PixelFormatEnum pixel_format() const { 
      return PixelFormatEnum(m_header.pixel_format()); 
    }

    virtual ChannelTypeEnum channel_type() const {
      return ChannelTypeEnum(m_header.channel_type());
    }

    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    virtual IndexRecord read_request(vw::int32 col, vw::int32 row, 
                                     vw::int32 depth, vw::int32 transaction_id = -1);
  
    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(vw::int32 size);

    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_complete(TileHeader const& header, IndexRecord const& record);

    // Clients are expected to make a transaction request whenever
    // they start a self-contained chunk of mosaicking work.  .
    virtual int32 transaction_request(std::string transaction_description);

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id);

    // Return the current location of the transaction cursor.  This
    // will be the last transaction id that refers to a coherent
    // version of the mosaic.
    virtual int32 transaction_cursor();

    /// Use only for debugging small trees.
    void print() { m_root->print(); }

    virtual void map(boost::shared_ptr<TreeMapFunc> func) { m_root->map(func); }

  };


}} // namespace vw::plate

#endif // __VW_PLATEFILE_INDEX_H__
