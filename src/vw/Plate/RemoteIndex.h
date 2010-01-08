// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_REMOTE_INDEX_H__
#define __VW_PLATE_REMOTE_INDEX_H__

#include <vw/Plate/Index.h>
#include <vw/Plate/Amqp.h>
#include <vw/Plate/RpcServices.h>

namespace vw {
namespace platefile {

  class RemoteIndex : public Index {
  
    int m_platefile_id;
    IndexHeader m_index_header;
    std::string m_short_plate_filename;
    std::string m_full_plate_filename;

    // Remote connection
    boost::shared_ptr<AmqpRpcClient> m_rpc_controller;
    boost::shared_ptr<IndexService> m_index_service;
  
  public:
    /// Constructor (for opening an existing index)
    RemoteIndex(std::string const& url);

    /// Constructor (for creating a new index)
    RemoteIndex(std::string const& url, IndexHeader new_index_info);

    /// destructor
    virtual ~RemoteIndex();

    /// Sync any unsaved data in the index to disk.
    virtual void sync() {
      vw_throw(NoImplErr() << "Error: sync() is not implemented for remote indices.");
    }
  
    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    ///
    /// By default, this call to read will return a tile with the MOST
    /// RECENT transaction_id <= to the transaction_id you specify
    /// here in the function arguments (if a tile exists).  However,
    /// setting exact_transaction_match = true will force the
    /// PlateFile to search for a tile that has the EXACT SAME
    /// transaction_id as the one that you specify.
    ///
    /// A transaction ID of -1 indicates that we should return the
    /// most recent tile, regardless of its transaction id.
    virtual IndexRecord read_request(int col, int row, int level, int transaction_id, bool exact_transaction_match = false);

    /// Return multiple index entries that match the specified
    /// transaction id range.  This range is inclusive of the first
    /// entry and the last entry: [ begin_transaction_id, end_transaction_id ]
    ///
    /// Results are return as a std::pair<int32, IndexRecord>.  The
    /// first value in the pair is the transaction id for that
    /// IndexRecord.
    virtual std::list<std::pair<int32, IndexRecord> > multi_read_request(int col, int row, int level, 
                                                                         int begin_transaction_id, 
                                                                         int end_transaction_id);

    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(int size);
  
    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_update(TileHeader const& header, IndexRecord const& record);

    /// Writing, pt. 3: Signal the completion 
    virtual void write_complete(int blob_id, uint64 blob_offset);

    // ----------------------- PROPERTIES  ----------------------

    /// Returns a list of tile headers for any valid tiles that exist
    /// at a the specified level and transaction_id.  The
    /// transaction_id is treated the same as it would be for
    /// read_request() above.  The region specifies a tile range of
    /// interest.
    virtual std::list<TileHeader> valid_tiles(int level, BBox2i const& region,
                                              int begin_transaction_id, 
                                              int end_transaction_id, 
                                              int min_num_matches) const;

    virtual IndexHeader index_header() const;
  
    virtual int32 version() const;
    virtual int32 num_levels() const;

    virtual std::string platefile_name() const;

    virtual int32 tile_size() const;
    virtual std::string tile_filetype() const;

    virtual PixelFormatEnum pixel_format() const;
    virtual ChannelTypeEnum channel_type() const;

    // --------------------- TRANSACTIONS ------------------------

    // Clients are expected to make a transaction request whenever
    // they start a self-contained chunk of mosaicking work.  .
    virtual int32 transaction_request(std::string transaction_description,
                                      int transaction_id_override);

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id, bool update_read_cursor);

    // If a transaction fails, we may need to clean up the mosaic.  
    virtual void transaction_failed(int32 transaction_id);

    virtual int32 transaction_cursor();

  };

}} // namespace vw::plate

#endif // __VW_PLATE_REMOTE_INDEX_H__
