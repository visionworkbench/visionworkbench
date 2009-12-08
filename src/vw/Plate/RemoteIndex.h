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
    virtual IndexRecord read_request(int col, int row, int depth, int transaction_id, bool exact_transaction_match = false);
  
    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(int size);
  
    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_complete(TileHeader const& header, 
                                IndexRecord const& record);

    virtual IndexHeader index_header() const;
  
    virtual int32 version() const;
    virtual int32 max_depth() const;

    virtual std::string platefile_name() const;

    virtual int32 tile_size() const;
    virtual std::string tile_filetype() const;

    virtual PixelFormatEnum pixel_format() const;
    virtual ChannelTypeEnum channel_type() const;

    // --------------------- TRANSACTIONS ------------------------

    // Clients are expected to make a transaction request whenever
    // they start a self-contained chunk of mosaicking work.  .
    virtual int32 transaction_request(std::string transaction_description,
                                      std::vector<TileHeader> const& tile_headers);

    /// Called right before the beginning of the mipmapping pass
    virtual void root_complete(int32 transaction_id,
                               std::vector<TileHeader> const& tile_headers);

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id);

    virtual int32 transaction_cursor();

  };

}} // namespace vw::plate

#endif // __VW_PLATE_REMOTE_INDEX_H__
