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

  class RemoteIndex {//: public IndexBase {
  
    std::string m_queue_name;
    
    // Remote connection
    boost::shared_ptr<AmqpRpcChannel> m_rpc_channel;
    boost::shared_ptr<AmqpRpcController> m_rpc_controller;
    boost::shared_ptr<IndexService> m_index_service;
  
    int32 m_secret;

  public:
    /// Constructor
    RemoteIndex(std::string const& requestor);

    /// destructor
    virtual ~RemoteIndex();
  
    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    virtual IndexRecord read_request(int platefile_id, int col, int row, 
                                     int depth, int transaction_id);
  
    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(int size);
  
    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_complete(TileHeader const& header, 
                                IndexRecord const& record);
  
    virtual int32 version() const;
    virtual int32 max_depth() const;

    virtual std::string platefile_name(int platefile_id) const;

    virtual int32 default_tile_size() const;
    virtual std::string default_tile_filetype() const;

    virtual PixelFormatEnum pixel_format() const;
    virtual ChannelTypeEnum channel_type() const;

    // --------------------- TRANSACTIONS ------------------------

    // Clients are expected to make a transaction request whenever
    // they start a self-contained chunk of mosaicking work.  .
    virtual int32 transaction_request(std::string transaction_description);

    // Once a chunk of work is complete, clients can "commit" their
    // work to the mosaic by issuding a transaction_complete method.
    virtual void transaction_complete(int32 transaction_id);

    virtual int32 transaction_cursor();

  };

}} // namespace vw::plate

#endif // __VW_PLATE_REMOTE_INDEX_H__
