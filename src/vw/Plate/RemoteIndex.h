// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
#ifndef __VW_PLATE_REMOTE_INDEX_H__
#define __VW_PLATE_REMOTE_INDEX_H__

#include <vw/Plate/Index.h>
#include <vw/Plate/Amqp.h>

namespace vw {
namespace platefile {

  class RemoteIndex : public IndexBase {
  
    std::string m_platefile;
    std::string m_queue_name;
    AmqpConnection m_conn;
  
    int32 m_platefile_id;
    int32 m_secret;
    IndexHeader m_index_header;


    // Wait for a response from the AMQP server.  If the routing key of
    // the response matches the expected_routing_key, then we parse and
    // return a protobuffer of type ProtoBufT.  If we get an error or
    // unexpeted reply, we throw an exception.
    //
    // TODO: Add the ability to set a wait timeout that throws a
    // AmqpTimeoutErr exception if we wait too long.
    //
    // ALSO TODO: Add a send_and_wait() method that catches these
    // exceptions and retries on failure.
    template <class ProtoBufT>
    ProtoBufT wait_for_response(std::string expected_routing_key) const {
      std::string routing_key;
      std::string response = m_conn.basic_consume(m_queue_name, routing_key, false); 
      if (routing_key == expected_routing_key) {
        ProtoBufT r;
        r.ParseFromString(response);
        return r;
      } else if (routing_key == m_queue_name + ".index_error") {
        IndexError r;
        r.ParseFromString(response);
        vw_throw(LogicErr() << "An index server error occured: " << r.message());
      } else {
        vw_throw(LogicErr() << "Unpected routing key: " 
                 << routing_key << " (expected IndexOpenReply)\n");
      }
    }

  public:
    /// Constructor
    RemoteIndex(std::string const& platefile, std::string const& requestor);

    /// destructor
    virtual ~RemoteIndex();
  
    /// Attempt to access a tile in the index.  Throws an
    /// TileNotFoundErr if the tile cannot be found.
    virtual IndexRecord read_request(int col, int row, int depth, int epoch = 0);
  
    // Writing, pt. 1: Locks a blob and returns the blob id that can
    // be used to write a tile.
    virtual int write_request(int size);
  
    // Writing, pt. 2: Supply information to update the index and
    // unlock the blob id.
    virtual void write_complete(TileHeader const& header, IndexRecord const& record);
  
    virtual int32 version() const;
    virtual int32 max_depth() const;

    virtual std::string platefile_name() const;

    virtual int32 default_tile_size() const;
    virtual std::string default_tile_filetype() const;

    virtual PixelFormatEnum pixel_format() const;
    virtual ChannelTypeEnum channel_type() const;

  };

}} // namespace vw::plate

#endif // __VW_PLATE_REMOTE_INDEX_H__
