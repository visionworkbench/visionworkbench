// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
#include <vw/Plate/RemoteIndex.h>
#include <vw/Plate/common.h>
#include <vw/Plate/ProtoBuffers.pb.h>

/// Constructor
vw::platefile::RemoteIndex::RemoteIndex(std::string const& platefile, std::string const& requestor) :
  m_platefile(platefile), m_queue_name(requestor) {

  // Create a queue and bind it to the index server exchange.
  m_conn.queue_declare(requestor, true, true, false);
  m_conn.queue_bind(requestor, INDEX_EXCHANGE, requestor + ".#");
    
  // Send an IndexOpenRequest to the AMQP index server.
  IndexOpenRequest req;
  req.set_requestor(m_queue_name);
  req.set_plate_filename(m_platefile);
  req.set_mode("readwrite");
  m_conn.basic_publish_protobuf(req, INDEX_EXCHANGE, "index.index_open_request");

  // Wait for a response and parse the result
  IndexOpenReply r = wait_for_response<IndexOpenReply>(m_queue_name + ".index_open_reply");
  m_platefile_id = r.platefile_id();
  m_secret = r.secret();
  m_index_header = r.index_header();
  vw_out(InfoMessage, "plate") << "Opened remote platefile \"" << platefile 
                               << "\"   ID: " << m_platefile_id << "  Secret: " 
                               << m_secret << "\n";
}
  
/// destructor
vw::platefile::RemoteIndex::~RemoteIndex() {
  IndexCloseRequest req;
  req.set_requestor(m_queue_name);
  req.set_platefile_id(m_platefile_id);
  req.set_secret(m_secret);
  m_conn.basic_publish_protobuf(req, INDEX_EXCHANGE, "index.index_close_request");
}
  
/// Attempt to access a tile in the index.  Throws an
/// TileNotFoundErr if the tile cannot be found.
vw::platefile::IndexRecord vw::platefile::RemoteIndex::read_request(int col, int row, 
                                                                    int depth, int epoch) {
  IndexReadRequest req;
  req.set_requestor(m_queue_name);
  req.set_platefile_id(m_platefile_id);
  req.set_secret(m_secret);
  req.set_col(col);
  req.set_row(row);
  req.set_depth(depth);
  req.set_epoch(epoch);
  m_conn.basic_publish_protobuf(req, INDEX_EXCHANGE, m_queue_name);

  IndexReadReply r = wait_for_response<IndexReadReply>(m_queue_name + ".read_reply");
  return r.index_record();
}
  
// Writing, pt. 1: Locks a blob and returns the blob id that can
// be used to write a tile.
int vw::platefile::RemoteIndex::write_request(int size) {
  IndexWriteRequest req;
  req.set_requestor(m_queue_name);
  req.set_platefile_id(m_platefile_id);
  req.set_secret(m_secret);
  req.set_size(size);
  m_conn.basic_publish_protobuf(req, INDEX_EXCHANGE, m_queue_name);

  IndexWriteReply r = wait_for_response<IndexWriteReply>(m_queue_name + ".write_reply");
  return r.blob_id();
}
  
// Writing, pt. 2: Supply information to update the index and
// unlock the blob id.
void vw::platefile::RemoteIndex::write_complete(TileHeader const& header, IndexRecord const& record) {
  IndexWriteComplete req;
  req.set_requestor(m_queue_name);
  req.set_platefile_id(m_platefile_id);
  req.set_secret(m_secret);
  *(req.mutable_header()) = header;
  *(req.mutable_record()) = record;
  m_conn.basic_publish_protobuf(req, INDEX_EXCHANGE, m_queue_name);    

  IndexSuccess r = wait_for_response<IndexSuccess>(m_queue_name + ".write_complete");
}
  
  
vw::int32 vw::platefile::RemoteIndex::version() const { 
  return m_index_header.platefile_version(); 
}

vw::int32 vw::platefile::RemoteIndex::max_depth() const { 
  IndexDepthRequest req;
  req.set_requestor(m_queue_name);
  req.set_platefile_id(m_platefile_id);
  req.set_secret(m_secret);
  m_conn.basic_publish_protobuf(req, INDEX_EXCHANGE, m_queue_name);    

  IndexDepthReply r = wait_for_response<IndexDepthReply>(m_queue_name + ".depth_reply");
  return r.depth();
}

std::string vw::platefile::RemoteIndex::platefile_name() const { 
  return m_platefile;
}

vw::int32 vw::platefile::RemoteIndex::default_tile_size() const { 
  return m_index_header.default_tile_size(); 
}

std::string vw::platefile::RemoteIndex::default_tile_filetype() const { 
  return m_index_header.default_file_type(); 
}

vw::PixelFormatEnum vw::platefile::RemoteIndex::pixel_format() const {
  return vw::PixelFormatEnum(m_index_header.pixel_format()); 
}

vw::ChannelTypeEnum vw::platefile::RemoteIndex::channel_type() const {
  return vw::ChannelTypeEnum(m_index_header.channel_type()); 
}



