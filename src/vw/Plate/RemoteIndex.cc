// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/RemoteIndex.h>
#include <vw/Plate/common.h>
#include <vw/Plate/ProtoBuffers.pb.h>


void null_closure() {}

/// Constructor
vw::platefile::RemoteIndex::RemoteIndex(std::string const& requestor) :
  m_queue_name(requestor) {

  // Set up the connection to the AmqpRpcService
  m_rpc_channel.reset( new AmqpRpcChannel(INDEX_EXCHANGE, "index", requestor) );
  m_rpc_controller.reset ( new AmqpRpcController() );
  m_index_service.reset ( new IndexService::Stub(m_rpc_channel.get()) );
  
  // Send an IndexOpenRequest to the AMQP index server.
  // IndexOpenRequest request;
  // request.set_plate_name(platefile);
  // request.set_mode("readwrite");

  // vw_out(0) << "ISSUING IndexOpenRequest\n" << request.DebugString() << "\n";

  // IndexOpenReply response;
  // m_index_service->OpenRequest(m_rpc_controller.get(), &request, &response, 
  //                              google::protobuf::NewCallback(&null_closure));

  // if (m_rpc_controller->Failed())
  //   vw_throw(IOErr() << "OpenRequest RPC failed: " << m_rpc_controller->ErrorText());

  // vw_out(0) << "RECEIVED IndexOpenReply\n" << response.DebugString() << "\n";

  //  m_platefile_id = response.platefile_id();
  // m_secret = response.secret();
  // vw_out(InfoMessage, "plate") << "Opened remote platefile \"" << platefile 
  //                              << "\"   ID: " << m_platefile_id << "  Secret: " 
  //                              << m_secret << "\n";
}
  
/// destructor
vw::platefile::RemoteIndex::~RemoteIndex() {
  // IndexCloseRequest request;
  // request.set_platefile_id(m_platefile_id);
  // request.set_secret(m_secret);
  //  m_conn.basic_publish_protobuf(req, INDEX_EXCHANGE, "index.index_close_request");
}
  
/// Attempt to access a tile in the index.  Throws an
/// TileNotFoundErr if the tile cannot be found.
vw::platefile::IndexRecord vw::platefile::RemoteIndex::read_request(int platefile_id, 
                                                                    int col, int row, 
                                                                    int depth, 
                                                                    int transaction_id) {
  IndexReadRequest request;
  request.set_platefile_id(platefile_id);
  request.set_secret(m_secret);
  request.set_col(col);
  request.set_row(row);
  request.set_depth(depth);
  request.set_transaction_id(transaction_id);
  vw_out(0) << "ISSUING IndexReadRequest\n" << request.DebugString() << "\n";

  IndexReadReply response;
  m_index_service->ReadRequest(m_rpc_controller.get(), &request, &response, 
                               google::protobuf::NewCallback(&null_closure));
  if (m_rpc_controller->Failed())
    vw_throw(IOErr() << "OpenRequest RPC failed: " << m_rpc_controller->ErrorText());

  vw_out(0) << "RECEIVED IndexReadReply\n" << response.DebugString() << "\n";

  return response.index_record();
}
  
// Writing, pt. 1: Locks a blob and returns the blob id that can
// be used to write a tile.
int vw::platefile::RemoteIndex::write_request(int size) {
  vw_throw(NoImplErr() << "write_request() not yet implemented.");
  return 0;

  // IndexWriteRequest request;
  // request.set_platefile_id(m_platefile_id);
  // request.set_secret(m_secret);
  // request.set_size(size);

  // IndexWriteReply response;
  // m_index_service->WriteRequest(m_rpc_controller.get(), request, response, 
  //                              protobuf::NewCallback(&null_closure));
  // if (m_rpc_controller.failed())
  //   vw_throw(IOErr() << "OpenRequest RPC failed: " << m_rpc_controller.ErrorText());

  // return r.blob_id();
}
  
// Writing, pt. 2: Supply information to update the index and
// unlock the blob id.
void vw::platefile::RemoteIndex::write_complete(TileHeader const& header, IndexRecord const& record) {
  vw_throw(NoImplErr() << "write_complete() not yet implemented.");


  // IndexWriteComplete request;
  // request.set_platefile_id(m_platefile_id);
  // request.set_secret(m_secret);
  // *(request.mutable_header()) = header;
  // *(request.mutable_record()) = record;
  // m_conn.basic_publish_protobuf(request, INDEX_EXCHANGE, "index.write_complete");    

  // IndexSuccess r = wait_for_response<IndexSuccess>(m_queue_name + ".write_complete");
}
  
  
vw::int32 vw::platefile::RemoteIndex::version() const { 
  vw_throw(NoImplErr() << "version() not yet implemented.");
  return 0;
}

vw::int32 vw::platefile::RemoteIndex::max_depth() const { 
  vw_throw(NoImplErr() << "max_depth() not yet implemented.");
  return 0;

  // IndexDepthRequest request;
  // request.set_platefile_id(m_platefile_id);
  // request.set_secret(m_secret);
  // m_conn.basic_publish_protobuf(request, INDEX_EXCHANGE, m_queue_name);    

  // IndexDepthReply r = wait_for_response<IndexDepthReply>(m_queue_name + ".depth_reply");
  // return r.depth();
}

std::string vw::platefile::RemoteIndex::platefile_name(int platefile_id) const { 
  IndexInfoRequest request;
  vw_out(0) << "Setting to ID = " << platefile_id << "\n";
   request.set_platefile_id(platefile_id);
   vw_out(0) << "ISSUING IndexInfoRequest\n" << request.DebugString() << "\n";

  IndexInfoReply response;
  m_index_service->InfoRequest(m_rpc_controller.get(), &request, &response, 
                               google::protobuf::NewCallback(&null_closure));

  vw_out(0) << "RECEIVED IndexInfoResponse\n" << response.DebugString() << "\n";

  if (m_rpc_controller->Failed())
    vw_throw(IOErr() << "InfoRequest RPC failed: " << m_rpc_controller->ErrorText());
  return response.plate_filename();
}

vw::int32 vw::platefile::RemoteIndex::default_tile_size() const { 
  vw_throw(NoImplErr() << "default_tile_size() not yet implemented.");
  return 0;
}

std::string vw::platefile::RemoteIndex::default_tile_filetype() const { 
  vw_throw(NoImplErr() << "default_tile_filetype() not yet implemented.");
  return "";
}

vw::PixelFormatEnum vw::platefile::RemoteIndex::pixel_format() const {
  vw_throw(NoImplErr() << "pixel_format() not yet implemented.");
  return VW_PIXEL_GRAY;
}

vw::ChannelTypeEnum vw::platefile::RemoteIndex::channel_type() const {
  vw_throw(NoImplErr() << "channel_type() not yet implemented.");
  return VW_CHANNEL_UINT8;
}


// --------------------- TRANSACTIONS ------------------------

// Clients are expected to make a transaction request whenever
// they start a self-contained chunk of mosaicking work.  .
vw::int32 vw::platefile::RemoteIndex::transaction_request(std::string transaction_description) {
  vw_throw(NoImplErr() << "transaction_request() not yet implemented.");
  return 0;
}
// Once a chunk of work is complete, clients can "commit" their
// work to the mosaic by issuding a transaction_complete method.
void vw::platefile::RemoteIndex::transaction_complete(int32 transaction_id) {
  vw_throw(NoImplErr() << "transaction_complete() not yet implemented.");
}

vw::int32 vw::platefile::RemoteIndex::transaction_cursor() {
  vw_throw(NoImplErr() << "transaction_cursor() not yet implemented.");
  return 0;
}


