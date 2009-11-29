// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Core/Exception.h>
#include <vw/Plate/RemoteIndex.h>
#include <vw/Plate/common.h>
#include <vw/Plate/ProtoBuffers.pb.h>
using namespace vw;
using namespace vw::platefile;

#include <unistd.h>
#include <boost/algorithm/string/split.hpp>
 
// A dummy method for passing to the RPC calls below.
static void null_closure() {}

// Parse a URL with the format: pf://<exchange>/<platefile name>.plate
void parse_url(std::string const& url, std::string &exchange, std::string &platefile_name) {
    if (url.find("pf://") != 0) {
      vw_throw(vw::ArgumentErr() << "RemoteIndex::parse_url() -- this does not appear to be a well-formed URL: " << url);
    } else {
      std::string substr = url.substr(5);

      std::vector<std::string> split_vec;
      boost::split( split_vec, substr, boost::is_any_of("/") );

      if (split_vec.size() != 2) 
        vw_throw(vw::ArgumentErr() << "RemoteIndex::parse_url() -- could not parse exchange and platefile_name from URL string: " << url);
      
      exchange = split_vec[0];
      platefile_name = split_vec[1];
    }

}

/// Constructor (for opening an existing Index)
vw::platefile::RemoteIndex::RemoteIndex(std::string const& url) {

  // Parse the URL string into a separate 'exchange' and 'platefile_name' field.
  std::string routing_key;
  std::string platefile_name;
  parse_url(url, routing_key, platefile_name);

  m_queue_name = AmqpRpcClient::UniqueQueueName(std::string("remote_index_") + platefile_name);

  // Set up the connection to the AmqpRpcService
  m_rpc_channel.reset( new AmqpRpcChannel(INDEX_EXCHANGE, routing_key, m_queue_name) );
  m_rpc_controller.reset ( new AmqpRpcClient() );
  m_index_service.reset ( new IndexService::Stub(m_rpc_channel.get()) );
  
  // Send an IndexOpenRequest to the AMQP index server.
  IndexOpenRequest request;
  request.set_plate_name(platefile_name);

  IndexOpenReply response;
  m_index_service->OpenRequest(m_rpc_controller.get(), &request, &response, 
                               google::protobuf::NewCallback(&null_closure));

  m_index_header = response.index_header();
  m_platefile_id = m_index_header.platefile_id();
  m_short_plate_filename = response.short_plate_filename();
  m_full_plate_filename = response.full_plate_filename();
  vw_out(InfoMessage, "plate") << "Opened remote platefile \"" << platefile_name
                               << "\"   ID: " << m_platefile_id << "\n";
}

/// Constructor (for creating a new Index)
vw::platefile::RemoteIndex::RemoteIndex(std::string const& url, IndexHeader index_header_info) :
  m_index_header(index_header_info) {

  // Parse the URL string into a separate 'exchange' and 'platefile_name' field.
  std::string routing_key;
  std::string platefile_name;
  parse_url(url, routing_key, platefile_name);

  m_queue_name = AmqpRpcClient::UniqueQueueName(std::string("remote_index_") + platefile_name);

  // Set up the connection to the AmqpRpcService
  m_rpc_channel.reset( new AmqpRpcChannel(INDEX_EXCHANGE, routing_key, m_queue_name) );
  m_rpc_controller.reset ( new AmqpRpcClient() );
  m_index_service.reset ( new IndexService::Stub(m_rpc_channel.get()) );
  
  // Send an IndexCreateRequest to the AMQP index server.
  IndexCreateRequest request;
  request.set_plate_name(platefile_name);
  index_header_info.set_platefile_id(0);  // this takes care of 'required' property, which is not set yet
  *(request.mutable_index_header()) = index_header_info;

  vw_out(0) << "ISSUING IndexCreateRequest\n" << request.DebugString() << "\n";

  IndexOpenReply response;
  m_index_service->CreateRequest(m_rpc_controller.get(), &request, &response, 
                                 google::protobuf::NewCallback(&null_closure));

  vw_out(0) << "RECEIVED IndexCreateReply\n" << response.DebugString() << "\n";

  m_index_header = response.index_header();
  m_platefile_id = m_index_header.platefile_id();
  m_short_plate_filename = response.short_plate_filename();
  m_full_plate_filename = response.full_plate_filename();
  vw_out(InfoMessage, "plate") << "Created remote platefile \"" << platefile_name
                               << "\"   ID: " << m_platefile_id << "\n";
}


/// Destructor
vw::platefile::RemoteIndex::~RemoteIndex() {}
  

/// Attempt to access a tile in the index.  Throws an
/// TileNotFoundErr if the tile cannot be found.
vw::platefile::IndexRecord vw::platefile::RemoteIndex::read_request(int col, int row, int depth, 
                                                                    int transaction_id) {
  IndexReadRequest request;
  request.set_platefile_id(m_platefile_id);
  request.set_col(col);
  request.set_row(row);
  request.set_depth(depth);
  request.set_transaction_id(transaction_id);

  IndexReadReply response;
  m_index_service->ReadRequest(m_rpc_controller.get(), &request, &response, 
                               google::protobuf::NewCallback(&null_closure));
  return response.index_record();
}
  
// Writing, pt. 1: Locks a blob and returns the blob id that can
// be used to write a tile.
int vw::platefile::RemoteIndex::write_request(int size) {
  IndexWriteRequest request;
  request.set_platefile_id(m_platefile_id);
  request.set_size(size);

  IndexWriteReply response;
  m_index_service->WriteRequest(m_rpc_controller.get(), &request, &response, 
                                google::protobuf::NewCallback(&null_closure));
  return response.blob_id();
}
  
// Writing, pt. 2: Supply information to update the index and
// unlock the blob id.
void vw::platefile::RemoteIndex::write_complete(TileHeader const& header, IndexRecord const& record) {
  IndexWriteComplete request;
  request.set_platefile_id(m_platefile_id);
  *(request.mutable_header()) = header;
  *(request.mutable_record()) = record;

  RpcNullMessage response;
  m_index_service->WriteComplete(m_rpc_controller.get(), &request, &response, 
                                 google::protobuf::NewCallback(&null_closure));
}
  
vw::int32 vw::platefile::RemoteIndex::max_depth() const { 
  IndexDepthRequest request;
  request.set_platefile_id(m_platefile_id);
  
  IndexDepthReply response;
  m_index_service->DepthRequest(m_rpc_controller.get(), &request, &response, 
                               google::protobuf::NewCallback(&null_closure));
  return response.depth();
}

vw::int32 vw::platefile::RemoteIndex::version() const { 
  return m_index_header.version();
}

std::string vw::platefile::RemoteIndex::platefile_name() const { 
  std::cout << "** Call to platefile_name() : " << m_full_plate_filename << "\n";
  return m_full_plate_filename;
}

IndexHeader RemoteIndex::index_header() const { 
  return m_index_header;
}

vw::int32 vw::platefile::RemoteIndex::tile_size() const { 
  return m_index_header.tile_size();
}

std::string vw::platefile::RemoteIndex::tile_filetype() const { 
  return m_index_header.tile_filetype();
}

vw::PixelFormatEnum vw::platefile::RemoteIndex::pixel_format() const {
  return vw::PixelFormatEnum(m_index_header.pixel_format());
}

vw::ChannelTypeEnum vw::platefile::RemoteIndex::channel_type() const {
  return vw::ChannelTypeEnum(m_index_header.channel_type());
}

// --------------------- TRANSACTIONS ------------------------

// Clients are expected to make a transaction request whenever
// they start a self-contained chunk of mosaicking work.  .
vw::int32 vw::platefile::RemoteIndex::transaction_request(std::string transaction_description) {
  IndexTransactionRequest request;
  request.set_platefile_id(m_platefile_id);
  request.set_description(transaction_description);
  
  IndexTransactionReply response;
  m_index_service->TransactionRequest(m_rpc_controller.get(), &request, &response, 
                                      google::protobuf::NewCallback(&null_closure));
  return response.transaction_id();
}

// Once a chunk of work is complete, clients can "commit" their
// work to the mosaic by issuding a transaction_complete method.
void vw::platefile::RemoteIndex::transaction_complete(int32 transaction_id) {
  IndexTransactionComplete request;
  request.set_platefile_id(m_platefile_id);
  request.set_transaction_id(transaction_id);
  
  RpcNullMessage response;
  m_index_service->TransactionComplete(m_rpc_controller.get(), &request, &response, 
                                       google::protobuf::NewCallback(&null_closure));
}

vw::int32 vw::platefile::RemoteIndex::transaction_cursor() {
  IndexTransactionCursorRequest request;
  request.set_platefile_id(m_platefile_id);
  
  IndexTransactionCursorReply response;
  m_index_service->TransactionCursor(m_rpc_controller.get(), &request, &response, 
                                     google::protobuf::NewCallback(&null_closure));
  return response.transaction_id();
}


