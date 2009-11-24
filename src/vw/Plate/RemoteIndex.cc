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
void null_closure() {}

// Parse a URL with the format: pf://<exchange>/<platefile name>.plate
void parse_url(std::string const& url, std::string &exchange, std::string &platefile_name) {
    if (url.find("pf://") != 0) {
      vw_throw(vw::ArgumentErr() << "RemoteIndex::parse_url() -- this does not appear to be a well-formed URL: " << url);
    } else {
      std::string substr = url.substr(5);
      std::cout << "REMOTE INDEX --> Parsing URL: " << substr << "\n";

      std::vector<std::string> split_vec;
      boost::split( split_vec, substr, boost::is_any_of("/") );

      if (split_vec.size() != 2) 
        vw_throw(vw::ArgumentErr() << "RemoteIndex::parse_url() -- could not parse exchange and platefile_name from URL string: " << url);
      
      exchange = split_vec[0];
      platefile_name = split_vec[1];

      std::cout << "REMOTE INDEX --> Opened " << exchange << " " << platefile_name << "\n";
    }

}

/// Constructor (for opening an existing Index)
vw::platefile::RemoteIndex::RemoteIndex(std::string const& url) {

  // Parse the URL string into a separate 'exchange' and 'platefile_name' field.
  std::string exchange;
  std::string platefile_name;
  parse_url(url, exchange, platefile_name);

  std::cout << "Parsed URL: " << exchange << " " << platefile_name << "\n";

  // Start by generating a unique queue name based on our hostname, PID, and thread ID.
  char hostname[255];
  gethostname(hostname, 255);
  std::ostringstream requestor;
  requestor << "remote_index_" << hostname << "_" << getpid() << "_" << Thread::id();
  m_queue_name = requestor.str();

  // Set up the connection to the AmqpRpcService
  m_rpc_channel.reset( new AmqpRpcChannel(INDEX_EXCHANGE, exchange, requestor.str()) );
  m_rpc_controller.reset ( new AmqpRpcController() );
  m_index_service.reset ( new IndexService::Stub(m_rpc_channel.get()) );
  
  // Send an IndexOpenRequest to the AMQP index server.
  IndexOpenRequest request;
  request.set_plate_name(platefile_name);

  vw_out(0) << "ISSUING IndexOpenRequest\n" << request.DebugString() << "\n";

  IndexOpenReply response;
  m_index_service->OpenRequest(m_rpc_controller.get(), &request, &response, 
                               google::protobuf::NewCallback(&null_closure));

  if (m_rpc_controller->Failed())
    vw_throw(IOErr() << "OpenRequest RPC failed: " << m_rpc_controller->ErrorText());

  vw_out(0) << "RECEIVED IndexOpenReply\n" << response.DebugString() << "\n";

  m_index_header = response.index_header();
  m_platefile_id = m_index_header.platefile_id();
  vw_out(InfoMessage, "plate") << "Opened remote platefile \"" << platefile_name
                               << " on exchange " << exchange
                               << "\"   ID: " << m_platefile_id << "\n";
}

/// Constructor (for creating a new Index)
vw::platefile::RemoteIndex::RemoteIndex(std::string const& url, IndexHeader index_header_info) :
  m_index_header(index_header_info) {

  // Parse the URL string into a separate 'exchange' and 'platefile_name' field.
  std::string exchange;
  std::string platefile_name;
  parse_url(url, exchange, platefile_name);

  // Start by generating a unique queue name
  char hostname[255];
  gethostname(hostname, 255);
  std::ostringstream requestor;
  requestor << "remote_index_" << hostname << "_" << getpid() << "_" << Thread::id() << "\n";
  m_queue_name = requestor.str();

  // Set up the connection to the AmqpRpcService
  m_rpc_channel.reset( new AmqpRpcChannel(exchange, "index", requestor.str()) );
  m_rpc_controller.reset ( new AmqpRpcController() );
  m_index_service.reset ( new IndexService::Stub(m_rpc_channel.get()) );
  
  // Send an IndexCreateRequest to the AMQP index server.
  IndexCreateRequest request;
  *(request.mutable_index_header()) = index_header_info;

  vw_out(0) << "ISSUING IndexCreateRequest\n" << request.DebugString() << "\n";

  IndexOpenReply response;
  m_index_service->CreateRequest(m_rpc_controller.get(), &request, &response, 
                                 google::protobuf::NewCallback(&null_closure));

  if (m_rpc_controller->Failed())
    vw_throw(IOErr() << "OpenRequest RPC failed: " << m_rpc_controller->ErrorText());

  vw_out(0) << "RECEIVED IndexCreateReply\n" << response.DebugString() << "\n";

  m_index_header = response.index_header();
  m_platefile_id = m_index_header.platefile_id();
  vw_out(InfoMessage, "plate") << "Opened remote platefile \"" << platefile_name
                               << " on exchange " << exchange
                               << "\"   ID: " << m_platefile_id << "\n";
}


/// Destructor
vw::platefile::RemoteIndex::~RemoteIndex() {}
  

/// Attempt to access a tile in the index.  Throws an
/// TileNotFoundErr if the tile cannot be found.
vw::platefile::IndexRecord vw::platefile::RemoteIndex::read_request(int col, int row, 
                                                                    int depth, 
                                                                    int transaction_id) {
  IndexReadRequest request;
  request.set_platefile_id(m_platefile_id);
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
  // *(request.mutable_header()) = header;
  // *(request.mutable_record()) = record;
  // m_conn.basic_publish_protobuf(request, INDEX_EXCHANGE, "index.write_complete");    

  // IndexSuccess r = wait_for_response<IndexSuccess>(m_queue_name + ".write_complete");
}
  
  
vw::int32 vw::platefile::RemoteIndex::version() const { 
  return m_index_header.version();
}

vw::int32 vw::platefile::RemoteIndex::max_depth() const { 
  IndexDepthRequest request;
  request.set_platefile_id(m_platefile_id);
  vw_out(0) << "[IndexDepthRequest] :\n" << request.DebugString() << "\n";
  
  IndexDepthReply response;
  m_index_service->DepthRequest(m_rpc_controller.get(), &request, &response, 
                               google::protobuf::NewCallback(&null_closure));
  
  vw_out(0) << "[IndexDepthResponse] :\n" << response.DebugString() << "\n";
  
  if (m_rpc_controller->Failed())
    vw_throw(IOErr() << "DepthRequest RPC failed: " << m_rpc_controller->ErrorText());
  return response.depth();
}

std::string vw::platefile::RemoteIndex::platefile_name() const { 
  IndexInfoRequest request;
  request.set_platefile_id(m_platefile_id);
  vw_out(0) << "[IndexInfoRequest] :\n" << request.DebugString() << "\n";
  
  IndexInfoReply response;
  m_index_service->InfoRequest(m_rpc_controller.get(), &request, &response, 
                               google::protobuf::NewCallback(&null_closure));
  
  vw_out(0) << "[IndexInfoResponse] :\n" << response.DebugString() << "\n";
  
  if (m_rpc_controller->Failed())
    vw_throw(IOErr() << "InfoRequest RPC failed: " << m_rpc_controller->ErrorText());
  return response.full_plate_filename();
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


