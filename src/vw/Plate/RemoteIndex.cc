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
#include <boost/lexical_cast.hpp>
 
// A dummy method for passing to the RPC calls below.
static void null_closure() {}

// Parse a URL with the format: pf://<exchange>/<platefile name>.plate
void parse_url(std::string const& url, std::string &hostname, int &port, 
               std::string &exchange, std::string &platefile_name) {
    if (url.find("pf://") != 0) {
      vw_throw(vw::ArgumentErr() << "RemoteIndex::parse_url() -- this does not appear to be a well-formed URL: " << url);
    } else {
      std::string substr = url.substr(5);

      std::vector<std::string> split_vec;
      boost::split( split_vec, substr, boost::is_any_of("/") );

      // No hostname was specified: pf://<routing_key>/<platefilename>.plate
      if (split_vec.size() == 2) {
        hostname = "localhost"; // default to localhost
        port = 5672;            // default rabbitmq port

        exchange = split_vec[0];
        platefile_name = split_vec[1];

      // No hostname was specified: pf://<ip address>:<port>/<routing_key>/<platefilename>.plate
      } else if (split_vec.size() == 3) {

        exchange = split_vec[1];
        platefile_name = split_vec[2];

        std::vector<std::string> host_port_vec;
        boost::split( host_port_vec, split_vec[0], boost::is_any_of(":") );
        
        if (host_port_vec.size() == 1) {

          hostname = host_port_vec[0];
          port = 5672;            // default rabbitmq port

        } else if (host_port_vec.size() == 2) {

          hostname = host_port_vec[0];
          port = boost::lexical_cast<int>(host_port_vec[1]);

        } else {
          vw_throw(vw::ArgumentErr() << "RemoteIndex::parse_url() -- " 
                   << "could not parse hostname and port from URL string: " << url);
        }

      } else {
        vw_throw(vw::ArgumentErr() << "RemoteIndex::parse_url() -- "
                 << "could not parse URL string: " << url);
      }
    }

}

/// Constructor (for opening an existing Index)
vw::platefile::RemoteIndex::RemoteIndex(std::string const& url) {

  // Parse the URL string into a separate 'exchange' and 'platefile_name' field.
  std::string hostname;
  int port;
  std::string routing_key;
  std::string platefile_name;
  parse_url(url, hostname, port, routing_key, platefile_name);

  std::string queue_name = AmqpRpcClient::UniqueQueueName(std::string("remote_index_") + platefile_name);

  // Set up the connection to the AmqpRpcService
  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection(hostname, port));
  m_rpc_controller.reset( new AmqpRpcClient(conn, INDEX_EXCHANGE, queue_name, routing_key) );
  m_index_service.reset ( new IndexService::Stub(m_rpc_controller.get() ) );
  m_rpc_controller->bind_service(m_index_service, queue_name);

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
  std::string hostname;
  int port;
  std::string routing_key;
  std::string platefile_name;
  parse_url(url, hostname, port, routing_key, platefile_name);

  std::string queue_name = AmqpRpcClient::UniqueQueueName(std::string("remote_index_") + platefile_name);

  // Set up the connection to the AmqpRpcService
  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection(hostname, port));
  m_rpc_controller.reset( new AmqpRpcClient(conn, INDEX_EXCHANGE, queue_name, routing_key) );
  m_index_service.reset ( new IndexService::Stub(m_rpc_controller.get() ) );
  m_rpc_controller->bind_service(m_index_service, queue_name);

  // Send an IndexCreateRequest to the AMQP index server.
  IndexCreateRequest request;
  request.set_plate_name(platefile_name);
  index_header_info.set_platefile_id(0);  // this takes care of this 'required' 
                                          // protobuf property, which is not set yet
  *(request.mutable_index_header()) = index_header_info;

  IndexOpenReply response;
  m_index_service->CreateRequest(m_rpc_controller.get(), &request, &response, 
                                 google::protobuf::NewCallback(&null_closure));

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
vw::platefile::IndexRecord vw::platefile::RemoteIndex::read_request(int col, int row, int level, 
                                                                    int transaction_id, bool exact_transaction_match) {
  IndexReadRequest request;
  request.set_platefile_id(m_platefile_id);
  request.set_col(col);
  request.set_row(row);
  request.set_level(level);
  request.set_transaction_id(transaction_id);
  request.set_exact_transaction_match(exact_transaction_match);

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
  
vw::int32 vw::platefile::RemoteIndex::num_levels() const { 
  IndexNumLevelsRequest request;
  request.set_platefile_id(m_platefile_id);
  
  IndexNumLevelsReply response;
  m_index_service->NumLevelsRequest(m_rpc_controller.get(), &request, &response, 
                                    google::protobuf::NewCallback(&null_closure));
  return response.num_levels();
}

vw::int32 vw::platefile::RemoteIndex::version() const { 
  return m_index_header.version();
}

std::string vw::platefile::RemoteIndex::platefile_name() const { 
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
vw::int32 vw::platefile::RemoteIndex::transaction_request(std::string transaction_description,
                                                          std::vector<TileHeader> const& tile_headers) {
  IndexTransactionRequest request;
  request.set_platefile_id(m_platefile_id);
  request.set_description(transaction_description);
  for (size_t i = 0; i < tile_headers.size(); ++i)
    *(request.mutable_tile_headers()->Add()) = tile_headers[i];
  
  IndexTransactionReply response;
  m_index_service->TransactionRequest(m_rpc_controller.get(), &request, &response, 
                                      google::protobuf::NewCallback(&null_closure));
  return response.transaction_id();
}

/// Called right before the beginning of the mipmapping pass
void vw::platefile::RemoteIndex::root_complete(int32 transaction_id,
                                               std::vector<TileHeader> const& tile_headers) {
  IndexRootComplete request;
  request.set_platefile_id(m_platefile_id);
  request.set_transaction_id(transaction_id);
  for (size_t i = 0; i < tile_headers.size(); ++i)
    *(request.mutable_tile_headers()->Add()) = tile_headers[i];
  
  RpcNullMessage response;
  m_index_service->RootComplete(m_rpc_controller.get(), &request, &response, 
                                google::protobuf::NewCallback(&null_closure));
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

// If a transaction fails, we may need to clean up the mosaic.  
void vw::platefile::RemoteIndex::transaction_failed(int32 transaction_id) {
  IndexTransactionFailed request;
  request.set_platefile_id(m_platefile_id);
  request.set_transaction_id(transaction_id);
  
  RpcNullMessage response;
  m_index_service->TransactionFailed(m_rpc_controller.get(), &request, &response, 
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


