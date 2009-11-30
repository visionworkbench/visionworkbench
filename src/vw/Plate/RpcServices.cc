// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Core/Exception.h>
#include <vw/Plate/Exception.h>
#include <vw/Plate/RpcServices.h>
#include <vw/Plate/ProtoBuffers.pb.h>
#include <vw/Core/Stopwatch.h>
using namespace vw;

#include <google/protobuf/descriptor.h>

// A dummy method for passing to the RPC calls below.
static void null_closure() {}

// -----------------------------------------------------------------------
//                           AmqpRpcChannel
// -----------------------------------------------------------------------
  
vw::platefile::AmqpRpcChannel::AmqpRpcChannel(std::string const& exchange, 
                                              std::string const& request_routing_key, 
                                              std::string const& response_queue) : 
  m_exchange(exchange),
  m_request_routing_key(request_routing_key),
  m_response_queue(response_queue) {
  
  // Create a queue and bind it to the index server exchange.
  m_conn.exchange_declare(exchange, "direct", true, false);
  m_conn.queue_declare(response_queue, true, true, false);
  m_conn.queue_bind(response_queue, exchange, response_queue);
}

// Client side RPC implementation using AMQP as the message
// passing interface.
void vw::platefile::AmqpRpcChannel::CallMethod(const google::protobuf::MethodDescriptor* method,
                                               google::protobuf::RpcController* controller,
                                               const google::protobuf::Message* request,
                                               google::protobuf::Message* response,
                                               google::protobuf::Closure* done) {

  Mutex::Lock lock(m_mutex);
    
  //  vw_out(0) << "[RPC --> " << method->name() << "]\n";
  //  vw_out(0) << "Request:\n" << request->DebugString() << "\n";

  // Serialize the message and pass it to AMQP to be transferred.
  // The WireMessege class saves the size of the message at the
  // beginning so that we know how many byets to parse through
  // when it reaches the other side.
  RpcRequestWrapper request_wrapper;
  request_wrapper.set_requestor(m_response_queue);
  request_wrapper.set_method(method->name());
  request_wrapper.set_payload(request->SerializeAsString());

  WireMessage wire_request(&request_wrapper);
  m_conn.basic_publish(wire_request.serialized_bytes(), 
                       wire_request.size(),
                       m_exchange, m_request_routing_key);
  
  // Wait for a response, and pass it along to the callback.
  std::string response_routing_key;
  boost::shared_array<uint8> response_bytes = m_conn.basic_consume(m_response_queue, 
                                                                   response_routing_key, 
                                                                   false); 

  WireMessage wire_response(response_bytes);
  RpcResponseWrapper response_wrapper;
  try {
    response_wrapper = wire_response.parse_as_message<RpcResponseWrapper>();
  } catch (RpcErr &e) {
    std::cout << "WARNING: An RPC error occurred: " << e.what() << "\n";
  }


  // Handle errors and exceptions. 
  //
  // I'll admit that this code is a little messy, and leaves something
  // to be desired.  Ideally we would write a general purpose method
  // that can regenerate an exception based on its name and
  // description strings so that we don't have to add additional code
  // here whenever we create a new PlatefileErr subclass.
  //
  if (response_wrapper.error()) {

    // For debugging:
    // vw_out(0) << "[RPC ERROR]  Type = " << response_wrapper.error_info().type()
    //           << "  Description = " << response_wrapper.error_info().message() << "\n";
    
    if (response_wrapper.error_info().type() == "TileNotFoundErr") {
      vw_throw(TileNotFoundErr() << response_wrapper.error_info().message());

    } else if (response_wrapper.error_info().type() == "InvalidPlatefileErr") {
      vw_throw(InvalidPlatefileErr() << response_wrapper.error_info().message());      

    } else if (response_wrapper.error_info().type() == "PlatefileCreationErr") {
      vw_throw(PlatefileCreationErr() << response_wrapper.error_info().message());      

    } else {
      vw_out(0) << "WARNING!! Unknown exception in RPC message:\n\t" 
                << "Type = " << response_wrapper.error_info().type() << "\n\t"
                << "Description = " << response_wrapper.error_info().message() << "\n";

      vw_throw(RpcErr() << "RPC error. Remote system threw this error: " 
               << response_wrapper.error_info().type() << " with this description: " 
               << response_wrapper.error_info().message());
    }

  } else {
    response->ParseFromString(response_wrapper.payload());
    //    vw_out(0) << "Response:\n" << response->DebugString() << "\n\n";
  }
  done->Run();
}

// -----------------------------------------------------------------------------
//                                AmqpRpcServer
// -----------------------------------------------------------------------------

void vw::platefile::AmqpRpcServer::run() {

  while(1) {
    // --------------------------------------
    // Step 1 : Wait for an incoming message.
    // --------------------------------------
    std::string routing_key;
    boost::shared_array<uint8> request_bytes = m_conn.basic_consume(m_queue, routing_key, false);
    WireMessage wire_request(request_bytes);
    RpcRequestWrapper request_wrapper = wire_request.parse_as_message<RpcRequestWrapper>();

    vw_out(0) << "[RPC: " << request_wrapper.method() 
              << " from " << request_wrapper.requestor() << "]\n";

    // -------------------------------------------------------------
    // Step 2 : Instantiate the proper messages and delegate them to
    // the proper method on the service.
    // -------------------------------------------------------------
    const google::protobuf::MethodDescriptor* method = 
      m_service->GetDescriptor()->FindMethodByName(request_wrapper.method());
      
    boost::shared_ptr<google::protobuf::Message> 
      request(m_service->GetRequestPrototype(method).New());
    boost::shared_ptr<google::protobuf::Message> 
      response(m_service->GetResponsePrototype(method).New());

    try {

      // Attempt to parse the actual request message from the
      // request_wrapper.
      if (!request->ParseFromString(request_wrapper.payload()))
        vw_throw(IOErr() << "Error parsing request from request_wrapper message.\n");

      // For debugging: 
      //      vw_out(0) << "Request:\n" << request->DebugString() << "\n";

      m_service->CallMethod(method, this, request.get(), response.get(), 
                            google::protobuf::NewCallback(&null_closure));

      // For debugging: 
      //      vw_out(0) << "Response:\n" << response->DebugString() << "\n\n";

      // ---------------------------
      // Step 3 : Return the result.
      // ---------------------------
      RpcResponseWrapper response_wrapper;
      response_wrapper.set_error(false);
      response_wrapper.set_payload(response->SerializeAsString());

      WireMessage wire_response(&response_wrapper);
      m_conn.basic_publish(wire_response.serialized_bytes(), 
                           wire_response.size(),
                           m_exchange, request_wrapper.requestor() );

    } catch (PlatefileErr &e) {

      RpcResponseWrapper response_wrapper;

      // If an exception occurred on the index_server, we pass it
      // along to the requestor by setting the RpcErrorInfo portion of
      // the RpcResponseWrapper.
      RpcErrorInfo err;
      err.set_type(e.name());
      err.set_message(e.what());
      response_wrapper.set_error(true);
      *(response_wrapper.mutable_error_info()) = err;

      WireMessage wire_response(&response_wrapper);
      m_conn.basic_publish(wire_response.serialized_bytes(), 
                           wire_response.size(),
                           m_exchange, request_wrapper.requestor() );

    } catch (vw::Exception &e) {

      vw_out(0) << "WARNING!! Uncaught Vision Workbench Exception:\n\t" << e.what() << "\n";

    }

  }    
}

std::string vw::platefile::AmqpRpcClient::UniqueQueueName(const std::string identifier) {
  // Start by generating a unique queue name based on our hostname, PID, and thread ID.
  char hostname[255];
  gethostname(hostname, 255);
  std::ostringstream requestor;
  requestor << identifier << "_" << hostname << "_" << getpid() << "_" << Thread::id() << vw::Stopwatch::microtime(false);
  return requestor.str();
}
