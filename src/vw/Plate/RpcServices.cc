// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Core/Exception.h>
#include <vw/Plate/RpcServices.h>
#include <vw/Plate/ProtoBuffers.pb.h>

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
    
  vw_out(0) << "[RPC] : " << method->name() << "\n";
  vw_out(0) << "Request:\n" << request->DebugString() << "\n";

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
  RpcResponseWrapper response_wrapper = wire_response.parse_as_message<RpcResponseWrapper>();
  if (response_wrapper.error()) {
    vw_out(0) << "WARNING: an RPC error occured.  Type = " 
              << response_wrapper.error_info().type()
              << "  Description = " << response_wrapper.error_info().message() << "\n";
    controller->SetFailed(response_wrapper.error_info().message());
  }
  response->ParseFromString(response_wrapper.payload());
  vw_out(0) << "Response:\n" << response->DebugString() << "\n\n";

  done->Run();
}

// -----------------------------------------------------------------------------
//                                AmqpRpcServer
// -----------------------------------------------------------------------------

void vw::platefile::AmqpRpcServer::run() {

  while(1) {
    
    // Step 1 : Wait for an incoming message.
    std::string routing_key;
    std::cout << "\n\nWaiting for message...\n";
    boost::shared_array<uint8> request_bytes = m_conn.basic_consume(m_queue, 
                                                                    routing_key, 
                                                                    false);
    WireMessage wire_request(request_bytes);
    RpcRequestWrapper request_wrapper = wire_request.parse_as_message<RpcRequestWrapper>();
    std::cout << "[RPC: " << request_wrapper.method() 
              << " from " << request_wrapper.requestor() << "]\n";

    // Step 2 : Delegate the message to the proper method on the service.
    const google::protobuf::MethodDescriptor* method = 
      m_service->GetDescriptor()->FindMethodByName(request_wrapper.method());
      
    boost::shared_ptr<google::protobuf::Message> 
      request(m_service->GetRequestPrototype(method).New());
      
    boost::shared_ptr<google::protobuf::Message> 
      response(m_service->GetResponsePrototype(method).New());
      
    try {

      request->ParseFromString(request_wrapper.payload());
      std::cout << "Request:\n" << request->DebugString() << "\n";
      m_service->CallMethod(method, this, request.get(), response.get(), 
                            google::protobuf::NewCallback(&null_closure));
      std::cout << "Response:\n" << response->DebugString() << "\n\n";
        
      // Step 3 : Return the result.
      RpcResponseWrapper response_wrapper;
      response_wrapper.set_payload(response->SerializeAsString());
      response_wrapper.set_error(false);

      // If the RPC generated an error, we pass it along here. 
      if (this->Failed()) {
        response_wrapper.set_error(true);
        RpcErrorInfo msg;
        msg.set_type("Rpc Error");
        msg.set_message(this->ErrorText());
        this->Reset();
      }
        
      WireMessage wire_response(&response_wrapper);
      m_conn.basic_publish(wire_response.serialized_bytes(), 
                           wire_response.size(),
                           m_exchange, request_wrapper.requestor() );
      

    } catch (vw::Exception &e) {

      RpcResponseWrapper response_wrapper;
      response_wrapper.set_payload("");

      // If an exception occurred on the plateindex_server, we pass
      // it along to the requestor as well.
      RpcErrorInfo msg;
      msg.set_type(e.name());
      msg.set_message(e.what());
      response_wrapper.set_error(true);
      *(response_wrapper.mutable_error_info()) = msg;

      WireMessage wire_response(&response_wrapper);
      m_conn.basic_publish(wire_response.serialized_bytes(), 
                           wire_response.size(),
                           m_exchange, request_wrapper.requestor() );
    } 

  }    
}
