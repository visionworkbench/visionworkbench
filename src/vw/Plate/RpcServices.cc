// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/RpcServices.h>
#include <vw/Plate/ProtoBuffers.pb.h>

#include <google/protobuf/descriptor.h>

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
