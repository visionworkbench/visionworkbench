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
using namespace vw::platefile;

#include <google/protobuf/descriptor.h>

// A dummy method for passing to the RPC calls below.
static void null_closure() {}

// -----------------------------------------------------------------------------
//                                AmqpRpcServer
// -----------------------------------------------------------------------------

void vw::platefile::AmqpRpcServer::run() {

  while(!m_terminate) {

    try {

      // --------------------------------------
      // Step 1 : Wait for an incoming message.
      // --------------------------------------
      RpcRequestWrapper request_wrapper;
      try {
        this->get_message(request_wrapper, -1);
        vw_out(DebugMessage, "platefile::rpc") << "[RPC: " << request_wrapper.method() 
                                               << " from " << request_wrapper.requestor() 
                                               << "  SEQ: " << request_wrapper.sequence_number() << "]\n";
      } catch (const vw::platefile::RpcErr&e) {
        vw_out() << "Invalid RPC, ignoring." << std::endl;
        continue;
      }

      // Record statistics
      m_stats.record_query(request_wrapper.ByteSize());

      // -------------------------------------------------------------
      // Step 2 : Instantiate the proper messages and delegate them to
      // the proper method on the service.
      // -------------------------------------------------------------
      
      try {

        const google::protobuf::MethodDescriptor* method =
          m_service->GetDescriptor()->FindMethodByName(request_wrapper.method());

        if (method == NULL)
          vw_throw(RpcErr() << "Unrecognized RPC method: " << request_wrapper.method());

        boost::shared_ptr<google::protobuf::Message>
          request(m_service->GetRequestPrototype(method).New());
        boost::shared_ptr<google::protobuf::Message>
          response(m_service->GetResponsePrototype(method).New());

        // Attempt to parse the actual request message from the
        // request_wrapper.
        RpcResponseWrapper response_wrapper;
        if (!request->ParseFromString(request_wrapper.payload()))
          vw_throw(IOErr() << "Error parsing request from request_wrapper message.\n");

        // For debugging:
        //        std::cout << "Request: " << request->DebugString() << "\n";

        m_service->CallMethod(method, this, request.get(), response.get(),
                              google::protobuf::NewCallback(&null_closure));

        // ---------------------------
        // Step 3 : Return the result.
        // ---------------------------
        response_wrapper.set_sequence_number(request_wrapper.sequence_number());
        response_wrapper.set_error(false);
        response_wrapper.set_payload(response->SerializeAsString());

        // If the requestor field was blank, then we simply don't send
        // a response at all.  Otherwise, we do.
        if (request_wrapper.requestor() != "")
          this->send_message(response_wrapper, request_wrapper.requestor());

      } catch (PlatefileErr &e) {

        RpcResponseWrapper response_wrapper;
        // If an exception occurred on the index_server, we pass it
        // along to the requestor by setting the RpcErrorInfo portion of
        // the RpcResponseWrapper.
        RpcErrorInfo err;
        err.set_type(e.name());
        err.set_message(e.what());
        response_wrapper.set_sequence_number(request_wrapper.sequence_number());
        response_wrapper.set_error(true);
        *(response_wrapper.mutable_error_info()) = err;

        this->send_message(response_wrapper, request_wrapper.requestor());
      }
    } catch (const AMQPErr &e) {
      vw_out() << "WARNING!! Uncaught AMQP error:\n\t" << e.what() << "\n";
    } catch (const vw::Exception &e) {
      vw_out() << "WARNING!! Uncaught Vision Workbench Exception:\n\t" << e.what() << "\n";
    }

  } // while
}

// -----------------------------------------------------------------------------
//                                AmqpRpcClient
// -----------------------------------------------------------------------------

// Client side RPC implementation using AMQP as the message
// passing interface.
void vw::platefile::AmqpRpcClient::CallMethod(const google::protobuf::MethodDescriptor* method,
                                              google::protobuf::RpcController* controller,
                                              const google::protobuf::Message* request,
                                              google::protobuf::Message* response,
                                              google::protobuf::Closure* done) {

  AmqpRpcEndpoint* real_controller = dynamic_cast<AmqpRpcEndpoint*>(controller);
  if (!real_controller)
    vw_throw(LogicErr() << "AmqpRpcClient::CallMethod(): Unknown RpcController");
  
  // For debugging:
  //  std::cout << "Request: " << request->DebugString() << "\n";

  // Serialize the message and pass it to AMQP to be transferred.
  RpcRequestWrapper request_wrapper;
  request_wrapper.set_requestor(real_controller->queue_name());
  request_wrapper.set_method(method->name());
  request_wrapper.set_payload(request->SerializeAsString());

  RpcResponseWrapper response_wrapper;
  unsigned ntries = 0;
  bool success = false;

  unsigned request_seq;
  while (!success && ntries < m_max_tries) {

    // Set the sequence number and send the message.
    request_seq = ++m_sequence_number;
    request_wrapper.set_sequence_number(request_seq);

    // For debugging:
    vw_out(DebugMessage, "platefile::rpc") << "[ RPC: " << request_wrapper.method() 
                                           << " from " << request_wrapper.requestor() 
                                           << "  SEQ: " << request_wrapper.sequence_number() << " ]\n";

    // Send the message
    real_controller->send_message(request_wrapper, m_request_routing_key);

    // Try to get the message.  If we succeed, check the sequence number.
    try {
        real_controller->get_message(response_wrapper, m_timeout);

      while (response_wrapper.sequence_number() != request_seq) {
        vw_out(WarningMessage)
          << "CallMethod() sequence number did not match.  ("
          << request_seq << " != " << response_wrapper.sequence_number()
          << ").  Retrying...\n";

        if (++ntries >= m_max_tries)
          break;

        real_controller->get_message(response_wrapper, m_timeout);
      }

      // If the sequence number does appear, then we declare victory
      // and continue onward to parse the result.
      if (response_wrapper.sequence_number() == request_seq)
        success = true;

    } catch (AMQPTimeout &e) {
      // If we timed out, increment the number of tries and loop back
      // to the beginning.
      ++ntries;
      if (ntries < m_max_tries)
        vw_out(WarningMessage) << "CallMethod() timed out.  Executing retry #" << ntries << ".\n";
    }
  }

  if (!success)
    vw_throw(AMQPTimeout() << "CallMethod timed out completely after "
                           << m_max_tries << " tries.");


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
    // vw_out() << "[RPC ERROR]  Type = " << response_wrapper.error_info().type()
    //           << "  Description = " << response_wrapper.error_info().message() << "\n";

    if (response_wrapper.error_info().type() == "TileNotFoundErr") {
      vw_throw(TileNotFoundErr() << response_wrapper.error_info().message());

    } else if (response_wrapper.error_info().type() == "InvalidPlatefileErr") {
      vw_throw(InvalidPlatefileErr() << response_wrapper.error_info().message());

    } else if (response_wrapper.error_info().type() == "PlatefileCreationErr") {
      vw_throw(PlatefileCreationErr() << response_wrapper.error_info().message());

    } else if (response_wrapper.error_info().type() == "BlobLimitErr") {
      vw_throw(BlobLimitErr() << response_wrapper.error_info().message());

    } else {
      vw_out() << "WARNING!! Unknown exception in RPC message:\n\t"
                << "Type = " << response_wrapper.error_info().type() << "\n\t"
                << "Description = " << response_wrapper.error_info().message() << "\n";

      vw_throw(RpcErr() << "RPC error. Remote system threw this error: "
               << response_wrapper.error_info().type() << " with this description: "
               << response_wrapper.error_info().message());
    }

  } else {
    response->ParseFromString(response_wrapper.payload());

    // For debugging:
    //    vw_out() << "Response:\n" << response->DebugString() << "\n\n";
  }
  done->Run();
}

std::string vw::platefile::AmqpRpcClient::UniqueQueueName(const std::string identifier) {
  // Start by generating a unique queue name based on our hostname, PID, and thread ID.
  char hostname[255];
  gethostname(hostname, 255);
  std::ostringstream requestor;
  requestor << identifier << "_" << hostname << "_" << getpid() << "_" << Thread::id() << vw::Stopwatch::microtime(false);
  return requestor.str();
}

AmqpRpcEndpoint::AmqpRpcEndpoint(boost::shared_ptr<AmqpConnection> conn, std::string exchange, std::string queue, uint32 exchange_count)
  : m_channel(new AmqpChannel(conn)), m_exchange(exchange), m_queue(queue), m_exchange_count(exchange_count), m_next_exchange(0) {

  for (uint32 i = 0; i < m_exchange_count; ++i)
    m_channel->exchange_declare(exchange + "_" + vw::stringify(i), "direct", false, false);
  m_channel->queue_declare(queue, false, true, true);

  this->Reset();
}

AmqpRpcEndpoint::~AmqpRpcEndpoint() {
  unbind_service();
}

void AmqpRpcEndpoint::serialize_message(const ::google::protobuf::Message& message, ByteArray& bytes) {
  bytes.resize(message.ByteSize(), false);
  message.SerializeToArray(bytes.begin(), message.ByteSize());
}

void AmqpRpcEndpoint::send_message(const ::google::protobuf::Message& message, std::string routing_key) {
  ByteArray raw;
  serialize_message(message, raw);
  send_bytes(raw, routing_key);
}

void AmqpRpcEndpoint::send_bytes(ByteArray const& message, std::string routing_key) {
  m_channel->basic_publish(message, m_exchange + "_" + vw::stringify(m_next_exchange++), routing_key);
  m_next_exchange %= m_exchange_count;
}

void AmqpRpcEndpoint::get_bytes(SharedByteArray& bytes, vw::int32 timeout) {
  if (timeout == -1)
    this->m_incoming_messages.wait_pop(bytes);
  else {
    if (!this->m_incoming_messages.timed_wait_pop(bytes, timeout)) {
      vw_throw(AMQPTimeout() << "Timeout");
    }
  }
}

void AmqpRpcEndpoint::bind_service(boost::shared_ptr<google::protobuf::Service> service,
                  std::string routing_key) {
  if (m_consumer)
    vw_throw(vw::ArgumentErr() << "AmqpRpcEndpoint::bind_service(): unbind your service before you start another one");

  m_service = service;
  m_routing_key = routing_key;
  for (uint32 i = 0; i < m_exchange_count; ++i)
    m_channel->queue_bind(m_queue, m_exchange + "_" + vw::stringify(i), routing_key);
  m_consumer = m_channel->basic_consume(m_queue, boost::bind(&vw::ThreadQueue<SharedByteArray>::push, boost::ref(m_incoming_messages), _1));
}

void AmqpRpcEndpoint::unbind_service() {
  if (m_consumer) {
    m_consumer.reset();
    for (uint32 i = 0; i < m_exchange_count; ++i)
      m_channel->queue_unbind(m_queue, m_exchange + "_" + vw::stringify(i), m_routing_key);
    m_routing_key = "";
    m_service.reset();
  }
}
