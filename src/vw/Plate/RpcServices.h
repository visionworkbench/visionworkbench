// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_RPC_SERVICES_H__
#define __VW_PLATE_RPC_SERVICES_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Thread.h>
#include <vw/Plate/ProtoBuffers.pb.h>
#include <vw/Plate/Amqp.h>
#include <vw/Plate/Exception.h>

#include <google/protobuf/service.h>

#include <boost/shared_array.hpp>

namespace vw {
namespace platefile {

  // -----------------------------------------------------------------------
  //                           AmqpRpcClient
  // -----------------------------------------------------------------------

  class AmqpRpcClient : public google::protobuf::RpcController {
    bool m_failed;
    std::string m_failed_reason;
    long m_timeout_millis;

  public: 

    AmqpRpcClient() : m_failed(false), m_failed_reason(""), m_timeout_millis(5000) {}
    virtual ~AmqpRpcClient() {}

    // Resets the RpcController to its initial state so that it may be reused in
    // a new call.  Must not be called while an RPC is in progress.
    virtual void Reset() {
      m_failed = false;
      m_failed_reason = "";
      m_timeout_millis = 5000; // 5 seconds
    }

    // After a call has finished, returns true if the call failed.  The possible
    // reasons for failure depend on the RPC implementation.  Failed() must not
    // be called before a call has finished.  If Failed() returns true, the
    // contents of the response message are undefined.
    virtual bool Failed() const {
      vw_throw(NoImplErr() << "AmqpRpcClient::Failed() has not been implemented.  "
               << "Use normal exception throwing/catching to handle RPC errors.");
      return m_failed;
    }

    // If Failed() is true, returns a human-readable description of the error.
    virtual std::string ErrorText() const {
      vw_throw(NoImplErr() << "AmqpRpcClient::ErrorText() has not been implemented.  "
               << "Use normal exception throwing/catching to handle RPC errors.");
      return m_failed_reason;
    }

    // Advises the RPC system that the caller desires that the RPC call be
    // canceled.  The RPC system may cancel it immediately, may wait awhile and
    // then cancel it, or may not even cancel the call at all.  If the call is
    // canceled, the "done" callback will still be called and the RpcController
    // will indicate that the call failed at that time.
    virtual void StartCancel() {}


    void set_timeout(long milliseconds) { 
      m_timeout_millis = milliseconds;
    }
  
    long timeout() const {
      return m_timeout_millis;
    }

    // Returns a good private queue name. Identifier should be something unique
    // to the service, like "remote_index" or "mod_plate"
    static std::string UniqueQueueName(const std::string identifier);


    // Server-side methods ---------------------------------------------

    // These calls may be made from the server side only.  Their results
    // are undefined on the client side (may crash).

    // Causes Failed() to return true on the client side.  "reason" will be
    // incorporated into the message returned by ErrorText().  If you find
    // you need to return machine-readable information about failures, you
    // should incorporate it into your response protocol buffer and should
    // NOT call SetFailed().
    virtual void SetFailed(const std::string& reason) {
      vw_throw(NoImplErr() << "AmqpRpcClient::SetFailed() has not been implemented.  "
               << "Use normal exception throwing/catching to handle RPC errors.");
      m_failed = true;
      m_failed_reason = reason;
    }

    // If true, indicates that the client canceled the RPC, so the server may
    // as well give up on replying to it.  The server should still call the
    // final "done" callback.
    virtual bool IsCanceled() const { return false; }

    // Asks that the given callback be called when the RPC is canceled.  The
    // callback will always be called exactly once.  If the RPC completes without
    // being canceled, the callback will be called after completion.  If the RPC
    // has already been canceled when NotifyOnCancel() is called, the callback
    // will be called immediately.
    //
    // NotifyOnCancel() must be called no more than once per request.
    virtual void NotifyOnCancel(google::protobuf::Closure* callback) {}

  };


  // -----------------------------------------------------------------------------
  //                                AmqpRpcServer
  // -----------------------------------------------------------------------------

  class NetworkMonitor {
    size_t m_total_bytes;
    int32 m_total_queries;
    vw::Mutex m_mutex;

  public:
    
    NetworkMonitor() : 
      m_total_bytes(0), m_total_queries(0) {}
    
    virtual ~NetworkMonitor() {}

    void record_query(size_t num_bytes) {
      Mutex::Lock lock(m_mutex);
      ++m_total_queries;
      m_total_bytes += num_bytes;
    }

    size_t bytes_processed() { 
      Mutex::Lock lock(m_mutex);
      return m_total_bytes; 
    }

    int queries_processed() { 
      Mutex::Lock lock(m_mutex);
      return m_total_queries; 
    };

    void reset() {
      Mutex::Lock lock(m_mutex);
      m_total_bytes = 0;
      m_total_queries = 0;
    }
  };

  class AmqpRpcServer : public google::protobuf::RpcController {
    std::string m_exchange;
    std::string m_queue;
    AmqpConnection m_conn;
    bool m_failed;
    std::string m_failed_reason;
    NetworkMonitor m_stats;
    bool m_terminate;
    bool m_debug;

    boost::shared_ptr<google::protobuf::Service> m_service;
    
  public:
    AmqpRpcServer(std::string exchange, std::string queue, 
                  std::string hostname = "localhost", int port = 5672, 
                  bool debug = false) :
      m_exchange(exchange), m_queue(queue), m_conn(hostname, port), 
      m_terminate(false), m_debug(debug) {
      
      m_conn.exchange_declare(exchange, "direct", false, false);
      m_conn.queue_declare(queue, false, true, false);
      
      this->Reset();
    }

    virtual ~AmqpRpcServer() {}

    void export_with_routing_key(boost::shared_ptr<google::protobuf::Service> service,
                                 std::string routing_key) {
      m_conn.queue_bind(m_queue, m_exchange, routing_key);
      std::cout << "\t    Bound \'" << m_queue << "\' to \'" << m_exchange 
                << "\' with the \'" << routing_key << "\' routing key.\n";
      m_service = service;
    }

    /// Run the server run loop
    void run();

    /// Shut down the server
    void shutdown() { m_terminate = true; }

    /// Return statistics about number of messages processed per second.
    int queries_processed() {
      return m_stats.queries_processed();
    }

    /// Return statistics about number of messages processed per second.
    size_t bytes_processed() {
      return m_stats.bytes_processed();
    }

    void reset_stats() {
      m_stats.reset();
    }

    // ------------------ RpcController Methods ---------------------

    // Resets the RpcController to its initial state so that it may be reused in
    // a new call.  Must not be called while an RPC is in progress.
    virtual void Reset() {
      m_failed = false;
      m_failed_reason = "";
    }

    // After a call has finished, returns true if the call failed.  The possible
    // reasons for failure depend on the RPC implementation.  Failed() must not
    // be called before a call has finished.  If Failed() returns true, the
    // contents of the response message are undefined.
    virtual bool Failed() const {
      vw_throw(NoImplErr() << "AmqpRpcServer::Failed() has not been implemented.  "
               << "Use normal exception throwing/catching to handle RPC errors.");
      return m_failed;
    }
  
    // If Failed() is true, returns a human-readable description of the error.
    virtual std::string ErrorText() const {
      vw_throw(NoImplErr() << "AmqpRpcServer::ErrorText() has not been implemented.  "
               << "Use normal exception throwing/catching to handle RPC errors.");
      return m_failed_reason;
    }

    // Advises the RPC system that the caller desires that the RPC call be
    // canceled.  The RPC system may cancel it immediately, may wait awhile and
    // then cancel it, or may not even cancel the call at all.  If the call is
    // canceled, the "done" callback will still be called and the RpcController
    // will indicate that the call failed at that time.
    virtual void StartCancel() {}  

    // Server-side methods ---------------------------------------------
  
    // Causes Failed() to return true on the client side.  "reason" will be
    // incorporated into the message returned by ErrorText().  If you find
    // you need to return machine-readable information about failures, you
    // should incorporate it into your response protocol buffer and should
    // NOT call SetFailed().
    virtual void SetFailed(const std::string& reason) {
      vw_throw(NoImplErr() << "AmqpRpcServer::SetFailed() has not been implemented.  "
               << "Use normal exception throwing/catching to handle RPC errors.");
      m_failed = true;
      m_failed_reason = reason;
    }

    void SetFailed(std::string type, std::string description) {

    }



    // If true, indicates that the client canceled the RPC, so the server may
    // as well give up on replying to it.  The server should still call the
    // final "done" callback.
    virtual bool IsCanceled() const { return false; }
  
    // Asks that the given callback be called when the RPC is canceled.  The
    // callback will always be called exactly once.  If the RPC completes without
    // being canceled, the callback will be called after completion.  If the RPC
    // has already been canceled when NotifyOnCancel() is called, the callback
    // will be called immediately.
    //
    // NotifyOnCancel() must be called no more than once per request.
    virtual void NotifyOnCancel(google::protobuf::Closure* callback) {}

  };


  // -----------------------------------------------------------------------
  //                           AmqpRpcChannel
  // -----------------------------------------------------------------------

  class AmqpRpcChannel : public google::protobuf::RpcChannel {

    std::string m_exchange, m_request_routing_key, m_response_queue;
    AmqpConnection m_conn;
    vw::Mutex m_mutex;

  public:
  
    AmqpRpcChannel(std::string const& exchange, 
                   std::string const& request_routing_key, 
                   std::string const& response_queue, 
                   std::string hostname = "localhost", 
                   int port = 5672 );

    ~AmqpRpcChannel() {}

    // Client side RPC implementation using AMQP as the message
    // passing interface.
    virtual void CallMethod(const google::protobuf::MethodDescriptor* method,
                            google::protobuf::RpcController* controller,
                            const google::protobuf::Message* request,
                            google::protobuf::Message* response,
                            google::protobuf::Closure* done);
  };

  /// A handy utility class for serializing/deserializing protocol
  /// buffers over the wire.  Store the buffer as a stream of bytes with
  /// the size of the message at the beginning of the stream.
  class WireMessage {

    typedef int32 size_type;

    boost::shared_array<uint8> m_serialized_bytes;
    size_type m_payload_size;

    uint8* message_start() const {
      return (m_serialized_bytes.get() + sizeof(size_type));
    }

    public:

    WireMessage(const google::protobuf::Message* message) {
      m_serialized_bytes.reset( new uint8[sizeof(size_type) + message->ByteSize()] );

      // Store the size at the beginning of the byte stream.
      m_payload_size = message->ByteSize();
      ((size_type*)(m_serialized_bytes.get()))[0] = m_payload_size;

      message->SerializeToArray((void*)(message_start()), message->ByteSize());
    }

    WireMessage(boost::shared_array<uint8> const& serialized_bytes) {
      m_serialized_bytes = serialized_bytes;
      m_payload_size = *( (size_type*)(m_serialized_bytes.get()) );
    }

    boost::shared_array<uint8> serialized_bytes() const { return m_serialized_bytes; }
    size_type size() const { return m_payload_size + sizeof(size_type); }

    template <class MessageT>
      MessageT parse_as_message() {
        MessageT message;
        bool success = message.ParseFromArray((void*)(message_start()), m_payload_size);
        if (!success)
          vw_throw(vw::platefile::RpcErr() << "Could not deserialize message.");
        return message;
      }

    void parse(google::protobuf::Message* message) {
      message->ParseFromArray((void*)(message_start()), m_payload_size);
    }

  };

}} // namespace vw::platefile

#endif // __VW_PLATE_RPC_SERVICES_H__
