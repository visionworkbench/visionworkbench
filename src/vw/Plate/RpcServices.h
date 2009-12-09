// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_RPC_SERVICES_H__
#define __VW_PLATE_RPC_SERVICES_H__

#include <arpa/inet.h>

#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Thread.h>
#include <vw/Core/ThreadQueue.h>
#include <vw/Plate/ProtoBuffers.pb.h>
#include <vw/Plate/Amqp.h>
#include <vw/Plate/Exception.h>

#include <google/protobuf/service.h>

#include <boost/shared_array.hpp>

namespace vw {
namespace platefile {

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

  class AmqpRpcEndpoint : public google::protobuf::RpcController {

    protected:

      boost::shared_ptr<AmqpChannel> m_channel;
      boost::shared_ptr<google::protobuf::Service> m_service;
      std::string m_exchange;
      std::string m_queue;
      std::string m_routing_key;
      boost::shared_ptr<AmqpConsumer> m_consumer;

      vw::ThreadQueue<boost::shared_ptr<ByteArray> > messages;

    public:
    AmqpRpcEndpoint(boost::shared_ptr<AmqpConnection> conn, std::string exchange, std::string queue)
      : m_channel(new AmqpChannel(conn)), 
        m_exchange(exchange), m_queue(queue) {

        m_channel->exchange_declare(exchange, "direct", false, false);
        m_channel->queue_declare(queue, false, true, true);

        this->Reset();
      }

      virtual ~AmqpRpcEndpoint() {
        unbind_service();
      }

      static void serialize_message(const ::google::protobuf::Message& message, ByteArray& bytes) {
        bytes.resize(message.ByteSize(), false);
        message.SerializeToArray(bytes.begin(), message.ByteSize());
      }

      template <typename MessageT>
      static void parse_message(const ByteArray& bytes, MessageT& message) {
        if (!message.ParseFromArray(bytes.begin(), bytes.size()))
          vw_throw(vw::platefile::RpcErr() << "Could not parse bytes into a message.");
      }

    void send_message(const ::google::protobuf::Message& message, std::string routing_key) {
        ByteArray raw;
        serialize_message(message, raw);
        // XXX: this flushes out the message queue. this might not be a good
        // idea- what if the caller just hasn't processed the message yet?
        messages.flush();
        m_channel->basic_publish(raw, m_exchange, routing_key);
      }

      template <typename MessageT>
      void get_message(MessageT& message, vw::int32 timeout = -1) {
        SharedByteArray response_bytes;

        if (timeout == -1)
          this->messages.wait_pop(response_bytes);
        else {
          if (!this->messages.timed_wait_pop(response_bytes, timeout)) {
            vw_throw(AMQPTimeout() << "Timeout");
          }
        }
        parse_message(*response_bytes.get(), message);
      }


      void bind_service(boost::shared_ptr<google::protobuf::Service> service,
                        std::string routing_key) {
        if (m_consumer)
          vw_throw(vw::ArgumentErr() << "AmqpRpcEndpoint::bind_service(): unbind your service before you start another one");

        m_service = service;
        m_routing_key = routing_key;
        m_channel->queue_bind(m_queue, m_exchange, routing_key);
        m_consumer = m_channel->basic_consume(m_queue, boost::bind(&vw::ThreadQueue<SharedByteArray>::push, boost::ref(messages), _1));
      }

      void unbind_service() {
        if (m_consumer) {
          m_consumer.reset();
          m_channel->queue_unbind(m_queue, m_exchange, m_routing_key);
          m_routing_key = "";
          m_service.reset();
        }
      }

      boost::shared_ptr<AmqpChannel> channel() {
        return m_channel;
      }

      const std::string queue_name() {
        return m_queue;
      }

      virtual void Reset() { }

    protected:
      // These functions are not used. Bail out if someone tries.
      virtual bool Failed() const {
        vw_throw(NoImplErr() << "AmqpRpcEndpoint::Failed(): Use exceptions instead.");
        return false;
      }

      virtual std::string ErrorText() const {
        vw_throw(NoImplErr() << "AmqpRpcEndpoint::ErrorText(): Use exceptions instead.");
        return "";
      }

      virtual void StartCancel() {
        vw_throw(NoImplErr() << "AmqpRpcEndpoint::StartCancel(): Use exceptions instead.");
      }

      virtual void SetFailed(const std::string& reason) {
        vw_throw(NoImplErr() << "AmqpRpcEndpoint::SetFailed(): Use exceptions instead.");
      }

      virtual bool IsCanceled() const { return false; }

      virtual void NotifyOnCancel(google::protobuf::Closure* callback) {
        vw_throw(NoImplErr() << "AmqpRpcEndpoint::NotifyOnCancel(): Use exceptions instead.");
      }
  };

  class AmqpRpcClient : public AmqpRpcEndpoint,
                        public google::protobuf::RpcChannel {

    // Routing key to use when callmethod is called.  Basically this
    // is the routing key to address message to the server.
    std::string m_request_routing_key;

  public:
    AmqpRpcClient(boost::shared_ptr<AmqpConnection> conn, std::string exchange, 
                  std::string queue, std::string request_routing_key) :
      AmqpRpcEndpoint(conn, exchange, queue), m_request_routing_key(request_routing_key) {}
    
    virtual ~AmqpRpcClient() {}

    virtual void CallMethod(const google::protobuf::MethodDescriptor* method,
                            google::protobuf::RpcController* controller,
                            const google::protobuf::Message* request,
                            google::protobuf::Message* response,
                            google::protobuf::Closure* done);

      // Returns a good private queue name. Identifier should be something unique
      // to the service, like "remote_index" or "mod_plate"
      static std::string UniqueQueueName(const std::string identifier);
  };


  class AmqpRpcServer : public AmqpRpcEndpoint {
    NetworkMonitor m_stats;
    std::string m_request_routing_key;
    bool m_terminate;
    bool m_debug;

  public:

    AmqpRpcServer(boost::shared_ptr<AmqpConnection> conn, std::string exchange, 
                  std::string queue, bool debug = false) :
      AmqpRpcEndpoint(conn, exchange, queue), m_terminate(false), m_debug(debug) {}
    
    virtual ~AmqpRpcServer() {}
    
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
  };

}} // namespace vw::platefile

#endif // __VW_PLATE_RPC_SERVICES_H__
