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

  class AmqpRpcChannel : public vw::platefile::AmqpChannel,  public google::protobuf::RpcChannel {
    public:
      AmqpRpcChannel(boost::shared_ptr<AmqpConnection> conn, int16 channel = -1):
        AmqpChannel(conn, channel) {};

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
    private:
      typedef int32 size_type;
      size_type m_payload_size;

      // XXX: VarArray only supports deep-copy. D'oh. Having a shared_array
      // inside a shared_ptr is lame. Revisit this.
      typedef vw::VarArray<uint8> ContainerT;
      NativeMessage data;

      uint8* message_start() const {
        return (data->begin() + sizeof(size_type));
      }
    public:
      WireMessage(const google::protobuf::Message* message)
        : m_payload_size(message->ByteSize()), data(new ContainerT(sizeof(size_type) + m_payload_size)) {

          reinterpret_cast<size_type*>(data->begin())[0] = htonl(m_payload_size);
          message->SerializeToArray((void*)(message_start()), message->ByteSize());
      }

      WireMessage(NativeMessage const& serialized_bytes) {
        data = serialized_bytes;
        m_payload_size = ntohl(reinterpret_cast<size_type*>(data->begin())[0]);
      }

      NativeMessage serialized_bytes() const { return data; }
      size_type size() const { return data->size(); }

      template <class MessageT>
        MessageT parse_as_message() {
          MessageT message;
          bool success = message.ParseFromArray(reinterpret_cast<void*>(message_start()), m_payload_size);
          if (!success)
            vw_throw(vw::platefile::RpcErr() << "Could not deserialize message.");
          return message;
        }

      void parse(google::protobuf::Message* message) {
        message->ParseFromArray(reinterpret_cast<void*>(message_start()), m_payload_size);
      }

  };

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

      boost::shared_ptr<AmqpRpcChannel> m_channel;
      boost::shared_ptr<google::protobuf::Service> m_service;
      std::string m_exchange;
      std::string m_queue;
      std::string m_routing_key;
      boost::shared_ptr<AmqpConsumer> m_consumer;

      vw::ThreadQueue<NativeMessage> messages;

    public:
      AmqpRpcEndpoint(boost::shared_ptr<AmqpConnection> conn, std::string exchange, std::string queue)
        : m_channel(new AmqpRpcChannel(conn)), m_exchange(exchange), m_queue(queue) {

        m_channel->exchange_declare(exchange, "direct", false, false);
        m_channel->queue_declare(queue, false, true, true);

        this->Reset();
      }

      virtual ~AmqpRpcEndpoint() {
        unbind_service();
      }

      void send_message(const ::google::protobuf::Message& message) {
        WireMessage msg(&message);
        messages.flush();
        m_channel->basic_publish( msg.serialized_bytes(), m_exchange, m_routing_key);
      }

      template <typename MessageT>
      void get_message(MessageT& message, vw::int32 timeout = -1) {
        NativeMessage response_bytes;

        if (timeout == -1)
          this->messages.wait_pop(response_bytes);
        else {
          if (!this->messages.timed_wait_pop(response_bytes, timeout)) {
            vw_throw(AMQPTimeout() << "Timeout");
          }
        }

        message = WireMessage(response_bytes).parse_as_message<MessageT>();
      }


      void bind_service(boost::shared_ptr<google::protobuf::Service> service,
          std::string routing_key) {
        if (m_consumer)
          vw_throw(vw::ArgumentErr() << "AmqpRpcEndpoint::bind_service(): unbind your service before you start another one");

        m_service = service;
        m_routing_key = routing_key;
        m_channel->queue_bind(m_queue, m_exchange, routing_key);
        m_consumer = m_channel->basic_consume(m_queue, boost::bind(&vw::ThreadQueue<NativeMessage>::push, boost::ref(messages), _1));
      }

      void unbind_service() {
        if (m_consumer) {
          m_consumer.reset();
          m_channel->queue_unbind(m_queue, m_exchange, m_routing_key);
          m_routing_key = "";
          m_service.reset();
        }
      }

      boost::shared_ptr<AmqpRpcChannel> channel() {
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

  class AmqpRpcClient : public AmqpRpcEndpoint {
    public:
      AmqpRpcClient(boost::shared_ptr<AmqpConnection> conn, std::string exchange, std::string queue) :
        AmqpRpcEndpoint(conn, exchange, queue) {}

      virtual ~AmqpRpcClient() {}

      // Returns a good private queue name. Identifier should be something unique
      // to the service, like "remote_index" or "mod_plate"
      static std::string UniqueQueueName(const std::string identifier);
  };


  class AmqpRpcServer : public AmqpRpcEndpoint {
      NetworkMonitor m_stats;
      bool m_terminate;
      bool m_debug;
    public:
      AmqpRpcServer(boost::shared_ptr<AmqpConnection> conn, std::string exchange, std::string queue, bool debug = false) :
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
