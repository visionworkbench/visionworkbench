// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_AMQPCONNECTION_H__
#define __VW_PLATE_AMQPCONNECTION_H__

#include <vw/Plate/RpcChannel.h>
#include <vw/Plate/Exception.h>
#include <vw/Core/ThreadQueue.h>

#include <boost/function.hpp>
#include <boost/shared_array.hpp>

struct amqp_rpc_reply_t_;
typedef amqp_rpc_reply_t_ amqp_rpc_reply_t;

namespace vw {
  class Thread;

namespace platefile {

  // To recover from this, you must reconnect
  VW_DEFINE_EXCEPTION(AMQPErr, NetworkErr);

  class AmqpConsumer;
  class AmqpConsumeTask;
  class AmqpConnection;
  struct AmqpData {
    boost::shared_ptr<std::vector<uint8> > data;
    std::string sender;
  };

  class AmqpChannel : public IChannel,
                      private boost::noncopyable
  {
    private:
      boost::shared_ptr<AmqpConnection> m_conn;
      boost::shared_ptr<AmqpConsumer> m_consumer;
      vw::ThreadQueue<AmqpData> m_incoming_messages;
      int16 m_channel_id;
      std::string m_human_name;
      std::string m_local_name;
      std::string m_remote_name;
      int32 m_timeout;
      uint32 m_seq;
      uint32 m_retries;
    protected:
      void locked_check_error(const amqp_rpc_reply_t& x, const std::string& context);
      void locked_exchange_declare(std::string const& exchange_name, std::string const& exchange_type, bool durable, bool auto_delete);
      void locked_queue_declare(std::string const& queue_name, bool durable, bool exclusive, bool auto_delete);
      void locked_queue_bind(std::string const& queue, std::string const& exchange, std::string const& routing_key);
      void locked_queue_unbind(std::string const& queue, std::string const& exchange, std::string const& routing_key);
      boost::shared_ptr<AmqpConsumer> locked_basic_consume(std::string const& queue, boost::function<void (AmqpData)> callback);

      void basic_publish(const uint8* message, uint64 len, std::string const& exchange, std::string const& routing_key) const;
      //bool basic_get(std::string const& queue, SharedByteArray& message) const;

      void create_endpoint(const std::string& rabbitmq, short port, const std::string& name);
      void remove_endpoint();
      void send_bytes(const uint8* message, size_t len);
      bool recv_bytes(std::vector<uint8>* bytes);

    public:
      AmqpChannel(const std::string& human_name)
        : m_human_name(human_name), m_timeout(DEFAULT_TIMEOUT), m_seq(0), m_retries(DEFAULT_RETRIES) {}
      virtual ~AmqpChannel() VW_NOTHROW;

      int32  timeout() const;
      uint32 retries() const;

      void set_timeout(vw::int32);
      void set_retries(vw::uint32);

      //uint64 queue_depth() const;
      std::string name() const;

      virtual void CallMethod(const google::protobuf::MethodDescriptor*,
                              google::protobuf::RpcController*,
                              const google::protobuf::Message*,
                              google::protobuf::Message*,
                              google::protobuf::Closure*);

      // url format:
      // amqp://${rabbit_ip}:${port}/${exchange}
      // or
      // pf://${rabbit_ip}:${port}/${exchange}/${platefile}.plate (DEPRECATED)
      //   which maps to amqp://${rabbit_ip}:${port}/${exchange}/index
      void conn(const Url& server);
      void bind(const Url& self);
  };

  class AmqpConsumer : private boost::noncopyable {
    private:
      boost::shared_ptr<AmqpConsumeTask> m_task;
      boost::shared_ptr<Thread> m_thread;
    public:
      AmqpConsumer(boost::shared_ptr<AmqpConsumeTask> task, boost::shared_ptr<vw::Thread> thread);
      ~AmqpConsumer();
  };

}} // namespace vw::platefile

#endif // __VW_PLATE_AMQP__
