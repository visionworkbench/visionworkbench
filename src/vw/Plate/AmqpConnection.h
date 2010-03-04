// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_AMQP__
#define __VW_PLATE_AMQP__

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/utility.hpp>
#include <boost/function.hpp>

#include <vw/Core/Log.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>
#include <vw/Core/VarArray.h>
#include <vw/Plate/Exception.h>

#include <set>

#include <amqp.h>

namespace vw {
namespace platefile {

  typedef vw::VarArray<uint8> ByteArray;
  typedef boost::shared_ptr<ByteArray> SharedByteArray;

  VW_DEFINE_EXCEPTION(AMQPErr,     PlatefileErr);
  VW_DEFINE_EXCEPTION(AMQPTimeout, AMQPErr);
  VW_DEFINE_EXCEPTION(AMQPEof,     AMQPErr);

  // This exception denotes a potentially desynchronizing AMQP error. Only
  // recovery mechanism is to recreate the connection.
  VW_DEFINE_EXCEPTION(AMQPAssertion, AMQPErr);

  // Forward declaration because amqp.h is gross
  class AmqpChannel;
  class AmqpConsumer;
  class AmqpConsumeTask;

  typedef boost::remove_pointer<amqp_connection_state_t>::type AmqpConnectionState;

  class AmqpConnection : private boost::noncopyable {
    private:
      // Any access (read or write) to m_state needs to hold the m_state_mutex lock
      boost::shared_ptr<AmqpConnectionState> m_state;
      vw::Mutex m_state_mutex;
      std::set<int16> m_used_channels;

      /// Allocates a communications channel number. Meant to be called by
      ///     AmqpChannel's constructor. Grab the m_state_mutex BEFORE you call this.
      /// @param channel Request a specific channel number. -1 for don't care.
      //                 Requesting a particular one will throw if it's already used.
      int16 get_channel(int16 channel = -1);

      friend class AmqpChannel;
      friend class AmqpConsumeTask;

    public:
      /// Open a new connection to the AMQP server.  This connection
      /// terminates automatically when this object is destroyed.
      AmqpConnection(std::string const& hostname = "localhost", int port = 5672);

      /// Closes the AMQP connection and destroys this object.
      ~AmqpConnection();

  };

  class AmqpChannel : private boost::noncopyable {
    private:
      boost::shared_ptr<AmqpConnection> m_conn;
      int16 m_channel;

    public:
      AmqpChannel(boost::shared_ptr<AmqpConnection> conn, int16 channel = -1);
      virtual ~AmqpChannel();

      void exchange_declare(std::string const& exchange_name, std::string const& exchange_type, bool durable, bool auto_delete);
      void queue_declare(std::string const& queue_name, bool durable, bool exclusive, bool auto_delete);
      void queue_bind(std::string const& queue, std::string const& exchange, std::string const& routing_key);
      void queue_unbind(std::string const& queue, std::string const& exchange, std::string const& routing_key);

      void basic_publish(ByteArray const& message, std::string const& exchange, std::string const& routing_key);

      boost::shared_ptr<AmqpConsumer> basic_consume(std::string const& queue, boost::function<void (SharedByteArray)> callback);
      bool basic_get(std::string const& queue, SharedByteArray& message);
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
