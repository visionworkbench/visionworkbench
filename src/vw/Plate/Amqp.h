// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_AMQP__
#define __VW_PLATE_AMQP__

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include <vw/Core/Log.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>

namespace vw {
namespace platefile {

  VW_DEFINE_EXCEPTION(AMQPErr,       IOErr);
  VW_DEFINE_EXCEPTION(AMQPTimeout,   AMQPErr);

  // This exception denotes a potentially desynchronizing AMQP error. Safest
  // recovery mechanism is to recreate the connection.
  VW_DEFINE_EXCEPTION(AMQPAssertion, AMQPErr);

  // Forward declaration
  struct AmqpConnectionState;

  class AmqpConnection {
    boost::shared_ptr<AmqpConnectionState> m_state;
    vw::Mutex m_mutex;
    vw::int64 m_timeout;

  public:

    // ------------------------------------------------------
    // Constructor / destructor
    // ------------------------------------------------------

    /// Open a new connection to the AMQP server.  This connection
    /// terminates automatically when this object is destroyed. 
    //Timeout is in ms, -1 means forever.
    AmqpConnection(std::string const& hostname = "localhost", int port = 5672, vw::int64 timeout = 1000);

    /// Closes the AMQP connection and destroys this object.
    ~AmqpConnection();

    // ------------------------------------------------------
    // Methods for modifying exchanges, queues, and bindings.  
    // ------------------------------------------------------
    void queue_declare(std::string const& queue_name, bool durable, 
                       bool exclusive, bool auto_delete);

    void exchange_declare(std::string const& exchange_name, 
                          std::string const& exchange_type = "direct", 
                          bool durable = false, bool auto_delete = false);

    void queue_bind(std::string const& queue, std::string const& exchange, 
                    std::string const& routing_key);

    void queue_unbind(std::string const& queue, std::string const& exchange, 
                      std::string const& routing_key);

    // ------------------------------------------------------
    // Methods for sending and receiving data
    // ------------------------------------------------------

    void basic_publish(boost::shared_array<uint8> const& message, 
                       int32 size,
                       std::string const& exchange, 
                       std::string const& routing_key);

    boost::shared_array<uint8> basic_get(std::string const& queue,
                                          bool no_ack = false);

    //boost::shared_array<uint8> basic_consume(std::string const& queue, 
    //                                         std::string &routing_key,
    //                                         bool no_ack);
  };

}} // namespace vw::platefile

#endif // __VW_PLATE_AMQP__
