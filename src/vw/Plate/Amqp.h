// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_AMQP__
#define __VW_PLATE_AMQP__

#include <string>
#include <boost/shared_ptr.hpp>

namespace vw {
namespace platefile {

  // Forward declaration
  struct AmqpConnectionState;

  class AmqpConnection {
    boost::shared_ptr<AmqpConnectionState> m_state;
    
  public: 

    // ------------------------------------------------------
    // Constructor / destructor
    // ------------------------------------------------------

    /// Open a new connection to the AMQP server.  This connection
    /// terminates automatically when this object is destroyed.
    AmqpConnection(std::string const& hostname = "localhost", int port = 5672);

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

    void basic_publish(std::string const& message, 
                       std::string const& exchange, 
                       std::string const& routing_key);

    template <class ProtoBufT>
    void basic_publish_protobuf(ProtoBufT const& protobuf, 
                                std::string const& exchange, 
                                std::string const& routing_key) {
      std::string message_bytes;
      protobuf.SerializeToString(&message_bytes);
      this->basic_publish(message_bytes, exchange, routing_key);
    }
    

    std::string basic_consume(std::string const& queue, 
                              std::string &routing_key,
                              bool no_ack);
  };

}} // namespace vw::platefile

#endif // __VW_PLATE_AMQP__
