// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/Amqp.h>

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#include <stdint.h>
#include <amqp.h>
#include <amqp_framing.h>

#include <unistd.h>
#include <assert.h>

#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>

#include <sstream>

#include <boost/shared_array.hpp>

using namespace vw;
using namespace vw::platefile;

// --------------------- Error Handling --------------------------
 
void die_on_error(int x, char const *context) {
  if (x < 0) {
    std::ostringstream msg;
    msg << "AMQP Error: " << context << " -- " << strerror(-x);
    vw::vw_throw(vw::IOErr() << msg.str());
  }
}

void die_on_amqp_error(amqp_rpc_reply_t x, char const *context) {
  std::ostringstream msg;

  switch (x.reply_type) {
  case AMQP_RESPONSE_NORMAL:
    return;
      
  case AMQP_RESPONSE_NONE:
    vw_throw(IOErr() << "AMQP Error: " << context << " -- missing RPC reply type.");

  case AMQP_RESPONSE_LIBRARY_EXCEPTION:
    vw_throw(IOErr() << "AMQP Error: " << context << " -- "
                     << (x.library_errno ? strerror(x.library_errno) : "(end-of-stream)"));

  case AMQP_RESPONSE_SERVER_EXCEPTION:
    switch (x.reply.id) {
    case AMQP_CONNECTION_CLOSE_METHOD: {
      amqp_connection_close_t *m = (amqp_connection_close_t *) x.reply.decoded;
      vw_throw(IOErr() << "AMQP Error: " << context << " -- "
                       << "server connection error " << m->reply_code << ", message: "
                       << (char *) m->reply_text.bytes);
    }
    case AMQP_CHANNEL_CLOSE_METHOD: {
      amqp_channel_close_t *m = (amqp_channel_close_t *) x.reply.decoded;
      vw_throw(IOErr()  << "AMQP Error: " << context << " -- "
                        << "server channel error " << m->reply_code << ", message: "
                        << (char *) m->reply_text.bytes);
    }
    default:
      vw_throw(IOErr() << "AMQP Error: " << context << " -- "
                       << "unknown server error, method id " << x.reply.id);
    }
    break;
  }
  vw_throw(IOErr() << "AMQP Error: unknown response type.");
}

// ------------------------------------------------------------------------------
//                       AmqpConnection State Structure
// ------------------------------------------------------------------------------

struct vw::platefile::AmqpConnectionState {
  int sockfd;
  amqp_connection_state_t conn;

  // -- Workaround --
  // 
  // These are standins fro AMQP_EMPTY_BYTES and AMQP_EMPTY_TABLE from
  // amqp.h, which don't seem to work properly in C++.  We sub in
  // these, which do the same thing, instead.
  amqp_table_t empty_table;
  amqp_bytes_t empty_bytes;

  AmqpConnectionState() {
    empty_table.num_entries=0;
    empty_table.entries = NULL;

    empty_bytes.len = 0;
    empty_bytes.bytes = NULL;
  }  
};

// ------------------------------------------------------------------------------
//                   AmqpConnection Implementation
// ------------------------------------------------------------------------------


// ------------------------
// Constructor / destructor
// ------------------------
 
AmqpConnection::AmqpConnection(std::string const& hostname, int port) {
  Mutex::Lock lock(m_mutex);
  
  m_state = boost::shared_ptr<AmqpConnectionState>(new AmqpConnectionState());

  // Open a socket and establish an amqp connection
  die_on_error(m_state->sockfd = amqp_open_socket(hostname.c_str(), port), "Opening socket");
  m_state->conn = amqp_new_connection();
  if (! m_state->conn)
    vw_throw(IOErr() << "Failed to open an amqp connection");
  amqp_set_sockfd(m_state->conn, m_state->sockfd);

  // Login to the Amqp Server
  die_on_amqp_error(amqp_login(m_state->conn, "/", 0, 131072, 0, AMQP_SASL_METHOD_PLAIN, 
                               "guest", "guest"), "Logging in");
  amqp_channel_open(m_state->conn, 1);
  die_on_amqp_error(amqp_rpc_reply, "Opening channel");
}

AmqpConnection::~AmqpConnection() {
  Mutex::Lock lock(m_mutex);

  vw_out(InfoMessage, "platefile::amqp") << "Closing AMQP connection.\n";
  die_on_amqp_error(amqp_channel_close(m_state->conn, 1, AMQP_REPLY_SUCCESS), "Closing channel");
  die_on_amqp_error(amqp_connection_close(m_state->conn, AMQP_REPLY_SUCCESS), "Closing connection");
  amqp_destroy_connection(m_state->conn);
  die_on_error(close(m_state->sockfd), "Closing socket");
}


// ------------------------------------------------------
// Methods for modifying exchanges, queues, and bindings.  
// ------------------------------------------------------

void AmqpConnection::queue_declare(std::string const& queue_name, bool durable, 
                                   bool exclusive, bool auto_delete) {
  Mutex::Lock lock(m_mutex);

  amqp_queue_declare_ok_t *r = amqp_queue_declare(m_state->conn, 1, 
                                                  amqp_cstring_bytes(queue_name.c_str()),
                                                  0, 0, 0, 1,
                                                  m_state->empty_table);
  die_on_amqp_error(amqp_rpc_reply, "Declaring queue");
}

void AmqpConnection::exchange_declare(std::string const& exchange_name, 
                                      std::string const& exchange_type, 
                                      bool durable, bool auto_delete) {
  Mutex::Lock lock(m_mutex);

  amqp_exchange_declare(m_state->conn, 1, 
                        amqp_cstring_bytes(exchange_name.c_str()), 
                        amqp_cstring_bytes(exchange_type.c_str()),
			0, durable, auto_delete, m_state->empty_table);
  die_on_amqp_error(amqp_rpc_reply, "Declaring exchange");
}

void AmqpConnection::queue_bind(std::string const& queue, std::string const& exchange, 
                                std::string const& routing_key) {
  Mutex::Lock lock(m_mutex);

  amqp_queue_bind(m_state->conn, 1, 
                  amqp_cstring_bytes(queue.c_str()), 
                  amqp_cstring_bytes(exchange.c_str()), 
                  amqp_cstring_bytes(routing_key.c_str()),
		  m_state->empty_table);
  die_on_amqp_error(amqp_rpc_reply, "Binding queue");
}

void AmqpConnection::queue_unbind(std::string const& queue, std::string const& exchange, 
                                  std::string const& routing_key) {
  Mutex::Lock lock(m_mutex);

  amqp_queue_unbind(m_state->conn, 1, 
                    amqp_cstring_bytes(queue.c_str()), 
                    amqp_cstring_bytes(exchange.c_str()), 
                    amqp_cstring_bytes(routing_key.c_str()),
                    m_state->empty_table);
  die_on_amqp_error(amqp_rpc_reply, "Unbinding queue");
}


// ------------------------------------------------------
// Methods for sending and receiving data
// ------------------------------------------------------

void AmqpConnection::basic_publish(boost::shared_array<uint8> const& message, 
                                   int32 size,
                                   std::string const& exchange, 
                                   std::string const& routing_key) {
  Mutex::Lock lock(m_mutex);

  // The delivery mode flag (below) can be used to select for
  // persistent delivery mode, which virtually gurantees that the
  // message will get through, even if the AMQP server crashes.
  // However, this comes at a performance penalty, so it is disabled
  // here for now.
  //
  // amqp_basic_properties_t props;
  // props._flags = AMQP_BASIC_CONTENT_TYPE_FLAG | AMQP_BASIC_DELIVERY_MODE_FLAG;
  // props.content_type = amqp_cstring_bytes("text/plain");
  // props.delivery_mode = 2; // persistent delivery mode
  amqp_bytes_t raw_data;
  raw_data.len = size;
  raw_data.bytes = (void*)(message.get());

  die_on_error(amqp_basic_publish(m_state->conn,
                                  1,
                                  amqp_cstring_bytes(exchange.c_str()),
                                  amqp_cstring_bytes(routing_key.c_str()),
                                  0,
                                  0,
                                  NULL,
                                  raw_data ),
               "Publishing");
}

bool operator<(const struct timeval& t1, const struct timeval& t2) {
  if (t1.tv_sec < t2.tv_sec)
    return true;
  if (t1.tv_sec > t2.tv_sec)
    return false;
  return (t1.tv_sec < t2.tv_sec);
}

static boost::shared_array<uint8> pump_queue(amqp_connection_state_t& conn, 
                                             amqp_frame_t &header) {

  amqp_maybe_release_buffers(conn);
  int result = amqp_simple_wait_frame(conn, &header);

  VW_ASSERT(result > 0,                             AMQPAssertion() << "AMQP error: unknown result waiting for frame.");
  VW_ASSERT(header.frame_type == AMQP_FRAME_HEADER, AMQPAssertion() << "Expected AMQP header!");

  size_t body_size = header.payload.properties.body_size;
  boost::shared_array<uint8> payload(new uint8[body_size]);

  size_t body_read = 0;
  amqp_frame_t frame;

  while (body_read < body_size) {
    result = amqp_simple_wait_frame(conn, &frame);

    VW_ASSERT(result > 0,
        AMQPAssertion() << "AMQP error: unknown result waiting for frame.");
    VW_ASSERT(frame.frame_type == AMQP_FRAME_BODY,
        AMQPAssertion() << "Expected AMQP message body!");
    VW_ASSERT(body_read + frame.payload.body_fragment.len <= body_size,
        AMQPAssertion() << "AMQP packet body size does not match header\'s body target.");

    // Copy the bytes out of the payload...
    memcpy(((char*)(payload.get()) + body_read),
        (const char*)(frame.payload.body_fragment.bytes),
        frame.payload.body_fragment.len);

    // ... and update the number of bytes we have received
    body_read += frame.payload.body_fragment.len;
  }

  return payload;
}

// timeout in ms
boost::shared_array<uint8> AmqpConnection::basic_get(std::string const& queue, 
                                                     bool no_ack,
                                                     vw::int64 timeout, 
                                                     int retries) {

  Mutex::Lock lock(m_mutex);

  struct timeval current, stop;
  int num_tries = 0;
  amqp_rpc_reply_t reply;

  const amqp_bytes_t queue_b = amqp_cstring_bytes(queue.c_str());

  if (timeout == -1) {
    stop.tv_sec = std::numeric_limits<time_t>::max();
  } else {
    gettimeofday(&stop, NULL);

    time_t sec = timeout / 1000;
    stop.tv_sec  +=  sec;
    stop.tv_usec += (timeout - (sec * 1000)) * 1000;
  }

  while (1) {
    gettimeofday(&current, NULL);

    if (! (current < stop) ) { // too lazy to define all cmp

      // Reset the timer and increment the number of tries
      gettimeofday(&stop, NULL);
      time_t sec = timeout / 1000;
      stop.tv_sec  +=  sec;
      stop.tv_usec += (timeout - (sec * 1000)) * 1000;
      ++num_tries;

      // Have we run out of retries?  If so, bail!
      if (num_tries == retries)
        vw_throw(AMQPTimeout() << "basic_get timed out");
    }

    amqp_maybe_release_buffers(m_state->conn);
    reply = amqp_basic_get(m_state->conn, 1, queue_b, no_ack);
    die_on_amqp_error(reply, "Getting");

    if (reply.reply.id == AMQP_BASIC_GET_EMPTY_METHOD) {
      usleep(1000); // yield for a bit
      continue;
    }
    if (reply.reply.id == AMQP_BASIC_GET_OK_METHOD)
      break;

    vw_throw(AMQPAssertion() << "Illegal AMQP response");
  }

  //amqp_basic_get_ok_t* get_ok = static_cast<amqp_basic_get_ok_t*>(reply.reply.decoded);

  amqp_frame_t header;
  return pump_queue(m_state->conn, header);
}

#if 1
// XXX: Fix this before uncommenting: consume starts a delivery process on the
// server every time it's called. It's designed to be used with a callback and
// a message queue pumper.
boost::shared_array<uint8> AmqpConnection::basic_consume(std::string const& queue, 
                                                         std::string &routing_key,
                                                         bool no_ack) {
  Mutex::Lock lock(m_mutex);

  amqp_basic_consume(m_state->conn, 
                     1,  // channel
                     amqp_cstring_bytes(queue.c_str()), 
                     m_state->empty_bytes, 
                     0,  // no_local
                     no_ack, 
                     0);  // exclusive

  // This last argument doesn't seem to exist in my version of rabbitmq-c.  Maybe it's new?
  // TODO: This is necessary to support AMQP 0.9.1, but we don't use that yet
  //                     m_state->empty_table); 
  die_on_amqp_error(amqp_rpc_reply, "Consuming");

  {
    amqp_frame_t frame;
    int result;

    amqp_basic_deliver_t *d;
    amqp_basic_properties_t *p;
    size_t body_target;
    size_t body_received;

    boost::shared_array<uint8> payload;

    amqp_maybe_release_buffers(m_state->conn);
    result = amqp_simple_wait_frame(m_state->conn, &frame);
    //    printf("Result %d\n", result);
    if (result <= 0)
      vw_throw(IOErr() << "AMQP error: unknown result code.");

    //    printf("Frame type %d, channel %d\n", frame.frame_type, frame.channel);
    if (frame.frame_type != AMQP_FRAME_METHOD)
      vw_throw(IOErr() << "AMQP error: unknown frame type.");

    //    printf("Method %s\n", amqp_method_name(frame.payload.method.id));
    if (frame.payload.method.id != AMQP_BASIC_DELIVER_METHOD)
      vw_throw(IOErr() << "AMQP error: unknown payload method.");

    d = (amqp_basic_deliver_t *) frame.payload.method.decoded;
    //    printf("Delivery %u, exchange %.*s routingkey %.*s\n",
           // (unsigned) d->delivery_tag,
           // (int) d->exchange.len, (char *) d->exchange.bytes,
           // (int) d->routing_key.len, (char *) d->routing_key.bytes);
    
    // copy the routing key to pass out to the caller
    boost::shared_array<char> rk_bytes(new char[d->routing_key.len + 1]);
    strncpy(rk_bytes.get(), (char*)(d->routing_key.bytes), d->routing_key.len);
    rk_bytes[d->routing_key.len] = '\0';
    routing_key = rk_bytes.get();

    result = amqp_simple_wait_frame(m_state->conn, &frame);
    if (result <= 0)
      vw_throw(IOErr() << "AMQP error: unknown result waiting for frame.");

    if (frame.frame_type != AMQP_FRAME_HEADER) {
      vw_out(ErrorMessage, "platefile::amqp") << "Expected AMQP header!";
      abort();
    }
    p = (amqp_basic_properties_t *) frame.payload.properties.decoded;
    if (p->_flags & AMQP_BASIC_CONTENT_TYPE_FLAG) {
      // printf("Content-type: %.*s\n",
      //        (int) p->content_type.len, (char *) p->content_type.bytes);
    }
    // printf("----\n");

    body_target = frame.payload.properties.body_size;
    body_received = 0;
    payload = boost::shared_array<uint8>( new uint8[body_target] );
    
    while (body_received < body_target) {
      result = amqp_simple_wait_frame(m_state->conn, &frame);
      if (result <= 0)
        vw_throw(IOErr() << "AMQP error: unknown result waiting for frame.");
      
      if (frame.frame_type != AMQP_FRAME_BODY) {
        vw_out(ErrorMessage, "platefile::amqp") << "Expected AMQP message body!";
        abort();
      }	  

      // Copy the bytes out of the payload...
      memcpy(((char*)(payload.get()) + body_received), 
             (const char*)(frame.payload.body_fragment.bytes),
             frame.payload.body_fragment.len);

      // ... and update the number of bytes we have received
      body_received += frame.payload.body_fragment.len;
      VW_ASSERT(body_received <= body_target, 
                IOErr() << "AMQP packet body size exceeded header\'s body target.");
    }
    
    // Check to make sure that we have received the right number of bytes.
    if (body_received != body_target) {
      vw_throw(IOErr() << "AMQP packet body size does not match header\'s body target.");
    }

    if (!no_ack)
      amqp_basic_ack(m_state->conn, 1, d->delivery_tag, 0);

    vw_out(VerboseDebugMessage, "platefile::amqp") << "Message payload: \"" 
                                                   << payload.get() << "\"\n";

    return payload;
  }

}
#endif
