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
#include <fcntl.h>

#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>

#include <sstream>

#include <boost/shared_array.hpp>
#include <cerrno>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/ip.h>
#include <netinet/tcp.h>

using namespace vw;
using namespace vw::platefile;

namespace {

  const amqp_bytes_t amqp_string(const std::string& s);
  std::string amqp_bytes(const amqp_bytes_t& s);
  void die_on_amqp_error(amqp_rpc_reply_t x, const std::string& context);
  void die_on_error(int x, const std::string& context);
  SharedByteArray pump_queue(amqp_connection_state_t conn);
  bool select_helper(int fd, vw::int32 timeout, const std::string& context);
  void vw_simple_wait_frame(amqp_connection_state_t conn, amqp_frame_t *frame,
                            vw::int32 timeout, const std::string& context);
}

AmqpConnection::AmqpConnection(std::string const& hostname, int port) {

  m_state.reset(amqp_new_connection(), amqp_destroy_connection);
  if (!m_state.get())
    vw_throw(IOErr() << "Failed to create amqp state object");

  // Open a socket and establish an amqp connection
  int fd = amqp_open_socket(hostname.c_str(), port);
  if (fd < 0)
    vw_throw(IOErr() << "Failed to open AMQP socket.");

  int flag = 1;
  die_on_error(setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &flag, sizeof(int)), "Setting TCP_NODELAY");

  amqp_set_sockfd(m_state.get(), fd);

  // Login to the Amqp Server
  die_on_amqp_error(amqp_login(m_state.get(), "/", 0, 131072, 0, AMQP_SASL_METHOD_PLAIN,
                               "guest", "guest"), "Logging in");

  // 0 is reserved
  m_used_channels.insert(0);
}

AmqpConnection::~AmqpConnection() {
  Mutex::Lock lock(m_state_mutex);

  die_on_amqp_error(amqp_connection_close(m_state.get(), AMQP_REPLY_SUCCESS), "Closing connection");
  die_on_error(close(amqp_get_sockfd(m_state.get())), "Closing socket");
}

// CALL THIS WITH THE m_state_mutex LOCK ALREADY HELD.
int16 AmqpConnection::get_channel(int16 channel) {
  if (channel == -1) {
    if(m_used_channels.empty())
      channel = 1;
    else
      channel = *m_used_channels.rbegin() + 1; // largest elt in set
  }

  if (m_used_channels.count(channel) != 0)
    vw_throw(LogicErr() << "Channel " << channel << "is already in use");

  m_used_channels.insert(channel);
  return channel;
}


AmqpChannel::AmqpChannel(boost::shared_ptr<AmqpConnection> conn, int16 channel)
    : m_conn(conn) {

  Mutex::Lock lock(m_conn->m_state_mutex);

  m_channel = m_conn->get_channel(channel);
  amqp_channel_open(m_conn->m_state.get(), m_channel);
  die_on_amqp_error(amqp_rpc_reply, "Opening channel");

}

AmqpChannel::~AmqpChannel() {
  Mutex::Lock lock(m_conn->m_state_mutex);
  die_on_amqp_error(amqp_channel_close(m_conn->m_state.get(), m_channel, AMQP_REPLY_SUCCESS), "Closing channel");
}

void AmqpChannel::exchange_declare(std::string const& exchange_name,
                                   std::string const& exchange_type,
                                   bool durable, bool auto_delete) {
  Mutex::Lock lock(m_conn->m_state_mutex);
  amqp_exchange_declare(m_conn->m_state.get(), m_channel, amqp_string(exchange_name),
                        amqp_string(exchange_type), 0, durable, auto_delete, amqp_table_t());
  die_on_amqp_error(amqp_rpc_reply, "Declaring Exchange");
}

void AmqpChannel::queue_declare(std::string const& queue_name, bool durable,
                                bool exclusive, bool auto_delete) {
  Mutex::Lock lock(m_conn->m_state_mutex);

  amqp_queue_declare(m_conn->m_state.get(), m_channel, amqp_string(queue_name), 0,
                     durable, exclusive, auto_delete, amqp_table_t());
  die_on_amqp_error(amqp_rpc_reply, "Declaring queue");
}

void AmqpChannel::queue_bind(std::string const& queue, std::string const& exchange,
                             std::string const& routing_key) {
  Mutex::Lock(m_conn->m_state_mutex);

  amqp_queue_bind(m_conn->m_state.get(), m_channel, amqp_string(queue),
                  amqp_string(exchange), amqp_string(routing_key), amqp_table_t());
  die_on_amqp_error(amqp_rpc_reply, "Binding queue");
}

void AmqpChannel::queue_unbind(std::string const& queue, std::string const& exchange,
                               std::string const& routing_key) {
  Mutex::Lock lock(m_conn->m_state_mutex);

  amqp_queue_unbind(m_conn->m_state.get(), m_channel, amqp_string(queue),
                    amqp_string(exchange), amqp_string(routing_key), amqp_table_t());
  die_on_amqp_error(amqp_rpc_reply, "Unbinding queue");
}

void AmqpChannel::basic_publish(ByteArray const& message,
                                std::string const& exchange, std::string const& routing_key) {
  Mutex::Lock lock(m_conn->m_state_mutex);

  amqp_bytes_t raw_data;
  raw_data.len   = message.size();
  raw_data.bytes = const_cast<void*>(reinterpret_cast<const void*>(message.begin()));

  int ret = amqp_basic_publish(m_conn->m_state.get(), m_channel, amqp_string(exchange),
                               amqp_string(routing_key), 0, 0, NULL, raw_data);
  die_on_error(ret, "Publishing");
}

bool AmqpChannel::basic_get(std::string const& queue, SharedByteArray& message) {

  Mutex::Lock lock(m_conn->m_state_mutex);

  amqp_rpc_reply_t reply;
  while (true) {
    amqp_maybe_release_buffers(m_conn->m_state.get());
    reply = amqp_basic_get(m_conn->m_state.get(), m_channel, amqp_string(queue), 1);
    die_on_amqp_error(reply, "Getting");

    if (reply.reply.id == AMQP_BASIC_GET_EMPTY_METHOD) {
      usleep(100000); // yield for a bit
      continue;
    }
    if (reply.reply.id == AMQP_BASIC_GET_OK_METHOD)
      break;

    vw_throw(AMQPAssertion() << "Illegal AMQP response");
  }

  message = pump_queue(m_conn->m_state.get());
  return true;
}

namespace vw {
  namespace platefile {

class AmqpConsumeTask {
  public:
    typedef boost::shared_ptr<AmqpConnection> Connection;
    typedef boost::function<void (SharedByteArray)> Callback;

    // XXX: This should really hold a channel, and not a connection...
    Connection m_conn;
    int16 m_chan;
    Callback m_callback;
    std::string m_queue;
    std::string m_consumer_tag;
    bool go;

    AmqpConsumeTask(Connection conn, int16 channel_id, const Callback& callback,
                    std::string queue_name, std::string consumer_tag)
        : m_conn(conn), m_chan(channel_id), m_callback(callback), m_queue(queue_name),
          m_consumer_tag(consumer_tag), go(true) {
    }

    void kill() {go = false;}

    void operator()() const {
      int fd;
      {
        // This isn't totally safe, because technically the socket could be
        // changed from under us. That would break so much of the rest of the
        // library, though, that I'm not too worried.
        Mutex::Lock(m_conn->m_state_mutex);
        fd = amqp_get_sockfd(m_conn->m_state.get());
      }

      while (go) {
        try {
          // Waiting for frames. We don't want to hold the lock while we do that, so select here.
          // Small timeout so it shuts down fast.
          if (!select_helper(fd, 100, "select() for a method frame"))
            continue;

          SharedByteArray msg;
          {
            // Okay, we should have some data. Lock!
            Mutex::Lock(m_conn->m_state_mutex);

            // XXX: Calling maybe_release a lot keeps our memory usage down, but
            // perhaps we don't need to call it so often. Not clear on tradeoff.
            amqp_maybe_release_buffers(m_conn->m_state.get());

            amqp_frame_t method;
            vw_simple_wait_frame(m_conn->m_state.get(), &method, 1000, "Waiting for a method frame");

            // Make sure we aren't confused somehow
            VW_ASSERT(method.frame_type == AMQP_FRAME_METHOD, AMQPAssertion() 
                      << "Expected a method frame");
            VW_ASSERT(method.payload.method.id == AMQP_BASIC_DELIVER_METHOD, AMQPAssertion() 
                      << "Expected a deliver method");

	    // For debugging:
            // vw_throw(AMQPAssertion() << "Test assertion");

            // Grab the rest of the message
            msg = pump_queue(m_conn->m_state.get());
          }

          // We got some data. Release the lock and call the callback.
          m_callback(msg);
        } catch (AmqpErr &e) {
          vw_out(0) << "WARNING -- an AMQP error occurred: " << e.what() << "\n";
        }
      }
    }
};

}} // namespace vw::platefile

boost::shared_ptr<AmqpConsumer> AmqpChannel::basic_consume(std::string const& queue,
                                                           boost::function<void (SharedByteArray)> callback) {
  Mutex::Lock lock(m_conn->m_state_mutex);

  amqp_basic_consume_ok_t *reply =
    amqp_basic_consume(m_conn->m_state.get(), m_channel, amqp_string(queue), amqp_string(""), 0, 1, 0);

  die_on_amqp_error(amqp_rpc_reply, "Starting Consumer");

  boost::shared_ptr<AmqpConsumeTask> task(new AmqpConsumeTask(m_conn, m_channel, callback, queue, amqp_bytes(reply->consumer_tag)));
  boost::shared_ptr<vw::Thread> thread(new vw::Thread(task));

  return boost::shared_ptr<AmqpConsumer>( new AmqpConsumer(task, thread));
}

AmqpConsumer::AmqpConsumer(boost::shared_ptr<AmqpConsumeTask> task, boost::shared_ptr<vw::Thread> thread)
        : m_task(task), m_thread(thread) {}

AmqpConsumer::~AmqpConsumer() {
  m_task->kill();
  m_thread->join();
}

namespace {

const amqp_bytes_t amqp_string(const std::string& s) {
  amqp_bytes_t ret;
  ret.len = s.size();
  ret.bytes = const_cast<void*>(reinterpret_cast<const void*>(s.c_str()));
  return ret;
}

std::string amqp_bytes(const amqp_bytes_t& s) {
  return std::string(reinterpret_cast<char*>(s.bytes), s.len);
}

void die_on_error(int x, const std::string& context) {
  if (x < 0) {
    std::ostringstream msg;
    msg << "AMQP Error: " << context << " -- " << strerror(-x);
    vw::vw_throw(vw::IOErr() << msg.str());
  }
}

void die_on_amqp_error(amqp_rpc_reply_t x, const std::string& context) {
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

// Called with the m_state_mutex lock already held
SharedByteArray pump_queue(amqp_connection_state_t conn) {

  amqp_frame_t header;
  vw_simple_wait_frame(conn, &header, 1000, "Waiting for a header frame");

  VW_ASSERT(header.frame_type == AMQP_FRAME_HEADER, AMQPAssertion() << "Expected AMQP header!");

  size_t body_size = header.payload.properties.body_size;
  SharedByteArray payload( new ByteArray(body_size) );

  size_t body_read = 0;
  amqp_frame_t frame;

  while (body_read < body_size) {
    vw_simple_wait_frame(conn, &frame, 1000, "Waiting for a body frame");

    VW_ASSERT(frame.frame_type == AMQP_FRAME_BODY,
        AMQPAssertion() << "Expected AMQP message body!");
    VW_ASSERT(body_read + frame.payload.body_fragment.len <= body_size,
        AMQPAssertion() << "AMQP packet body size does not match header\'s body target.");

    // Copy the bytes out of the payload...
    memcpy(payload->begin() + body_read,
        frame.payload.body_fragment.bytes,
        frame.payload.body_fragment.len);

    // ... and update the number of bytes we have received
    body_read += frame.payload.body_fragment.len;
  }

  return payload;
}

bool select_helper(int fd, vw::int32 timeout, const std::string& context) {
  // Darn, out of data. Let's try to get some.
  fd_set fds;
  struct timeval tv;

  FD_ZERO(&fds);
  FD_SET(fd, &fds);
  tv.tv_sec = timeout / 1000;
  tv.tv_usec = (timeout - (tv.tv_sec * 1000)) * 1000;

  int result = select(fd+1, &fds, NULL, NULL, &tv);

  if (result == -1)
    die_on_error(errno, context + ":" "select() failed");
  else if (result == 0)
    return false;
  return true;
}

} // namespace

// These definitions are copied from amqp_private.h. They have to be in the
// global scope here. This is a terrible, evil thing to do. This should be
// removed as soon as a non-blocking solution is available from upstream.
typedef enum amqp_connection_state_enum_ {
  CONNECTION_STATE_IDLE = 0,
  CONNECTION_STATE_WAITING_FOR_HEADER,
  CONNECTION_STATE_WAITING_FOR_BODY,
  CONNECTION_STATE_WAITING_FOR_PROTOCOL_HEADER
} amqp_connection_state_enum;

typedef struct amqp_link_t_ {
      struct amqp_link_t_ *next;
        void *data;
} amqp_link_t;

struct amqp_connection_state_t_ {
  amqp_pool_t frame_pool;
  amqp_pool_t decoding_pool;

  amqp_connection_state_enum state;

  int channel_max;
  int frame_max;
  int heartbeat;
  amqp_bytes_t inbound_buffer;

  size_t inbound_offset;
  size_t target_size;

  amqp_bytes_t outbound_buffer;

  int sockfd;
  amqp_bytes_t sock_inbound_buffer;
  size_t sock_inbound_offset;
  size_t sock_inbound_limit;

  amqp_link_t *first_queued_frame;
  amqp_link_t *last_queued_frame;
};


namespace {

void vw_simple_wait_frame(amqp_connection_state_t state, amqp_frame_t *frame, vw::int32 timeout, const std::string& context) {

  // If we already have frames in the queue (perhaps processed on a previous
  // call) return one of those
  if (state->first_queued_frame != NULL) {
    amqp_frame_t *f = (amqp_frame_t *) state->first_queued_frame->data;
    state->first_queued_frame = state->first_queued_frame->next;
    if (state->first_queued_frame == NULL) {
      state->last_queued_frame = NULL;
    }
    *frame = *f;
    return;
  }

  while (1) {
    // We're out of frames. Let's get some!
    // Do we have unframed data from a previous read? Deal with it.
    while (amqp_data_in_buffer(state)) {
      amqp_bytes_t buffer;
      buffer.len = state->sock_inbound_limit - state->sock_inbound_offset;
      buffer.bytes = ((char *) state->sock_inbound_buffer.bytes) + state->sock_inbound_offset;

      // Try to frame the data
      int result = amqp_handle_input(state, buffer, frame);

      die_on_error(result, "Processing Read Frame");
      if (result == 0)
        vw_throw(AmqpErr() << "AMQP Error: EOF on socket");

      state->sock_inbound_offset += result;

      if (frame->frame_type != 0) {
        /* Complete frame was read. Return it. */
        return;
      }
    }

    // Darn, out of data. Let's try to get some.
    if (!select_helper(state->sockfd, timeout, context))
      vw_throw(AMQPTimeout() << context << ": timed out");

    // Won't block because we select()'d, and we know it's primed.
    int result = read(state->sockfd,
                      state->sock_inbound_buffer.bytes,
                      state->sock_inbound_buffer.len);

    if (result < 0)
      die_on_error(errno, context + ":" + "read() failed");
    else if (result == 0)
      vw_throw(AMQPEof() << "Socket closed");

    state->sock_inbound_limit = result;
    state->sock_inbound_offset = 0;
  }
}

} // namespace anonymous
