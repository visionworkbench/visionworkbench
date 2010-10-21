// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/AmqpConnection.h>
#include <vw/Core/Debugging.h>

#include <cstdlib>
#include <cstdio>
#include <cerrno>

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
  SharedByteArray read_content(amqp_connection_state_t conn);
  bool select_helper(int fd, vw::int32 timeout, const std::string& context);
  bool vw_simple_wait_frame(amqp_connection_state_t conn, amqp_frame_t *frame,
                            vw::int32 timeout, const std::string& context);
}

AmqpConnection::AmqpConnection(std::string const& hostname, int port) {

  m_state.reset(amqp_new_connection(), amqp_destroy_connection);
  if (!m_state.get())
    vw_throw(AMQPErr() << "Failed to create amqp state object");

  // Open a socket and establish an amqp connection
  int fd = amqp_open_socket(hostname.c_str(), port);
  if (fd < 0)
    vw_throw(AMQPErr() << "Failed to open AMQP socket.");

  int flag = 1;
  die_on_error(setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &flag, sizeof(int)), "setting TCP_NODELAY");

  amqp_set_sockfd(m_state.get(), fd);

  // Login to the Amqp Server
  die_on_amqp_error(amqp_login(m_state.get(), "/", 0, 131072, 0, AMQP_SASL_METHOD_PLAIN,
                               "guest", "guest"), "logging in");

  // 0 is reserved
  m_used_channels.insert(0);
}

AmqpConnection::~AmqpConnection() {
  Mutex::Lock lock(m_state_mutex);

  try {
    die_on_amqp_error(amqp_connection_close(m_state.get(), AMQP_REPLY_SUCCESS), "closing connection");
    die_on_error(close(amqp_get_sockfd(m_state.get())), "closing socket");
  } catch (const AMQPErr& e) {
    vw_out(ErrorMessage, "plate.AMQP")
      << "Caught AMQPErr in " << VW_CURRENT_FUNCTION << ": "
      << e.what() << std::endl;
  }
}

vw::Mutex& AmqpConnection::get_mutex(AmqpConnectionState** state) {
  *state = m_state.get();
  return m_state_mutex;
}

// CALL THIS WITH THE m_state_mutex LOCK ALREADY HELD.
int16 AmqpConnection::get_channel(int16 channel) {
  if (channel == -1) {
    if(m_used_channels.empty())
      channel = 1;
    else
      channel = boost::numeric_cast<int16>(*m_used_channels.rbegin() + 1); // largest elt in set
  }

  if (m_used_channels.count(channel) != 0)
    vw_throw(LogicErr() << "Channel " << channel << "is already in use");

  m_used_channels.insert(channel);
  return channel;
}


#define ASSERT_CHANNEL_OPEN() do { \
  if (!is_open)\
    vw_throw(AMQPChannelErr() << "Channel is already closed due to a previous channel error");\
} while(0)

void AmqpChannel::check_error(amqp_rpc_reply_t x, const std::string& context) {
  if (x.reply_type == AMQP_RESPONSE_SERVER_EXCEPTION &&
      x.reply.id   == AMQP_CHANNEL_CLOSE_METHOD)
    is_open = false;
  die_on_amqp_error(x, context);
}

AmqpChannel::AmqpChannel(boost::shared_ptr<AmqpConnection> conn, int16 channel)
    : m_conn(conn), is_open(false) {

  amqp_connection_state_t state;

  Mutex::Lock lock(m_conn->get_mutex(&state));
  m_channel = m_conn->get_channel(channel);

  amqp_channel_open(state, m_channel);
  check_error(amqp_get_rpc_reply(state), "opening channel");
  is_open = true;
}

AmqpChannel::~AmqpChannel() {
  amqp_connection_state_t state;
  Mutex::Lock lock(m_conn->get_mutex(&state));
  if (is_open) {
    try {
      die_on_amqp_error(amqp_channel_close(state, m_channel, AMQP_REPLY_SUCCESS), "closing connection");
    } catch (const AMQPErr& e) {
      vw_out(ErrorMessage, "plate.AMQP")
        << "Caught AMQPErr in " << VW_CURRENT_FUNCTION << ": "
        << e.what() << std::endl;
    }
  }
}

void AmqpChannel::exchange_declare(std::string const& exchange_name,
                                   std::string const& exchange_type,
                                   bool durable, bool auto_delete) {
  amqp_connection_state_t state;
  Mutex::Lock lock(m_conn->get_mutex(&state));
  ASSERT_CHANNEL_OPEN();

  amqp_exchange_declare(state, m_channel, amqp_string(exchange_name),
                        amqp_string(exchange_type), 0, durable, auto_delete, amqp_table_t());

  check_error(amqp_get_rpc_reply(state), "declaring exchange");
}

void AmqpChannel::queue_declare(std::string const& queue_name, bool durable,
                                bool exclusive, bool auto_delete) {
  amqp_connection_state_t state;
  Mutex::Lock lock(m_conn->get_mutex(&state));
  ASSERT_CHANNEL_OPEN();

  amqp_queue_declare(state, m_channel, amqp_string(queue_name), 0,
                     durable, exclusive, auto_delete, amqp_table_t());

  check_error(amqp_get_rpc_reply(state), "declaring queue");
}

void AmqpChannel::queue_bind(std::string const& queue, std::string const& exchange,
                             std::string const& routing_key) {
  amqp_connection_state_t state;
  Mutex::Lock lock(m_conn->get_mutex(&state));
  ASSERT_CHANNEL_OPEN();

  amqp_queue_bind(state, m_channel, amqp_string(queue),
                  amqp_string(exchange), amqp_string(routing_key), amqp_table_t());

  check_error(amqp_get_rpc_reply(state), "binding queue");
}

void AmqpChannel::queue_unbind(std::string const& queue, std::string const& exchange,
                               std::string const& routing_key) {
  amqp_connection_state_t state;
  Mutex::Lock lock(m_conn->get_mutex(&state));
  ASSERT_CHANNEL_OPEN();

  amqp_queue_unbind(state, m_channel, amqp_string(queue),
                    amqp_string(exchange), amqp_string(routing_key), amqp_table_t());

  check_error(amqp_get_rpc_reply(state), "unbinding queue");
}

void AmqpChannel::basic_publish(ByteArray const& message,
                                std::string const& exchange, std::string const& routing_key) {
  amqp_connection_state_t state;
  Mutex::Lock lock(m_conn->get_mutex(&state));
  ASSERT_CHANNEL_OPEN();

  amqp_bytes_t raw_data;
  raw_data.len   = message.size();
  raw_data.bytes = const_cast<void*>(reinterpret_cast<const void*>(message.begin()));

  int ret = amqp_basic_publish(state, m_channel, amqp_string(exchange),
                               amqp_string(routing_key), 0, 0, NULL, raw_data);
  die_on_error(ret, "doing a basic.publish");
}

bool AmqpChannel::basic_get(std::string const& queue, SharedByteArray& message) {

  amqp_connection_state_t state;
  Mutex::Lock lock(m_conn->get_mutex(&state));
  ASSERT_CHANNEL_OPEN();

  amqp_rpc_reply_t reply;
  while (true) {
    amqp_maybe_release_buffers(state);
    reply = amqp_basic_get(state, m_channel, amqp_string(queue), 1);
    check_error(reply, "doing a basic.get");

    if (reply.reply.id == AMQP_BASIC_GET_EMPTY_METHOD) {
      usleep(100000); // yield for a bit
      continue;
    }
    if (reply.reply.id == AMQP_BASIC_GET_OK_METHOD)
      break;

    vw_throw(AMQPAssertion() << "Illegal AMQP response. Expected GET_OK or GET_EMPTY, got: "
                             << amqp_method_name(reply.reply.id));
  }

  message = read_content(state);
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
        amqp_connection_state_t state;
        Mutex::Lock lock(m_conn->get_mutex(&state));
        fd = amqp_get_sockfd(state);
      }

      while (go) {
        // Waiting for frames. We don't want to hold the lock while we do that, so select here.
        // Small timeout so it shuts down fast.
        if (!select_helper(fd, 100, "select() for a method frame"))
          continue;

        amqp_frame_t method;

        SharedByteArray msg;
        {
          // Okay, we should have some data. Lock!
          amqp_connection_state_t state;
          Mutex::Lock lock(m_conn->get_mutex(&state));

          // XXX: Calling maybe_release a lot keeps our memory usage down, but
          // perhaps we don't need to call it so often. Not clear on tradeoff.
          amqp_maybe_release_buffers(state);

          while (1) {
            if (!go) return;
            if (vw_simple_wait_frame(state, &method, 1000, "Waiting for a method frame"))
              break;
          }

          // Make sure we aren't confused somehow
          VW_ASSERT(method.frame_type == AMQP_FRAME_METHOD,
              AMQPAssertion() << "Expected a method frame, got: " << method.frame_type);

          // Grab the rest of the message (if there is any)
          if (amqp_method_has_content(method.payload.method.id))
            msg = read_content(state);
        }

        // We got some data. Release the lock and call the callback if it was a delivery.
        if (method.payload.method.id == AMQP_BASIC_DELIVER_METHOD)
          m_callback(msg);
        else {
          vw_out(WarningMessage, "plate.amqp") << "Dropped "
            << amqp_method_name(method.payload.method.id) << "on the floor."
            << std::endl;
        }

      }
    }
};

}} // namespace vw::platefile

boost::shared_ptr<AmqpConsumer> AmqpChannel::basic_consume(std::string const& queue,
                                                           boost::function<void (SharedByteArray)> callback) {
  amqp_connection_state_t state;
  Mutex::Lock lock(m_conn->get_mutex(&state));

  amqp_basic_consume_ok_t *reply =
    amqp_basic_consume(state, m_channel, amqp_string(queue), amqp_string(""), 0, 1, 0);

  check_error(amqp_get_rpc_reply(state), "starting consumer");

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
    vw_throw(AMQPErr() << strerror(-x) << " while " << context);
  }
}

void die_on_amqp_error(const amqp_rpc_reply_t x, const std::string& context) {
  switch (x.reply_type) {
    case AMQP_RESPONSE_NORMAL: return;
    case AMQP_RESPONSE_NONE:
      vw_throw(AMQPErr() << "missing RPC reply type while " << context);
    case AMQP_RESPONSE_LIBRARY_EXCEPTION: {
      // amqp_error_string returns a malloc'd string. use free as a custom deleter
      boost::shared_ptr<char> error(amqp_error_string(x.library_error), free);
      vw_throw(AMQPErr() << *error << " while " << context);
    }
    case AMQP_RESPONSE_SERVER_EXCEPTION:
      switch (x.reply.id) {
        case AMQP_CONNECTION_CLOSE_METHOD: {
          amqp_connection_close_t *m = (amqp_connection_close_t *) x.reply.decoded;
          vw_throw(AMQPConnectionErr()
                     << static_cast<const char*>(m->reply_text.bytes)
                     << " while " << context);
        }
        case AMQP_CHANNEL_CLOSE_METHOD: {
          amqp_channel_close_t *m = (amqp_channel_close_t *) x.reply.decoded;
          vw_throw(AMQPChannelErr()
                     << static_cast<const char*>(m->reply_text.bytes)
                     << " while " << context);
        }
        default:
          vw_throw(AMQPErr() << "unknown server error [id " << x.reply.id
                             << "] while " << context);
    }
    default:
      vw_throw(AMQPErr() << "AMQP Error: unknown response type while "
                         << context);
  }
}

// Called with the m_state_mutex lock already held
SharedByteArray read_content(amqp_connection_state_t conn) {

  amqp_frame_t header;
  vw_simple_wait_frame(conn, &header, 1000, "Waiting for a header frame");

  VW_ASSERT(header.frame_type == AMQP_FRAME_HEADER,
      AMQPAssertion() << "Expected AMQP header, got: " << header.frame_type);

  size_t body_size = header.payload.properties.body_size;
  SharedByteArray payload( new ByteArray(body_size) );

  size_t body_read = 0;
  amqp_frame_t frame;

  while (body_read < body_size) {
    vw_simple_wait_frame(conn, &frame, 1000, "Waiting for a body frame");

    VW_ASSERT(frame.frame_type == AMQP_FRAME_BODY,
        AMQPAssertion() << "Expected AMQP message body, got: " << frame.frame_type);
    VW_ASSERT(body_read + frame.payload.body_fragment.len <= body_size,
        AMQPAssertion() << "AMQP packet body size does not match header's body target.");

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
    die_on_error(-errno, context + ":" "select() failed");
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

bool vw_simple_wait_frame(amqp_connection_state_t state, amqp_frame_t *frame, vw::int32 timeout, const std::string& context) {

  // If we already have frames in the queue (perhaps processed on a previous
  // call) return one of those
  if (state->first_queued_frame != NULL) {
    amqp_frame_t *f = (amqp_frame_t *) state->first_queued_frame->data;
    state->first_queued_frame = state->first_queued_frame->next;
    if (state->first_queued_frame == NULL) {
      state->last_queued_frame = NULL;
    }
    *frame = *f;
    return true;
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

      die_on_error(result, "processing read frame");
      if (result == 0)
        vw_throw(AMQPErr() << "EOF on socket");

      state->sock_inbound_offset += result;

      if (frame->frame_type != 0) {
        /* Complete frame was read. Return it. */
        return true;
      }
    }

    // Darn, out of data. Let's try to get some.
    if (!select_helper(state->sockfd, timeout, context))
      return false;

    // Won't block because we select()'d, and we know it's primed.
    ssize_t result = read(state->sockfd,
                          state->sock_inbound_buffer.bytes,
                          state->sock_inbound_buffer.len);

    if (result < 0)
      die_on_error(-errno, context + ":" + "read() failed");
    else if (result == 0)
      vw_throw(AMQPEof() << "Socket closed");

    state->sock_inbound_limit = result;
    state->sock_inbound_offset = 0;
  }
}

} // namespace anonymous
