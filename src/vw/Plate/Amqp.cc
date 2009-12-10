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
  void die_on_amqp_error(amqp_rpc_reply_t x, char const *context);
  void die_on_error(int x, char const *context);
  SharedByteArray pump_queue(amqp_connection_state_t conn);
  void set_nonblock(int fd, bool yes);
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
        // XXX: This is so evil. Turning on nonblock here prevents rabbitmq-c
        // from failing on rpc, but will break if someone tries to start
        // another capture thread
        Mutex::Lock lock(m_conn->m_state_mutex);
        fd = amqp_get_sockfd(m_conn->m_state.get());
      }

      fd_set fds;
      struct timeval tv;
      int ret;
      amqp_frame_t method;
      SharedByteArray msg;

      while(go) {
        FD_ZERO(&fds);
        FD_SET(fd, &fds);
        tv.tv_sec = 0;
        tv.tv_usec = 100000;

        ret = select(fd+1, &fds, NULL, NULL, &tv);

        if (ret == -1)
          vw_throw(IOErr() << "select() failed: " << strerror(errno));
        else if (ret == 0)
          continue; // timeout
        else {
          Mutex::Lock lock(m_conn->m_state_mutex);
          set_nonblock(fd, true);

          amqp_maybe_release_buffers(m_conn->m_state.get());
          ret = amqp_simple_wait_frame(m_conn->m_state.get(), &method);
          if (ret == EAGAIN || ret == EWOULDBLOCK)
            continue;

          VW_ASSERT(method.frame_type == AMQP_FRAME_METHOD, AMQPAssertion() << "Expected a method frame");
          VW_ASSERT(method.payload.method.id == AMQP_BASIC_DELIVER_METHOD, AMQPAssertion() << "Expected a deliver method");

          msg = pump_queue(m_conn->m_state.get());
          set_nonblock(fd, false);
        }
        m_callback(msg);
      }

      {
        Mutex::Lock lock(m_conn->m_state_mutex);

        int flags = fcntl(fd,F_GETFL,0);
        if (flags == -1)
          vw_throw(IOErr() << "Couldn't read socket flags");

        if (fcntl(fd, F_SETFL, flags & (~O_NONBLOCK)) < 0)
          vw_throw(IOErr() << "Couldn't make amqp socket blocking again");
      }

    }
};

}} // namespace vw::platefile

boost::shared_ptr<AmqpConsumer> AmqpChannel::basic_consume(std::string const& queue,
                                                           boost::function<void (SharedByteArray)> callback) {
  Mutex::Lock lock(m_conn->m_state_mutex);

  amqp_basic_consume_ok_t *reply = amqp_basic_consume(m_conn->m_state.get(), m_channel, amqp_string(queue), amqp_string(""), 0, 1, 0);
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

SharedByteArray pump_queue(amqp_connection_state_t conn) {

  amqp_frame_t header;

  amqp_maybe_release_buffers(conn);
  int result = amqp_simple_wait_frame(conn, &header);

  VW_ASSERT(result > 0,                             AMQPAssertion() << "AMQP error: unknown result waiting for frame.");
  VW_ASSERT(header.frame_type == AMQP_FRAME_HEADER, AMQPAssertion() << "Expected AMQP header!");

  size_t body_size = header.payload.properties.body_size;
  SharedByteArray payload( new ByteArray(body_size) );

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
    memcpy(payload->begin() + body_read,
        frame.payload.body_fragment.bytes,
        frame.payload.body_fragment.len);

    // ... and update the number of bytes we have received
    body_read += frame.payload.body_fragment.len;
  }

  return payload;
}

void set_nonblock(int fd, bool yes) {
  int flags = fcntl(fd,F_GETFL,0);
  if (flags == -1)
    vw_throw(IOErr() << "Couldn't read socket flags");

  if (yes)
    flags |= O_NONBLOCK;
  else
    flags &= ~O_NONBLOCK;

  if (fcntl(fd, F_SETFL, flags) < 0)
    vw_throw(IOErr() << "Couldn't make amqp socket non-blocking");
}

} // namespace anonymous
