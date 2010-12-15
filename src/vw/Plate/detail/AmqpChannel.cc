// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/detail/AmqpChannel.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Core/Debugging.h>
#include <vw/Core/Log.h>
#include <vw/Core/Thread.h>
#include <vw/Core/ThreadQueue.h>
#include <vw/config.h>
#include <vw/Plate/FundamentalTypes.h>

#include <google/protobuf/descriptor.h>

#include <amqp.h>
#include <amqp_framing.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/ip.h>
#include <netinet/tcp.h>

#include <set>

using namespace vw;
using namespace vw::platefile;
namespace pb = ::google::protobuf;

namespace {

  void check_error(const amqp_rpc_reply_t& x, const std::string& context);
  const amqp_bytes_t amqp_string(const std::string& s);
  std::string amqp_bytes(const amqp_bytes_t& s);
  void die_on_amqp_error(amqp_rpc_reply_t x, const std::string& context);
  void die_on_error(int x, const std::string& context);
  void read_content(amqp_connection_state_t conn, AmqpData& msg);
  bool select_helper(int fd, int32 timeout, const std::string& context);
  bool vw_simple_wait_frame(amqp_connection_state_t conn, amqp_frame_t *frame,
                            int32 timeout, const std::string& context);
}

namespace vw {
namespace platefile {

// This is an implementation detail.
typedef boost::remove_pointer<amqp_connection_state_t>::type AmqpConnectionState;
class AmqpConnection : private boost::noncopyable {
  private:
    // Any access (read or write) to m_state needs to hold the m_state_mutex lock
    boost::shared_ptr<AmqpConnectionState> m_state;
    Mutex m_state_mutex;
    std::set<int16> m_used_channels;

  public:
    /// Open a new connection to the AMQP server.  This connection
    /// terminates automatically when this object is destroyed.
    AmqpConnection(std::string const& hostname = "localhost", int port = 5672);

    Mutex& get_mutex(AmqpConnectionState** state) VW_WARN_UNUSED;
    void locked_get_state(AmqpConnectionState** state);

    /// Allocates a communications channel number. Meant to be called by
    ///     AmqpChannel's constructor. Grab the mutex BEFORE you call this.
    /// @param channel Request a specific channel number. -1 for don't care.
    //                 Requesting a particular one will throw if it's already used.
    int16 locked_get_channel(int16 channel = -1);


    /// Closes the AMQP connection and destroys this object.
    ~AmqpConnection();
};

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
  amqp_connection_state_t state;
  Mutex::Lock lock(get_mutex(&state));

  try {
    die_on_amqp_error(amqp_connection_close(state, AMQP_REPLY_SUCCESS), "closing connection");
    die_on_error(close(amqp_get_sockfd(state)), "closing socket");
  } catch (const AMQPErr& e) {
    vw_out(ErrorMessage, "plate.AMQP")
      << "Caught AMQPErr in " << VW_CURRENT_FUNCTION << ": "
      << e.what() << std::endl;
  }
}

Mutex& AmqpConnection::get_mutex(AmqpConnectionState** state) {
  *state = m_state.get();
  return m_state_mutex;
}

void AmqpConnection::locked_get_state(AmqpConnectionState** state) {
  *state = m_state.get();
}

int16 AmqpConnection::locked_get_channel(int16 channel) {
  VW_ASSERT(channel >= -1 && channel < 9 && m_used_channels.size() < 9, LogicErr() << "Only 8 channels available");
  if (channel == -1) {
    if(m_used_channels.empty())
      channel = 1;
    else
      channel = static_cast<int16>(*m_used_channels.rbegin() + 1); // largest elt in set
  }

  if (m_used_channels.count(channel) != 0)
    vw_throw(LogicErr() << "Channel " << channel << "is already in use");

  m_used_channels.insert(channel);
  return channel;
}


#define ASSERT_CHANNEL_OPEN() do { \
  if (!m_conn)\
    vw_throw(AMQPErr() << "Channel is already closed due to a previous channel error");\
} while(0)

void AmqpChannel::create_endpoint(const std::string& rabbitmq, short port, const std::string& name) {
  VW_ASSERT(!m_conn, LogicErr() << "please do not reuse AmqpChannels. create a new one instead");
  VW_ASSERT(!name.empty() && name[0] == '/', ArgumentErr() << "endpoint name should look like a pathname (got " << name << ")");

  vw_out(VerboseDebugMessage, "plate.AMQP") << m_human_name << " creating endpoint " << name << std::endl;
  m_conn.reset(new AmqpConnection(rabbitmq, port));

  amqp_connection_state_t state;
  Mutex::Lock lock(m_conn->get_mutex(&state));
  m_channel_id = m_conn->locked_get_channel();

  amqp_channel_open(state, m_channel_id);
  check_error(amqp_get_rpc_reply(state), "opening channel");

  m_local_name = name;

  // declare the local side of the pipe
  this->locked_queue_declare(m_local_name, false, true, true);
  this->locked_exchange_declare(m_local_name, "fanout", false, true);
  this->locked_queue_bind(m_local_name, m_local_name, m_local_name);

  m_consumer = this->locked_basic_consume(m_local_name,
      boost::bind(&ThreadQueue<AmqpData>::push, boost::ref(m_incoming_messages), _1));
}

void AmqpChannel::remove_endpoint() {
  if (!m_conn)
    return;

  m_consumer.reset();

  if (!m_local_name.empty()) {
    amqp_connection_state_t state;
    Mutex::Lock lock(m_conn->get_mutex(&state));
    this->locked_queue_unbind(m_local_name, m_local_name, m_local_name);
    die_on_amqp_error(amqp_channel_close(state, m_channel_id, AMQP_REPLY_SUCCESS), "closing connection");
    m_channel_id = -1;
  }
  m_conn.reset();
}

void AmqpChannel::conn(const Url& endpoint) {
  const std::string host = endpoint.hostname().empty() ? "localhost" : endpoint.hostname();
  const short port       = short(endpoint.port() == 0  ? 5672        : endpoint.port());

  this->create_endpoint(host, port, unique_name(std::string("/") + m_human_name, "/"));

  vw_out(VerboseDebugMessage, "plate.AMQP") << m_human_name << " connecting to " << endpoint.path() << std::endl;

  // TODO: Make sure the server side exists, but don't declare it
  m_remote_name = endpoint.path();
}

void AmqpChannel::bind(const Url& endpoint) {
  const std::string host = endpoint.hostname().empty() ? "localhost" : endpoint.hostname();
  const short port       = short(endpoint.port() == 0  ? 5672        : endpoint.port());

  this->create_endpoint(host, port, endpoint.path());
  m_remote_name = "";
}


AmqpChannel::~AmqpChannel() VW_NOTHROW {
  try {
    remove_endpoint();
  } catch (const AMQPErr& e) {
    vw_out(ErrorMessage, "plate.AMQP") << "Failed to close amqp channel: " << e.what() << std::endl;
  }
}

void AmqpChannel::send_bytes(const uint8* message, size_t len) {
  VW_ASSERT(!m_remote_name.empty(), LogicErr() << m_human_name << ": Don't know who to send_bytes to.");
  vw_out(VerboseDebugMessage, "plate.AMQP") << m_human_name << " sending data to " << m_remote_name << std::endl;
  this->basic_publish(message, len, m_remote_name, m_remote_name);
}

bool AmqpChannel::recv_bytes(std::vector<uint8>* bytes) {
  vw_out(VerboseDebugMessage, "plate.AMQP") << m_human_name << " waiting for data for " << this->timeout() << "ms" << std::endl;
  AmqpData d;
  if (this->timeout() == -1)
    this->m_incoming_messages.wait_pop(d);
  else {
    if (!this->m_incoming_messages.timed_wait_pop(d, this->timeout()))
      return false;
  }

  bytes->clear();
  bytes->reserve(d.data->size());
  std::copy(d.data->begin(), d.data->end(), std::back_inserter(*bytes));
  m_remote_name = d.sender;
  return true;
}

void AmqpChannel::CallMethod(const pb::MethodDescriptor* method,
                             pb::RpcController* /*controller*/,
                             const pb::Message* request,
                             pb::Message* response,
                             pb::Closure* done)
{
  detail::RequireCall call(done);
  RpcWrapper q_wrap, a_wrap;

  q_wrap.set_method(method->name());
  q_wrap.set_payload(request->SerializeAsString());
  q_wrap.set_requestor(this->name());

  for (uint32 trial = 0; trial <= m_retries; ++trial) {
    if (trial > 0)
      vw_out(WarningMessage) << "Retry (" << trial << "/" << m_retries << ")" << std::endl;
    q_wrap.set_seq(++m_seq);

    send_message(q_wrap);

CallMethod_receive_again:
    switch (recv_message(a_wrap)) {
      case 0:
        vw_out(WarningMessage) << "CallMethod Timeout. ";
        continue;
      case -1:
        vw_out(WarningMessage) << "CallMethod(): corrupted message. ";
        continue;
      default:
        if (a_wrap.seq() != q_wrap.seq()) {
          vw_out(WarningMessage) << "Sequence mismatch on \"" << this->name() << "\" (expected " << m_seq << ", got " << a_wrap.seq() << ") ";
          // If we get a seq less than we were expecting, one of the messages
          // we timed out waiting for finally arrived. Drop it on the floor and
          // wait again rather than retrying completely.
          if (a_wrap.seq() < q_wrap.seq())
            goto CallMethod_receive_again;
          continue;
        }
        throw_rpc_error(a_wrap.error());
        response->ParseFromString(a_wrap.payload());
        return;
    }
  }
  vw_out(WarningMessage) << "No more retries." << std::endl;
  vw_throw(RpcErr() << "CallMethod timed out completely");
}

void AmqpChannel::locked_exchange_declare(std::string const& exchange_name,
                                          std::string const& exchange_type,
                                          bool durable, bool auto_delete) {
  ASSERT_CHANNEL_OPEN();
  amqp_connection_state_t state;
  m_conn->locked_get_state(&state);

  amqp_exchange_declare(state, m_channel_id, amqp_string(exchange_name),
                        amqp_string(exchange_type), 0, durable, auto_delete, amqp_table_t());
  check_error(amqp_get_rpc_reply(state), "declaring exchange");
}

void AmqpChannel::locked_queue_declare(std::string const& queue_name, bool durable,
                                       bool exclusive, bool auto_delete) {
  ASSERT_CHANNEL_OPEN();
  amqp_connection_state_t state;
  m_conn->locked_get_state(&state);

  amqp_queue_declare(state, m_channel_id, amqp_string(queue_name), 0,
                     durable, exclusive, auto_delete, amqp_table_t());

  check_error(amqp_get_rpc_reply(state), "declaring queue");
}

void AmqpChannel::locked_queue_bind(std::string const& queue, std::string const& exchange,
                                    std::string const& routing_key) {
  ASSERT_CHANNEL_OPEN();
  amqp_connection_state_t state;
  m_conn->locked_get_state(&state);

  amqp_queue_bind(state, m_channel_id, amqp_string(queue),
                  amqp_string(exchange), amqp_string(routing_key), amqp_table_t());

  check_error(amqp_get_rpc_reply(state), "binding queue");
}

void AmqpChannel::locked_queue_unbind(std::string const& queue, std::string const& exchange,
                                      std::string const& routing_key) {
  ASSERT_CHANNEL_OPEN();
  amqp_connection_state_t state;
  m_conn->locked_get_state(&state);

  amqp_queue_unbind(state, m_channel_id, amqp_string(queue),
                    amqp_string(exchange), amqp_string(routing_key), amqp_table_t());

  check_error(amqp_get_rpc_reply(state), "unbinding queue");
}

void AmqpChannel::basic_publish(const uint8* message, uint64 len,
                                std::string const& exchange, std::string const& routing_key) const {
  ASSERT_CHANNEL_OPEN();
  amqp_connection_state_t state;
  Mutex::Lock lock(m_conn->get_mutex(&state));

  amqp_bytes_t raw_data;
  raw_data.len   = len;
  raw_data.bytes = const_cast<void*>(reinterpret_cast<const void*>(message));

  amqp_basic_properties_t props;
  props._flags = AMQP_BASIC_REPLY_TO_FLAG;
  props.reply_to   = amqp_string(m_local_name);
  // props._flags |= |AMQP_MESSAGE_ID_FLAG;
  // props.message_id =

  int ret = amqp_basic_publish(state, m_channel_id, amqp_string(exchange),
                               amqp_string(routing_key), 0, 0, &props, raw_data);
  die_on_error(ret, "doing a basic.publish");
}

#if 0
bool AmqpChannel::basic_get(std::string const& queue, SharedByteArray& message) const {

  ASSERT_CHANNEL_OPEN();
  amqp_connection_state_t state;
  Mutex::Lock lock(m_conn->get_mutex(&state));

  amqp_rpc_reply_t reply;
  while (true) {
    amqp_maybe_release_buffers(state);
    reply = amqp_basic_get(state, m_channel_id, amqp_string(queue), 1);
    check_error(reply, "doing a basic.get");

    if (reply.reply.id == AMQP_BASIC_GET_EMPTY_METHOD) {
      Thread::sleep_ms(1); // yield for a bit
      continue;
    }
    if (reply.reply.id == AMQP_BASIC_GET_OK_METHOD)
      break;

    vw_throw(AMQPErr() << "Illegal AMQP response. Expected GET_OK or GET_EMPTY, got: "
                       << amqp_method_name(reply.reply.id));
  }

  AmqpData msg;
  read_content(state, msg);
  message = msg.data;
  return true;
}
#endif

int32 AmqpChannel::timeout() const {
  return m_timeout;
}

void AmqpChannel::set_timeout(int32 x) {
  m_timeout = x;
}

uint32 AmqpChannel::retries() const {
  return m_retries;
}

void AmqpChannel::set_retries(uint32 x) {
  m_retries = x;
}

#if 0
uint64 AmqpChannel::queue_depth() const {
  return m_incoming_messages.size();
}
#endif

std::string AmqpChannel::name() const {
  return m_human_name;
}

class AmqpConsumeTask {
  public:
    typedef boost::shared_ptr<AmqpConnection> Connection;
    typedef boost::function<void (AmqpData)> Callback;

    // XXX: This should really hold a channel, and not a connection...
    Connection m_conn;
    int16 m_chan;
    Callback m_callback;
    std::string m_queue;
    std::string m_consumer_tag;
    bool go;
    int m_fd;

    AmqpConsumeTask(Connection conn, int16 channel_id, const Callback& callback,
                    std::string queue_name, std::string consumer_tag)
        : m_conn(conn), m_chan(channel_id), m_callback(callback), m_queue(queue_name),
          m_consumer_tag(consumer_tag), go(true) {

        // This isn't totally safe, because technically the socket could be
        // changed from under us. That would break so much of the rest of the
        // library, though, that I'm not too worried.
        amqp_connection_state_t state;
        m_conn->locked_get_state(&state);
        m_fd = amqp_get_sockfd(state);
    }

    void kill() {go = false;}

    void operator()() const {
      while (go) {
        // Waiting for frames. We don't want to hold the lock while we do that, so select here.
        // Small timeout so it shuts down fast.
        if (!select_helper(m_fd, 100, "select() for a method frame"))
          continue;

        amqp_frame_t method;

        AmqpData msg;
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
              AMQPErr() << "Expected a method frame, got: " << method.frame_type);

          // Grab the rest of the message (if there is any)
          if (amqp_method_has_content(method.payload.method.id))
            read_content(state, msg);
        }

        switch (method.payload.method.id) {
          case AMQP_BASIC_DELIVER_METHOD:
            m_callback(msg);
            break;
          case AMQP_CONNECTION_CLOSE_METHOD: {
            amqp_connection_close_t *m = reinterpret_cast<amqp_connection_close_t *>(method.payload.method.decoded);
            vw_throw(AMQPErr() << "Unexpected connection close method: " << reinterpret_cast<const char*>(m->reply_text.bytes));
          }
          case AMQP_CHANNEL_CLOSE_METHOD: {
            amqp_channel_close_t *m = reinterpret_cast<amqp_channel_close_t *>(method.payload.method.decoded);
            vw_throw(AMQPErr() << "Unexpected channel close method: " << reinterpret_cast<const char*>(m->reply_text.bytes));
          }
          default:
            vw_out(WarningMessage, "plate.AMQP")
              << "Dropped " << amqp_method_name(method.payload.method.id) << " on the floor."
              << std::endl;
        }
      }
    }
};

boost::shared_ptr<AmqpConsumer>
AmqpChannel::locked_basic_consume(std::string const& queue, boost::function<void (AmqpData)> callback) {

  ASSERT_CHANNEL_OPEN();
  amqp_connection_state_t state;
  m_conn->locked_get_state(&state);

  amqp_basic_consume_ok_t *reply =
    amqp_basic_consume(state, m_channel_id, amqp_string(queue), amqp_string(""), 0, 1, 0);

  check_error(amqp_get_rpc_reply(state), "starting consumer");

  boost::shared_ptr<AmqpConsumeTask> task(new AmqpConsumeTask(m_conn, m_channel_id, callback, queue, amqp_bytes(reply->consumer_tag)));
  boost::shared_ptr<Thread> thread(new Thread(task));

  return boost::shared_ptr<AmqpConsumer>( new AmqpConsumer(task, thread));
}

AmqpConsumer::AmqpConsumer(boost::shared_ptr<AmqpConsumeTask> task, boost::shared_ptr<Thread> thread)
        : m_task(task), m_thread(thread) {}

AmqpConsumer::~AmqpConsumer() {
  m_task->kill();
  m_thread->join();
}

}} // namespace vw::platefile

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
          vw_throw(AMQPErr()
                     << static_cast<const char*>(m->reply_text.bytes)
                     << " while " << context);
        }
        case AMQP_CHANNEL_CLOSE_METHOD: {
          amqp_channel_close_t *m = (amqp_channel_close_t *) x.reply.decoded;
          vw_throw(AMQPErr()
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
void read_content(amqp_connection_state_t conn, AmqpData& msg) {

  amqp_frame_t header;
  vw_simple_wait_frame(conn, &header, 1000, "Waiting for a header frame");

  VW_ASSERT(header.frame_type == AMQP_FRAME_HEADER,
      AMQPErr() << "Expected AMQP header, got: " << header.frame_type);

  amqp_basic_properties_t *p = reinterpret_cast<amqp_basic_properties_t*>(header.payload.properties.decoded);

  VW_ASSERT(p->_flags & AMQP_BASIC_REPLY_TO_FLAG, LogicErr() << "Got a message without a reply-to property");
  msg.sender = amqp_bytes(p->reply_to);

  size_t body_size = header.payload.properties.body_size;
  msg.data.reset(new std::vector<uint8>(body_size) );

  size_t body_read = 0;
  amqp_frame_t frame;

  while (body_read < body_size) {
    vw_simple_wait_frame(conn, &frame, 1000, "Waiting for a body frame");

    VW_ASSERT(frame.frame_type == AMQP_FRAME_BODY,
        AMQPErr() << "Expected AMQP message body, got: " << frame.frame_type);
    VW_ASSERT(body_read + frame.payload.body_fragment.len <= body_size,
        AMQPErr() << "AMQP packet body size does not match header's body target.");

    // Copy the bytes out of the payload...
    memcpy(&(msg.data->operator[](0)) + body_read,
        frame.payload.body_fragment.bytes,
        frame.payload.body_fragment.len);

    // ... and update the number of bytes we have received
    body_read += frame.payload.body_fragment.len;
  }
}

bool select_helper(int fd, int32 timeout, const std::string& context) {
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

void check_error(const amqp_rpc_reply_t& x, const std::string& context) {
  if (x.reply_type == AMQP_RESPONSE_SERVER_EXCEPTION &&
      x.reply.id   == AMQP_CHANNEL_CLOSE_METHOD)
  die_on_amqp_error(x, context);
}


bool vw_simple_wait_frame(amqp_connection_state_t state, amqp_frame_t *frame, int32 timeout, const std::string& context) {

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
      vw_throw(AMQPErr() << "Socket closed unexpectedly");

    state->sock_inbound_limit = result;
    state->sock_inbound_offset = 0;
  }
}

} // namespace anonymous
