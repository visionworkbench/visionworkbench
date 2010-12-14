// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/detail/ZeroMQChannel.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Plate/FundamentalTypes.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Core/Log.h>

#include <boost/algorithm/string/trim.hpp>
#include <google/protobuf/descriptor.h>
#include <zmq.hpp>
#include <cerrno>

using namespace vw;
using namespace vw::platefile;
namespace pb = ::google::protobuf;

static const uint8_t THREAD_COUNT = 10;

namespace {
  vw::RunOnce zmq_init_once = VW_RUNONCE_INIT;
  boost::shared_ptr<zmq::context_t> zmq_ctx_ptr;
  void init_zmq() {
    zmq_ctx_ptr.reset(new zmq::context_t(THREAD_COUNT));
  }

  boost::shared_ptr<zmq::context_t> get_ctx() {
    zmq_init_once.run(init_zmq);
    return zmq_ctx_ptr;
  }

  const char* zmq_err() {
    return zmq_strerror(errno);
  }
}

ZeroMQChannel::ZeroMQChannel(const std::string& human_name)
  : m_ctx(get_ctx()), m_human_name(human_name), m_timeout(DEFAULT_TIMEOUT), m_retries(DEFAULT_RETRIES) {}

ZeroMQChannel::~ZeroMQChannel() VW_NOTHROW {
  if (m_sock) {
    VW_ASSERT(m_id == Thread::id(), LogicErr() << "ZeroMQChannels must not be shared between threads");
    m_sock.reset();
  }
}

void ZeroMQChannel::send_bytes(const uint8* message, size_t len) {
  VW_ASSERT(m_id == Thread::id(), LogicErr() << "ZeroMQChannels must not be shared between threads");
  // send() just queues the message. need to copy it.
  zmq::message_t rmsg(len);
  ::memcpy(rmsg.data(), message, len);

  // false return means EAGAIN
  while (!m_sock->send(rmsg)) {
    Thread::sleep_ms(10);
  }
}

bool ZeroMQChannel::recv_bytes(std::vector<uint8>* bytes) {
  VW_ASSERT(m_id == Thread::id(), LogicErr() << "ZeroMQChannels must not be shared between threads");

  zmq::message_t rmsg;

  if (this->timeout() > -1) {
    zmq::pollitem_t item[1];
    item[0].socket = *m_sock;
    item[0].events = ZMQ_POLLIN;
    item[0].revents = 0;
    item[0].fd = -1;

    // zeromq timeouts are in us
    int32 timeout_us = this->timeout() * 1000;
    const uint64_t stop = Stopwatch::microtime() + timeout_us;
    int rc;
    while (true) {
      rc = ::zmq_poll(item, 1, timeout_us);
      if (rc > 0)
        break;
      if (rc < 0 && !(errno == EINTR || errno == EAGAIN))
        vw_throw(ZeroMQErr() << "ZeroMQChannel::recv_bytes failed: poll error: " << zmq_err());

      //recalculate time and go again. poll can return 0 before timeout has
      //completely expired, though zmq revision
      //e2167cecaefec6557c7a5712fb75e51487ff69a6 might fix that. anyway, assume
      //it might have returned early and deal with it.
      timeout_us = stop - Stopwatch::microtime();
      if (timeout_us <= 0)
        return false;
    }
  }

  // false return means EAGAIN
  while (!m_sock->recv(&rmsg))
    Thread::sleep_ms(10);

  bytes->clear();
  bytes->reserve(rmsg.size());
  std::copy(reinterpret_cast<const char*>(rmsg.data()),
            reinterpret_cast<const char*>(rmsg.data()) +rmsg.size(),
            std::back_inserter(*bytes));

  return true;
}

void ZeroMQChannel::CallMethod(const pb::MethodDescriptor* method,
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

    send_message(q_wrap);
    switch (recv_message(a_wrap)) {
      case 0:
        // we can't recover from a timeout, since we can't send twice in a row
        vw_throw(NetworkErr() << "CallMethod Timeout.");
      case -1:
        vw_out(WarningMessage) << "CallMethod(): corrupted message. ";
        continue;
      default:
        throw_rpc_error(a_wrap.error());
        response->ParseFromString(a_wrap.payload());
        return;
    }
  }
  vw_out(WarningMessage) << "No more retries." << std::endl;
  vw_throw(RpcErr() << "CallMethod timed out completely");
}

// Updates the url in place
void ZeroMQChannel::init_endpoint(Url& endpoint) {
  VW_ASSERT(!m_sock, LogicErr() << "Please don't reuse ZeroMQ channels");

  if (endpoint.scheme() == "zmq")
    endpoint.scheme("zmq+tcp");

  VW_ASSERT(endpoint.scheme().substr(0,4) == "zmq+", LogicErr() << "Expected a zmq+ url");

  if (endpoint.scheme() == "zmq+tcp") {
    VW_ASSERT(endpoint.path() == "/", ArgumentErr() << "zmq+tcp urls do not support paths");
    VW_ASSERT(endpoint.port() != 0,   ArgumentErr() << "zmq+tcp urls must have a port");
  }
  if (endpoint.scheme() == "zmq+ipc") {
    VW_ASSERT(endpoint.path() != "/", ArgumentErr() << "zmq+ipc urls must have a path");
  }

  // remove zmq+ prefix to the scheme
  endpoint.scheme(endpoint.scheme().substr(4));

  endpoint.query().clear();
  endpoint.fragment("");

  m_id = Thread::id();
}

void ZeroMQChannel::conn(const Url& endpoint_) {
  Url endpoint(endpoint_);
  this->init_endpoint(endpoint);

  // We don't really know what the underlying resolver can handle, but a single
  // wildcard is definitely wrong for a client url.
  if (endpoint.scheme() == "tcp") {
    VW_ASSERT(endpoint.hostname() != "*", ArgumentErr() << "zmq+tcp clients cannot connect to *");
  }

  // remove trailing / on the path, if it's there
  std::string url = endpoint.string();
  boost::algorithm::trim_right_if(url, boost::is_any_of("/"));

  m_sock.reset(new zmq::socket_t(*m_ctx, ZMQ_REQ));

  // We use the c interface rather than the c++ one here to avoid having to
  // catch an exception and rethrow
  if (zmq_connect(*m_sock, url.c_str()) != 0)
    vw_throw(ZeroMQErr() << "Failed to connect to " << url.c_str() << ": " << zmq_err);
}

void ZeroMQChannel::bind(const Url& endpoint_) {
  Url endpoint(endpoint_);
  this->init_endpoint(endpoint);

  // remove trailing / on the path, if it's there
  std::string url = endpoint.string();
  boost::algorithm::trim_right_if(url, boost::is_any_of("/"));

  m_sock.reset(new zmq::socket_t(*m_ctx, ZMQ_REP));

  // We use the c interface rather than the c++ one here to avoid having to
  // catch an exception and rethrow
  if (zmq_bind(*m_sock, url.c_str()) != 0)
    vw_throw(ZeroMQErr() << "Failed to bind to " << url.c_str() << ": " << zmq_strerror(zmq_errno()));
}

int32 ZeroMQChannel::timeout() const {
  return m_timeout;
}

void ZeroMQChannel::set_timeout(int32 x) {
  m_timeout = x;
}

uint32 ZeroMQChannel::retries() const {
  return m_retries;
}

void ZeroMQChannel::set_retries(uint32 x) {
  m_retries = x;
}

std::string ZeroMQChannel::name() const {
  return m_human_name;
}
