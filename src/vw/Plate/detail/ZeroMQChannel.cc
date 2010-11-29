// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/detail/ZeroMQChannel.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Core/Log.h>

#include <boost/algorithm/string/trim.hpp>
#include <google/protobuf/descriptor.h>
#include <zmq.hpp>

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
}

ZeroMQChannel::ZeroMQChannel(const std::string& human_name)
  : m_ctx(get_ctx()), m_human_name(human_name), m_timeout(DEFAULT_TIMEOUT), m_retries(10) {}


ZeroMQChannel::~ZeroMQChannel() VW_NOTHROW {
  if (m_sock) {
    VW_ASSERT(m_id == Thread::id(), LogicErr() << "ZeroMQChannels must not be shared between threads");
    m_sock.reset();
  }
}

void ZeroMQChannel::send_bytes(const uint8* message, size_t len) {
  VW_ASSERT(m_id == Thread::id(), LogicErr() << "ZeroMQChannels must not be shared between threads");
  zmq::message_t rmsg(len);
  ::memcpy(rmsg.data(), message, len);
  // send() just queues the message. need to copy it.
  if (!m_sock->send(rmsg))
    vw_throw(ZeroMQErr() << "Failed to send bytes");
}

bool ZeroMQChannel::recv_bytes(SharedByteArray& bytes) {
  VW_ASSERT(m_id == Thread::id(), LogicErr() << "ZeroMQChannels must not be shared between threads");

  zmq::message_t rmsg;

  zmq::pollitem_t item[1];
  item[0].socket = *m_sock.get();
  item[0].events = ZMQ_POLLIN;
  item[0].revents = 0;
  item[0].fd = -1;

  if (this->timeout() == -1) {
    if (!m_sock->recv(&rmsg))
      vw_throw(ZeroMQErr() << "Failed to recv bytes");
  } else {
    // zeromq timeouts are in us
    int32 timeout_us = this->timeout() * 1000;

    uint64_t stop = Stopwatch::microtime() + timeout_us;

    // This is only necessary for a while
    // zmq revision e2167cecaefec6557c7a5712fb75e51487ff69a6 prevents poll from
    // returning early (so maybe v2.1.0?)
    int rc;
    do {
      rc = zmq::poll(item, 1, timeout_us);
    } while (rc == 0 && Stopwatch::microtime() < stop);

    if (rc == 0)
      return false;

    if (!m_sock->recv(&rmsg))
      vw_throw(ZeroMQErr() << "Failed to recv bytes");
  }

  bytes.reset(new ByteArray(reinterpret_cast<const char*>(rmsg.data()),
                            reinterpret_cast<const char*>(rmsg.data()) + rmsg.size()));
  return true;
}

void ZeroMQChannel::CallMethod(const pb::MethodDescriptor* method,
                               pb::RpcController* /*controller*/,
                               const pb::Message* request,
                               pb::Message* response,
                               pb::Closure* done)
{
  RpcWrapper q_wrap, a_wrap;

  q_wrap.set_method(method->name());
  q_wrap.set_payload(request->SerializeAsString());
  q_wrap.set_requestor(this->name());

  for (uint32 trial = 0; trial <= m_retries; ++trial) {
    if (trial > 0)
      vw_out(WarningMessage) << "Retry (" << trial << "/" << m_retries << ")" << std::endl;

    send_message(q_wrap);
    if (!recv_message(a_wrap)) {
      vw_out(WarningMessage) << "CallMethod Timeout. ";
      continue;
    }

    throw_rpc_error(a_wrap.error());

    response->ParseFromString(a_wrap.payload());
    done->Run();
    return;
  }
  vw_out(WarningMessage) << "No more retries." << std::endl;
  vw_throw(RpcErr() << "CallMethod timed out completely");
}

void ZeroMQChannel::conn(const Url& endpoint_) {
  Url endpoint(endpoint_);
  if (endpoint.scheme() == "zmq")
    endpoint.scheme("zmq+tcp");

  VW_ASSERT(endpoint.scheme().substr(0,4) == "zmq+", LogicErr() << "Expected a zmq+ url");
  VW_ASSERT(!m_sock, LogicErr() << "Please don't reuse ZeroMQ channels");

  std::string url = endpoint.string();
  // remove zmq+ prefix
  url = url.substr(4);

  // remove trailing /, if it's there
  boost::algorithm::trim_right_if(url, boost::is_any_of("/"));

  m_id = Thread::id();
  m_sock.reset(new zmq::socket_t(*m_ctx, ZMQ_REQ));
  m_sock->connect(url.c_str());
}

void ZeroMQChannel::bind(const Url& endpoint_) {
  Url endpoint(endpoint_);
  if (endpoint.scheme() == "zmq")
    endpoint.scheme("zmq+tcp");

  VW_ASSERT(endpoint.scheme().substr(0,4) == "zmq+", LogicErr() << "Expected a zmq+ url");
  VW_ASSERT(!m_sock, LogicErr() << "Please don't reuse ZeroMQ channels");

  std::string url = endpoint.string();
  // remove zmq+ prefix
  url = url.substr(4);

  // remove trailing /, if it's there
  boost::algorithm::trim_right_if(url, boost::is_any_of("/"));

  m_id = Thread::id();
  m_sock.reset(new zmq::socket_t(*m_ctx, ZMQ_REP));
  m_sock->bind(url.c_str());
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

#if 0
uint64 ZeroMQChannel::queue_depth() const {
  vw_throw(NoImplErr() << "Dang. Can't do this with zeromq yet.");
}
#endif

std::string ZeroMQChannel::name() const {
  return m_human_name;
}
