// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/Exception.h>
#include <vw/Plate/Rpc.h>
#include <vw/Plate/RpcChannel.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Plate/Rpc.pb.h>
#include <vw/Core/Log.h>
#include <vw/Core/Debugging.h>
#include <google/protobuf/descriptor.h>
#include <boost/scoped_ptr.hpp>

using namespace vw;
using namespace vw::platefile;
namespace pb = ::google::protobuf;

#define NOIMPL { vw_throw(vw::NoImplErr() << "Not implemented: " << VW_CURRENT_FUNCTION); }

void RpcBase::Reset() NOIMPL
bool RpcBase::Failed() const NOIMPL
std::string RpcBase::ErrorText() const NOIMPL
void RpcBase::StartCancel() NOIMPL
void RpcBase::SetFailed(const std::string& /*reason*/) NOIMPL
bool RpcBase::IsCanceled() const NOIMPL
void RpcBase::NotifyOnCancel(::google::protobuf::Closure* /*callback*/) NOIMPL

class RpcServerBase::Task {
    RpcBase* m_rpc;
    const Url m_url;
    boost::shared_ptr<IChannel> m_chan;
    ThreadMap& m_stats;
    bool m_go;
  public:
    Task(RpcBase* rpc, const Url& u, ThreadMap& stats)
      : m_rpc(rpc), m_url(u), m_stats(stats), m_go(true) {}

    void operator()();
    void stop() {m_go = false;}
  protected:
    // return = false means timeout
    bool handle_one_request();
};

ThreadMap::Locked::Locked(ThreadMap& m)
  : m_map(m), m_lock(new Mutex::Lock(m_map.m_mutex)) {}

ThreadMap::Locked::~Locked() {
  m_lock.reset();
}

void ThreadMap::Locked::add(const std::string& key, vw::int64 val) {
  m_map.m_data[key] += val;
}

int64 ThreadMap::Locked::get(const std::string& key) const {
  return m_map.m_data[key];
}

void ThreadMap::Locked::clear() {
  m_map.m_data.clear();
}

void ThreadMap::add(const std::string& key, vw::int64 val) {
  Mutex::Lock lock(m_mutex);
  m_data[key] += val;
}

vw::int64 ThreadMap::get(const std::string& key) const {
  Mutex::Lock lock(m_mutex);
  map_t::const_iterator i = m_data.find(key);
  if (i == m_data.end())
    return 0;
  return i->second;
}

void ThreadMap::clear() {
  Mutex::Lock lock(m_mutex);
  m_data.clear();
}

void RpcServerBase::launch_thread(const Url& url) {
  m_task.reset(new Task(this, url, m_stats));
  m_thread.reset(new Thread(m_task));
}

RpcServerBase::RpcServerBase(const Url& url) {
  launch_thread(url);
}

RpcServerBase::~RpcServerBase() {
  stop();
  m_thread.reset();
  m_task.reset();
}

void RpcServerBase::stop() {
  if (m_task)
    m_task->stop();
  if (m_thread)
    m_thread->join();
}

void RpcServerBase::bind(const Url& url) {
  launch_thread(url);
}

ThreadMap::Locked RpcServerBase::stats() {
  return ThreadMap::Locked(m_stats);
}

void RpcBase::set_debug(bool on) {
  m_debug = on;
}

bool RpcBase::debug() const {
  return m_debug;
}

void RpcServerBase::Task::operator()() {
  try {
    // TODO: pass something more useful than u.string()
    m_chan.reset(IChannel::make_bind(m_url, m_url.string()));
    m_chan->set_timeout(250);

    while (m_go) {
      try {
        handle_one_request();
      } catch (const NetworkErr& err) {
        vw_out(ErrorMessage) << "Network error! This is probably fatal. Recreating channel." << std::endl;
        m_stats.add("fatal_error");
        // TODO: pass something more useful than u.string()
        m_chan.reset(IChannel::make_bind(m_url, m_url.string()));
        m_chan->set_timeout(250);
      }
    }
    m_chan.reset();
  } catch (const std::exception& e) {
    std::cerr << VW_CURRENT_FUNCTION << ": caught exception: " << e.what() << std::endl;
    std::abort();
  }
}

bool RpcServerBase::Task::handle_one_request() {
  RpcWrapper q_wrap, a_wrap;

  switch (m_chan->recv_message(q_wrap)) {
    case 0:
      m_stats.add("timeout");
      return false;
    case -1:
      // message corruption?
      m_stats.add("client_error");
      return false;
    default:
      break;
  }

  const pb::MethodDescriptor* method = m_rpc->service()->GetDescriptor()->FindMethodByName(q_wrap.method());
  VW_ASSERT(method, RpcErr() << "Unrecognized RPC method: " << q_wrap.method());

  typedef boost::scoped_ptr<pb::Message> msg_t;

  msg_t q(m_rpc->service()->GetRequestPrototype(method).New());
  msg_t a(m_rpc->service()->GetResponsePrototype(method).New());

  a_wrap.set_method(q_wrap.method());

  // Attempt to parse the actual request message from the request_wrapper.
  if (!q->ParseFromString(q_wrap.payload())) {
    // they're not speaking our protocol. just ignore it.
    m_stats.add("client_error");
    return false;
  }

  // copy the metadata over
  a_wrap.set_requestor(q_wrap.requestor());
  if (q_wrap.seq() != a_wrap.seq())
    a_wrap.set_seq(q_wrap.seq());

  try {
    m_rpc->service()->CallMethod(method, m_rpc, q.get(), a.get(), null_callback());
    a_wrap.set_payload(a->SerializeAsString());
    a_wrap.mutable_error()->set_code(RpcErrorMsg::SUCCESS);
    m_stats.add("msgs");
  } catch (const NetworkErr &e) {
    // we can't respond to these, unfortunately... channel might be down. just
    // rethrow, let higher-level logic work it out.
    throw;
  } catch (const PlatefileErr &e) {
    // These exceptions should be reported back to the far side as a specific
    // error code
    a_wrap.mutable_error()->set_code(e.code());
    a_wrap.mutable_error()->set_msg(e.what());
    m_stats.add("client_error");
  } catch (const std::exception &e) {
    // These exceptions should be reported back to the far side as a general
    // server error- unless we're in debug mode, in which case... rethrow!
    a_wrap.mutable_error()->set_code(RpcErrorMsg::REMOTE_ERROR);
    a_wrap.mutable_error()->set_msg(e.what());
    m_stats.add("server_error");
    vw_out(WarningMessage) << "Server Error: " << e.what() << std::endl;
    if (m_rpc->debug()) {
      m_chan->send_message(a_wrap);
      throw;
    }
  }

  m_chan->send_message(a_wrap);
  return true;
}

// TODO: Make the clientname settable here.
RpcClientBase::RpcClientBase(const Url& u)
  : m_chan(IChannel::make_conn(u, u.string())) {}

// TODO: Make the clientname settable here.
RpcClientBase::RpcClientBase(const Url& u, int32 timeout, uint32 retries)
  : m_chan(IChannel::make_conn(u, u.string())) {
  m_chan->set_timeout(timeout);
  m_chan->set_retries(retries);
}

pb::RpcChannel* RpcClientBase::base_channel() {
  return m_chan.get();
}

void RpcClientBase::set_timeout(int32 t) {
  m_chan->set_timeout(t);
}

void RpcClientBase::set_retries(uint32 t) {
  m_chan->set_retries(t);
}
