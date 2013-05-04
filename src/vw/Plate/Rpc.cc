// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
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
    bool m_error;
    std::string m_error_msg;
    vw::Mutex     m_startup_mutex;
    vw::Condition m_startup_cond;
  public:
    Task(RpcBase* rpc, const Url& u, ThreadMap& stats)
      : m_rpc(rpc), m_url(u), m_stats(stats), m_go(true), m_error(false), m_error_msg("") {}

    void operator()();
    void stop() {m_go = false;}
    bool error() const { return m_error; }
    Mutex&     mutex() {return m_startup_mutex;}
    Condition& cond()  {return m_startup_cond;}
    const std::string& error_msg() const { return m_error_msg; }
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
  Mutex::Lock lock(m_task->mutex());
  m_thread.reset(new Thread(m_task));
  m_task->cond().wait(m_task->mutex());
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

const char* RpcServerBase::error() const {
  if (!m_task)
    return "Server message task is gone!";
  if (!m_task->error())
    return 0;
  return m_task->error_msg().c_str();
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
  Mutex::Lock lock(mutex());
  try {
    // TODO: pass something more useful than u.string()
    m_chan.reset(IChannel::make_bind(m_url, m_url.string()));
    m_chan->set_timeout(250);
  } catch (const std::exception& e) {
    m_chan.reset();
    m_error = true;
    m_error_msg = e.what();
  }

  cond().notify_all();
  lock.unlock();

  try {
    if (!m_error) {
      while (m_go) {
        // we could try to recreate the channel here, but it's better if we die
        // and let the clients die, too
        handle_one_request();
      }
      m_chan.reset();
    }
  } catch (const std::exception& e) {
    m_chan.reset();
    m_error_msg = e.what();
    m_error = true;
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
    VW_OUT(WarningMessage) << "Server Error: " << e.what() << std::endl;
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
  VW_ASSERT(m_chan, LogicErr() << "Cannot set timeout before constructing channel");
  m_chan->set_timeout(t);
}

void RpcClientBase::set_retries(uint32 t) {
  VW_ASSERT(m_chan, LogicErr() << "Cannot set retries before constructing channel");
  m_chan->set_retries(t);
}
