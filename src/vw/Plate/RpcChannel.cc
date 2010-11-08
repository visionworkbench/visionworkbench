#include <vw/Plate/RpcChannel.h>
#include <vw/Plate/Rpc.h>
#include <vw/Core/Features.h>

#if defined(VW_HAVE_PKG_RABBITMQ_C) && VW_HAVE_PKG_RABBITMQ_C==1
#include <vw/Plate/AmqpChannel.h>
#endif

#if defined(VW_HAVE_PKG_ZEROMQ) && VW_HAVE_PKG_ZEROMQ==1
#include <vw/Plate/ZeroMQChannel.h>
#endif

#include <vw/Plate/HTTPUtils.h>
#include <vw/Core/Stopwatch.h>
#include <google/protobuf/descriptor.h>
#include <boost/numeric/conversion/cast.hpp>

using namespace vw;
using namespace vw::platefile;
namespace pb = ::google::protobuf;

std::string vw::platefile::unique_name(const std::string& identifier, const std::string& join) {
  // Start by generating a unique name based on our hostname, PID, and thread ID.
  char hostname[255];
  gethostname(hostname, 255);
  std::ostringstream uuid;
  uuid << identifier
       << join << hostname
       << join << getpid()
       << join << Thread::id()
       << join << vw::Stopwatch::microtime(false);
  return uuid.str();
}

void IChannel::send_message(const RpcWrapper& message) {
  const size_t len = message.ByteSize();

  boost::scoped_array<uint8> bytes(new uint8[len]);
  message.SerializeToArray(bytes.get(), boost::numeric_cast<int>(len));
  send_bytes(bytes.get(), len);
}

bool IChannel::recv_message(RpcWrapper& message) {
  SharedByteArray bytes;
  if (!recv_bytes(bytes))
    return false;
  bool worked = message.ParseFromArray(bytes->begin(), boost::numeric_cast<int>(bytes->size()));
  VW_ASSERT(worked, RpcErr() << "IChannel::recv_message failed to parse message");
  return true;
}

IChannel* IChannel::make(const std::string& scheme, const std::string& clientname) {

#if defined(VW_HAVE_PKG_RABBITMQ_C) && VW_HAVE_PKG_RABBITMQ_C==1
  if (scheme == "pf" || scheme == "amqp")
    return new AmqpChannel(clientname);
#endif

#if defined(VW_HAVE_PKG_ZEROMQ) && VW_HAVE_PKG_ZEROMQ==1
  if (scheme == "zmq+ipc" || scheme == "zmq+tcp" || scheme == "zmq+inproc")
    return new ZeroMQChannel(clientname);
#endif

  vw_throw(NoImplErr() << "Unsupported channel scheme: " << scheme);
}

IChannel* IChannel::make_conn(const Url& u, const std::string& clientname) {
  std::auto_ptr<IChannel> chan(make(u.scheme(), clientname));
  chan->set_timeout(u.query().get("timeout", int32(-1)));
  chan->conn(u);
  return chan.release();
}

IChannel* IChannel::make_bind(const Url& u, const std::string& clientname) {
  std::auto_ptr<IChannel> chan(make(u.scheme(), clientname));
  chan->set_timeout(u.query().get("timeout", int32(-1)));
  chan->bind(u);
  return chan.release();
}
