#include <vw/Plate/RpcChannel.h>
#include <vw/Plate/Rpc.h>
#include <vw/Core/Features.h>

#if defined(VW_HAVE_PKG_RABBITMQ_C) && VW_HAVE_PKG_RABBITMQ_C==1
#include <vw/Plate/detail/AmqpChannel.h>
#endif

#if defined(VW_HAVE_PKG_ZEROMQ) && VW_HAVE_PKG_ZEROMQ==1
#include <vw/Plate/detail/ZeroMQChannel.h>
#endif

#include <vw/Plate/HTTPUtils.h>
#include <vw/Core/Stopwatch.h>
#include <google/protobuf/descriptor.h>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/scoped_array.hpp>

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

void IChannel::send_message(RpcWrapper& message) {
  message.set_checksum(this->checksum(message));
  const size_t len = message.ByteSize();

  boost::scoped_array<uint8> bytes(new uint8[len]);
  message.SerializeToArray(bytes.get(), boost::numeric_cast<int>(len));
  send_bytes(bytes.get(), len);
}

int32 IChannel::recv_message(RpcWrapper& message) {
  std::vector<uint8> bytes;
  if (!recv_bytes(&bytes))
    return 0;
  if (!message.ParseFromArray(&bytes[0], boost::numeric_cast<int>(bytes.size())))
    return -1;
  if (this->checksum(message) != message.checksum())
    return -1;
  return 1;
}

IChannel* IChannel::make(const std::string& scheme, const std::string& clientname) {

#if defined(VW_HAVE_PKG_RABBITMQ_C) && VW_HAVE_PKG_RABBITMQ_C==1
  if (scheme == "pf" || scheme == "amqp")
    return new AmqpChannel(clientname);
#endif

#if defined(VW_HAVE_PKG_ZEROMQ) && VW_HAVE_PKG_ZEROMQ==1
  if (scheme == "zmq" || scheme == "zmq+ipc" || scheme == "zmq+tcp" || scheme == "zmq+inproc")
    return new ZeroMQChannel(clientname);
#endif

  vw_throw(NoImplErr() << "Unsupported channel scheme: " << scheme);
}

IChannel* IChannel::make_conn(const Url& u, const std::string& clientname) {
  std::auto_ptr<IChannel> chan(make(u.scheme(), clientname));

  if (u.query().has("timeout"))
    chan->set_timeout(u.query().get<int32>("timeout"));
  if (u.query().has("retries"))
    chan->set_retries(u.query().get<uint32>("retries"));

  chan->conn(u);
  return chan.release();
}

IChannel* IChannel::make_bind(const Url& u, const std::string& clientname) {
  std::auto_ptr<IChannel> chan(make(u.scheme(), clientname));
  if (u.query().has("timeout"))
    chan->set_timeout(u.query().get<int32>("timeout"));
  if (u.query().has("retries"))
    chan->set_retries(u.query().get<uint32>("retries"));

  chan->bind(u);
  return chan.release();
}

namespace {
  // This is basically fletcher16, except not quite (since it preserves the
  // full c0 and c1)
  void fletcher_add( uint32& sum, const char *data, size_t len ) {
    uint16 c0 = (sum >> 16), c1 = (sum & 0xffff);
    for (size_t i = 0; i < len; ++i) {
      c0 += *data++;
      c1 += c0;
    }
    sum = c0;
    sum <<= 16;
    sum += c1;
  }
}

uint32 IChannel::checksum(const RpcWrapper& message) {
  uint32 sum = 0;
  fletcher_add(sum, message.method().c_str(), message.method().size());
  fletcher_add(sum, message.payload().c_str(), message.payload().size());
  return sum;
}
