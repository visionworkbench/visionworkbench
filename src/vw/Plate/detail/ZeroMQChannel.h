// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_ZEROMQCHANNEL_H__
#define __VW_PLATE_ZEROMQCHANNEL_H__

#include <vw/Plate/RpcChannel.h>
#include <vw/Plate/Exception.h>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

namespace zmq {
  class context_t;
  class socket_t;
}

namespace vw {
namespace platefile {

  // To recover from this, you may need to reconnect
  VW_DEFINE_EXCEPTION(ZeroMQErr, NetworkErr);

  class ZeroMQChannel: public IChannel,
                       private boost::noncopyable
  {
      boost::shared_ptr<zmq::context_t> m_ctx;
      boost::shared_ptr<zmq::socket_t>  m_sock;
      std::string m_human_name;
      uint64 m_id;
      int32 m_timeout;
      uint32 m_retries;

    protected:
      void send_bytes(const uint8* message, size_t len);
      bool recv_bytes(std::vector<uint8>* bytes);
      void init_endpoint(Url& u);

    public:
      ZeroMQChannel(const std::string& human_name);
      virtual ~ZeroMQChannel() VW_NOTHROW;

      int32  timeout() const;
      uint32 retries() const;

      void set_timeout(vw::int32);
      void set_retries(vw::uint32);

      //uint64 queue_depth() const;
      std::string name() const;

      virtual void CallMethod(const google::protobuf::MethodDescriptor*,
                              google::protobuf::RpcController*,
                              const google::protobuf::Message*,
                              google::protobuf::Message*,
                              google::protobuf::Closure*);

      // url format:
      // zmq://x -> maps directly to zmq urls
      void conn(const Url& server);
      void bind(const Url& self);
  };

}} // namespace vw::platefile

#endif
