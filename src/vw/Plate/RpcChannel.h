// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_RPCCHANNEL_H__
#define __VW_PLATE_RPCCHANNEL_H__

#include <vw/Core/FundamentalTypes.h>
#include <google/protobuf/service.h>
#include <vector>

namespace vw {
namespace platefile {

class RpcWrapper;
class Url;
class RpcBase;

// Returns a unique name based on the identifier
// (uuid built from id, hostname, pid, thread, time)
std::string unique_name(const std::string& identifier, const std::string& join = "_");

// IChannel has thread affinity, and must not be called from any thread other
// than the constructing one
class IChannel : public ::google::protobuf::RpcChannel {
  public:
    virtual void send_bytes(const uint8* message, size_t len) = 0;
    virtual bool recv_bytes(std::vector<uint8>* bytes) VW_WARN_UNUSED = 0;

    void send_message(RpcWrapper& message);
    // This returns -1 for recoverable error, 0 for timeout, 1 for success
    // It throws for nonrecoverable errors (that's sort of ugly).
    int32 recv_message(RpcWrapper& message) VW_WARN_UNUSED;

    // Calculate checksum.
    static uint32 checksum(const RpcWrapper& message);

    virtual void CallMethod(const google::protobuf::MethodDescriptor*,
                            google::protobuf::RpcController*,
                            const google::protobuf::Message*,
                            google::protobuf::Message*,
                            google::protobuf::Closure*) = 0;

    // -1 means "never timeout", other values in ms
    virtual void set_timeout(int32 val) = 0;
    virtual void set_retries(uint32 r) = 0;

    virtual int32  timeout() const = 0;
    virtual uint32 retries() const = 0;

    // Defaults for all channels, if none is provided
    static const int32  DEFAULT_TIMEOUT = 10000; // ms
    static const uint32 DEFAULT_RETRIES = 10;

    //virtual uint64 queue_depth() const = 0;
    virtual std::string name() const = 0;

    IChannel() {}
    virtual ~IChannel() {}

    // Factories
    static IChannel* make(const std::string& scheme, const std::string& clientname);
    static IChannel* make_conn(const Url& url, const std::string& clientname);
    static IChannel* make_bind(const Url& url, const std::string& clientname);

    // Form connections
    virtual void conn(const Url& server) = 0;
    virtual void bind(const Url& self) = 0;
};


}} //  vw::platefile

#endif
