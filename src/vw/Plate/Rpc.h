// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_RPC_H__
#define __VW_PLATE_RPC_H__

#include <vw/Core/Thread.h>
#include <google/protobuf/service.h>
#include <boost/lambda/construct.hpp>
#include <boost/lambda/bind.hpp>

namespace vw {
namespace platefile {

  class IChannel;
  class Url;

  inline google::protobuf::Closure* null_callback() {
    return google::protobuf::NewCallback(&google::protobuf::DoNothing);
  }

  class ThreadMap : private boost::noncopyable {
    private:
      typedef std::map<std::string, int64> map_t;
      map_t m_data;
      mutable Mutex m_mutex;
      friend class Locked;
    public:
      class Locked {
        private:
          ThreadMap& m_map;
          boost::shared_ptr<Mutex::Lock> m_lock;
        public:
          Locked(ThreadMap& m);
          ~Locked();
          void add(const std::string& key, int64 val = 1);
          int64 get(const std::string& key) const;
          void clear();
      };

      void add(const std::string& key, int64 val = 1);
      int64 get(const std::string& key) const;
      void clear();
  };

  class RpcBase : public ::google::protobuf::RpcController, private boost::noncopyable {
    private:
      // Client-side methods ---------------------------------------------
      void Reset();
      bool Failed() const;
      std::string ErrorText() const;
      void StartCancel();

      // Server-side methods ---------------------------------------------
      void SetFailed(const std::string& /*reason*/);
      bool IsCanceled() const;
      void NotifyOnCancel(::google::protobuf::Closure* /*callback*/);

    protected:
      bool m_debug;

    public:
      // retrieve the abstract service
      virtual ::google::protobuf::Service* service() = 0;
      RpcBase() : m_debug(false) {}
      virtual ~RpcBase() {}

      void set_debug(bool on);
      bool debug() const;
  };

  class RpcServerBase : public RpcBase {
    private:
      class Task;
      boost::shared_ptr<Task> m_task;
      boost::shared_ptr<Thread> m_thread;
      ThreadMap m_stats;

      void launch_thread(const Url& u);

    public:
      RpcServerBase() {}
      RpcServerBase(const Url& url);
      virtual ~RpcServerBase();

      void bind(const Url& url);
      void stop();

      ThreadMap::Locked stats();
  };

  class RpcClientBase : public RpcBase {
    private:
      boost::shared_ptr<IChannel> m_chan;
    protected:
      // This exists to allow a fwd-decl of IChannel
      ::google::protobuf::RpcChannel* base_channel();
    public:
      RpcClientBase(const Url& u);
      RpcClientBase(const Url& u, int32 timeout, uint32 retries);
      virtual ~RpcClientBase() {}

      // -1 means "never", other values in ms
      void set_timeout(int32 t);
      void set_retries(uint32 t);
  };

  template <typename ServiceT>
  class RpcServer : public RpcServerBase {
    private:
      boost::shared_ptr<ServiceT> m_service;
    public:
      RpcServer() {}
      // Shares ownership of service
      RpcServer(const Url& url, boost::shared_ptr<ServiceT> service)
        : RpcServerBase(url), m_service(service) {}
      // Takes ownership of service
      RpcServer(const Url& url, ServiceT* service)
        : RpcServerBase(url), m_service(service) {}

      void bind(const Url& url, boost::shared_ptr<ServiceT> service);

      ::google::protobuf::Service* service() {return m_service.get();}
      ServiceT* impl() {return m_service.get();}
  };

  template <typename ServiceT>
  class RpcClient : public RpcClientBase, public ServiceT::Stub {
    public:
      ::google::protobuf::Service* service() {return this;}

      RpcClient(const Url& url)
        : RpcClientBase(url), ServiceT::Stub(base_channel()) {}

      RpcClient(const Url& url, int32 timeout, uint32 retries)
        : RpcClientBase(url, timeout, retries), ServiceT::Stub(base_channel()) {}

      typedef boost::function<boost::shared_ptr<RpcClient> ()> Factory;

      static Factory make_factory(const Url& url) {
        using namespace boost;
        typedef boost::shared_ptr<RpcClient> ptr_t;
        return lambda::bind(lambda::constructor<ptr_t>(), lambda::bind(lambda::new_ptr<RpcClient>(), url));
      }
      static Factory make_factory(const Url& url, int32 timeout, uint32 retries) {
        using namespace boost;
        typedef boost::shared_ptr<RpcClient> ptr_t;
        return lambda::bind(lambda::constructor<ptr_t>(), lambda::bind(lambda::new_ptr<RpcClient>(), url, timeout, retries));
      }
  };

}} // namespace vw::platefile

#endif
