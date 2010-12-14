// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>
#include <vw/Plate/Rpc.h>
#include <vw/Plate/RpcChannel.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Plate/Exception.h>
#include <vw/Plate/tests/TestRpcService.pb.h>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <cstdlib>

using namespace std;
using namespace vw;
using namespace vw::platefile;
using namespace vw::test;

namespace pb = ::google::protobuf;

#if defined(VW_HAVE_PKG_ZEROMQ) && (VW_HAVE_PKG_ZEROMQ==1)
#define HAS_ZEROMQ(x) x
#else
#define HAS_ZEROMQ(x) DISABLED_ ## x
#endif

class TestServiceImpl : public TestService {
  public:
    static const int CLIENT_ERROR = 900001;
    static const int SERVER_ERROR = 900002;

    virtual void DoubleRequest(pb::RpcController* /*controller*/, const DoubleMessage* request, DoubleMessage* response, pb::Closure* done) {
      switch (request->num()) {
        case CLIENT_ERROR:
          done->Run();
          vw_throw(PlatefileErr() << "Sad panda");
        case SERVER_ERROR:
          done->Run();
          vw_throw(LogicErr() << "It broke! [<--- this is expected! Please ignore.]");
        default:
          response->set_num(request->num() * 2);
          done->Run();
      }
    }
};

#if 0
class TestServiceRude : public TestServiceImpl {
    const uint32 m_skip;
    uint32 m_seq;
  public:
    TestServiceRude(uint32 skip) : m_skip(skip), m_seq(0) {}
    virtual void DoubleRequest(pb::RpcController* /*controller*/, const DoubleMessage* request, DoubleMessage* response, pb::Closure* done) {
      if (m_seq++ < m_skip)
        vw_throw(NetworkErr() <<
      response->set_num(request->num() * 2);
      done->Run();
    }
};
#endif

typedef RpcClient<TestService> TestClient;
typedef RpcServer<TestService> TestServer;
typedef boost::shared_ptr<TestClient> Client;
typedef boost::shared_ptr<TestServer> Server;

const uint64 TIMEOUT = 1000;

struct RpcTest : public ::testing::TestWithParam<Url> {
  Server server;
  vector<Client> clients;

  void SetUp() {
    make_server();
    make_clients(0);
  }

  void make_clients(uint64 count) {
    clients.resize(count);
    ASSERT_NO_THROW(std::generate(clients.begin(), clients.end(), TestClient::make_factory(GetParam(), TIMEOUT, 10)));
  }

  void make_server() {
    ASSERT_NO_THROW(server.reset(new RpcServer<TestService>(GetParam(), new TestServiceImpl())));
    // XXX: This is pretty lame. Unlike every other type of 0mq socket, inproc
    // sockets must bind() before connect(). This delay makes it likely the
    // server starts before the client.
    Thread::sleep_ms(50);
  }

  void TearDown() {
    clients.resize(0);
    server.reset();
  }
};

TEST_P(RpcTest, Basic) {
  make_clients(1);

  DoubleMessage q, a;
  for (uint32 i = 0; i < 1000; ++i) {
    q.set_num(i);
    EXPECT_NO_THROW(clients[0]->DoubleRequest(clients[0].get(), &q, &a, null_callback()));
    EXPECT_EQ(i*2, a.num());
  }
  EXPECT_EQ(1000, server->stats().get("msgs"));
  EXPECT_EQ(0, server->stats().get("client_error"));
}

TEST_P(RpcTest, Err) {
  make_clients(1);

  EXPECT_EQ(0, server->stats().get("msgs"));
  EXPECT_EQ(0, server->stats().get("client_error"));
  EXPECT_EQ(0, server->stats().get("server_error"));

  DoubleMessage q, a;

  q.set_num(TestServiceImpl::CLIENT_ERROR);
  EXPECT_THROW(clients[0]->DoubleRequest(clients[0].get(), &q, &a, null_callback()), PlatefileErr);
  EXPECT_EQ(0, server->stats().get("msgs"));
  EXPECT_EQ(1, server->stats().get("client_error"));
  EXPECT_EQ(0, server->stats().get("server_error"));

  q.set_num(TestServiceImpl::SERVER_ERROR);
  EXPECT_THROW(clients[0]->DoubleRequest(clients[0].get(), &q, &a, null_callback()), RpcErr);
  EXPECT_EQ(0, server->stats().get("msgs"));
  EXPECT_EQ(1, server->stats().get("client_error"));
  EXPECT_EQ(1, server->stats().get("server_error"));
}

TEST(TestRpc, HAS_ZEROMQ(KillServerDeathTest)) {
  // This test only works if the death test starts in a different process
  // (so we don't get a channel collision)
  Url u("zmq+inproc://unittest2");
  Server server;
  Client client;

  ASSERT_NO_THROW(server.reset(new RpcServer<TestService>(u, new TestServiceImpl())));
  Thread::sleep_ms(50);
  ASSERT_NO_THROW(client.reset(new TestClient(u)));

  DoubleMessage q, a;

  server->set_debug(true);
  q.set_num(TestServiceImpl::SERVER_ERROR);
  EXPECT_DEATH(
      client->DoubleRequest(client.get(), &q, &a, null_callback()),
      "broke");
}

class ClientTask {
    TestClient::Factory f;
    Mutex &mutex;
    Condition &cond;
  public:
    static const uint64 MSG_COUNT = 250;
    ClientTask(TestClient::Factory f, Mutex& m, Condition& c) : f(f), mutex(m), cond(c) {}

    void operator()() {
      uint32 offset = Thread::id() * MSG_COUNT;
      Client c = f();
      DoubleMessage q, a;

      {
        Mutex::Lock lock(mutex);
        cond.timed_wait(lock, 1000);
      }

      for (uint32 i = 0; i < MSG_COUNT; ++i) {
        q.set_num(i + offset);
        EXPECT_NO_THROW(c->DoubleRequest(c.get(), &q, &a, null_callback()));
        EXPECT_EQ(q.num() * 2, a.num());
      }
      c.reset();
    }
};

TEST_P(RpcTest, MultiThreadTorture) {
  uint64 THREAD_COUNT = 20;

  typedef boost::shared_ptr<ClientTask> task_t;
  typedef boost::shared_ptr<Thread> thread_t;

  vector<task_t>   tasks(THREAD_COUNT);
  vector<thread_t> threads(THREAD_COUNT);

  // no retries
  TestClient::Factory f = TestClient::make_factory(GetParam(), TIMEOUT, 0);

  EXPECT_EQ(0, server->stats().get("msgs"));
  EXPECT_EQ(0, server->stats().get("client_error"));
  EXPECT_EQ(0, server->stats().get("server_error"));

  Mutex m;
  Condition c;

  for (uint64 i = 0; i < THREAD_COUNT; ++i) {
    tasks[i]   = task_t(new ClientTask(f, m, c));
    threads[i] = thread_t(new Thread(tasks[i]));
  }

  Thread::sleep_ms(100);
  // start everyone at the same time
  c.notify_all();

  BOOST_FOREACH(thread_t& t, threads)
    t->join();

  EXPECT_EQ(THREAD_COUNT * ClientTask::MSG_COUNT, server->stats().get("msgs"));
  EXPECT_EQ(0, server->stats().get("client_error"));
  EXPECT_EQ(0, server->stats().get("server_error"));
}

Url amqp_url(string hostname = "", short port = -1) {
  if (hostname.empty())
    hostname = getenv2("AMQP_TEST_HOSTNAME", "localhost");

  string sport = getenv2("AMQP_TEST_PORT", "5672");
  if (port != -1) {
    sport = boost::lexical_cast<string>(port);
  }
  return string("amqp://") + hostname + ":" + sport + "/unittest/server";
}

std::vector<Url> test_urls() {
  std::vector<Url> v;
#if defined(VW_HAVE_PKG_RABBITMQ_C) && VW_HAVE_PKG_RABBITMQ_C==1
  v.push_back(amqp_url());
#endif
#if defined(VW_HAVE_PKG_ZEROMQ) && VW_HAVE_PKG_ZEROMQ==1
    v.push_back(Url("zmq+ipc://" TEST_OBJDIR "/unittest"));
    v.push_back(Url("zmq+tcp://127.0.0.1:54321"));
    v.push_back(Url("zmq+inproc://unittest"));
#endif
    return v;
}

INSTANTIATE_TEST_CASE_P(URLs, RpcTest, ::testing::ValuesIn(test_urls()));
