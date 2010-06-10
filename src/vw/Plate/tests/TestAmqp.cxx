// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <vw/Plate/AmqpConnection.h>
#include <vw/Plate/ProtoBuffers.pb.h>
#include <vw/Plate/RpcServices.h>
#include <vw/Plate/IndexService.h>
#include <google/protobuf/stubs/common.h>

using namespace std;
using namespace vw;
using namespace vw::platefile;

#define EXCHANGE "unittest_exchange"
#define QUEUE "unittest_queue"
#define ROUTING_KEY "unittest"

const std::string default_hostname = "localhost";
const int default_port = 5672;


class AMQPTest : public ::testing::Test {
  protected:
    typedef boost::shared_ptr<AmqpConnection> Connection;
    typedef boost::shared_ptr<AmqpRpcServer> Server;
    typedef boost::shared_ptr<AmqpRpcClient> Client;
    typedef boost::shared_ptr<google::protobuf::Service> Service;

    Connection doconn(std::string hostname = "", int port = -1) {
      Connection c;
      try {
        if (hostname.empty()) {
          const char *val = getenv("AMQP_TEST_HOSTNAME");
          hostname = val ? std::string(val) : default_hostname;
        }
        if (port == -1) {
          const char *val = getenv("AMQP_TEST_PORT");
          port = val ? atoi(val) : default_port;
        }
        c.reset(new AmqpConnection(hostname, port));
      } catch (const AMQPErr& e) {
        // If we can't open the AMQP socket, treat this as a disabled test.
        // XXX: This is icky. We need a way to disable tests at runtime (or mock AMQP)
        std::cerr << "Could not open AMQP socket. Test disabled." << std::endl;
        c.reset();
      }
      return c;
    }

};

TEST_F(AMQPTest, BasicIO) {
  Connection conn = doconn();
  if (!conn)
    return;

  AmqpChannel chan(conn);

  chan.exchange_declare(EXCHANGE, "direct", true, false);
  chan.queue_declare(QUEUE, true, true, false);
  chan.queue_bind(QUEUE, EXCHANGE, ROUTING_KEY);

  IndexReadRequest request;
  request.set_platefile_id(23);

  request.set_col(1);
  request.set_row(2);
  request.set_level(3);
  request.set_transaction_id(32);
  request.set_exact_transaction_match(true);

  ByteArray msg;
  AmqpRpcEndpoint::serialize_message(request, msg);

  chan.basic_publish(msg, EXCHANGE, ROUTING_KEY);

  SharedByteArray response_bytes;
  if (!chan.basic_get(QUEUE, response_bytes))
    vw_throw(IOErr() << "Basic.Get failed");

  IndexReadRequest response;
  AmqpRpcEndpoint::parse_message(*response_bytes.get(), response);

  EXPECT_EQ(request.platefile_id(), response.platefile_id());
  EXPECT_EQ(request.col(), response.col());
  EXPECT_EQ(request.row(), response.row());
  EXPECT_EQ(request.level(), response.level());
}

TEST_F(AMQPTest, RPCSerialization) {
  Connection c1 = doconn(), c2 = doconn();
  if (!c1 || !c2)
    return;

  Server server( new AmqpRpcServer(c1, EXCHANGE, QUEUE"1") );
  Client client( new AmqpRpcClient(c2, EXCHANGE, QUEUE"2", "server") );

  Service srv_s( new IndexServiceImpl(TEST_SRCDIR) );
  Service srv_c( new IndexService::Stub(client.get() ) );

  server->bind_service(srv_s, "server");
  client->bind_service(srv_c, "client");

  IndexListRequest r1, r2;
  r1.set_type("toast");

  ASSERT_NO_THROW(client->send_message(r1, "server"));
  ASSERT_NO_THROW(server->get_message(r2));

  EXPECT_EQ(r1.type(), r2.type());
  EXPECT_EQ(r1.ByteSize(), r2.ByteSize());
}

TEST_F(AMQPTest, ChannelFailure) {
  Connection c1 = doconn(), c2 = doconn();
  if (!c1 || !c2)
    return;

  AmqpChannel ch1(c1), ch2(c2);
  EXPECT_NO_THROW(ch1.queue_declare(QUEUE, false, true, true));

  // Declare another queue on another connection with the same name
  EXPECT_THROW(ch2.queue_declare(QUEUE, false, true, true), AMQPChannelErr);

  // All further accesses to the illegal channel should be poisoned
  EXPECT_THROW(ch2.queue_declare(QUEUE".super-unique", false, true, true), AMQPChannelErr);

  // But the connection still be okay, and we should be able to grab another channel on it
  AmqpChannel ch3(c2);
  EXPECT_NO_THROW(ch3.queue_declare(QUEUE".unique1", false, true, true));
}

TEST_F(AMQPTest, SharedQueueServer) {
  Connection c1 = doconn(), c2 = doconn();

  if (!c1 || !c2)
    return;

  Server server1( new AmqpRpcServer(c1, EXCHANGE, QUEUE) ),
         server2;

  // Use the same queue twice. It should fail, since we're using exclusive queues
  EXPECT_THROW(server2.reset( new AmqpRpcServer(c2, EXCHANGE, QUEUE) ), AMQPChannelErr);
}
