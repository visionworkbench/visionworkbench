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

TEST(AMQP, BasicIO) {
  boost::shared_ptr<AmqpConnection> conn;
  try {
    conn.reset(new AmqpConnection());
  } catch (const IOErr& e) {
#warning Need a better solution for this
    // If we can't open the AMQP socket, treat this as a disabled test.
    std::cerr << "Could not open AMQP socket. Test disabled." << std::endl;
    return;
  }
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

TEST(AMQP, RPCSerialization) {
  boost::shared_ptr<AmqpConnection> c1, c2;
  try {
    c1.reset( new AmqpConnection() );
    c2.reset( new AmqpConnection() );
  } catch (const IOErr& e) {
#warning Need a better solution for this
    std::cerr << "Could not open AMQP socket. Test disabled." << std::endl;
    // If we can't open the AMQP socket, treat this as a disabled test.
    return;
  }

  boost::shared_ptr<AmqpRpcServer> server( new AmqpRpcServer(c1, EXCHANGE, QUEUE"1") );
  boost::shared_ptr<AmqpRpcClient> client( new AmqpRpcClient(c2, EXCHANGE, QUEUE"2", "server") );

  boost::shared_ptr<google::protobuf::Service> srv_s( new IndexServiceImpl(TEST_SRCDIR) );
  boost::shared_ptr<google::protobuf::Service> srv_c( new IndexService::Stub(client.get() ) );

  server->bind_service(srv_s, "server");
  client->bind_service(srv_c, "client");

  IndexListRequest r1, r2;
  r1.set_type("toast");

  ASSERT_NO_THROW(client->send_message(r1, "server"));
  ASSERT_NO_THROW(server->get_message(r2));

  EXPECT_EQ(r1.type(), r2.type());
  EXPECT_EQ(r1.ByteSize(), r2.ByteSize());
}
