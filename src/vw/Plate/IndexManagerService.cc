// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/IndexManagerService.h>
using namespace vw::platefile;

IndexManagerServiceImpl::IndexManagerServiceImpl() {

  // std::string queue_name = AmqpRpcClient::UniqueQueueName("index_manager");

  //m_client.reset  ( new AmqpRpcClient() );
  //m_index_service.reset ( new IndexService::Stub(
  //                          new AmqpRpcChannel(INDEX_EXCHANGE, "index", queue_name),
  //                          google::protobuf::Service::STUB_OWNS_CHANNEL) );
}

void IndexManagerServiceImpl::ListRequest(::google::protobuf::RpcController* controller,
                                          const IndexListRequest* request,
                                          IndexListReply* response,
                                          ::google::protobuf::Closure* done) {

#if 0

  std::string queue_name = AmqpClient::UniqueQueueName("index_mgr");

  AmqpConnection conn;
  m_conn.exchange_declare("index_mgr", "fanout", true, false);
  m_conn.queue_declare(queue_name, true, true, false);
  m_conn.queue_bind(queue_name, "index", "index");
#endif

#warning This is fake data
  response->add_platefile_id(0);
  response->add_platefile_id(1);
  done->Run();
}
