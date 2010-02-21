// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/Amqp.h>
#include <vw/Plate/RpcServices.h>
#include <vw/Plate/IndexService.h>
#include <vw/Plate/common.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Core/ThreadQueue.h>
#include <csignal>

using namespace vw;

#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace vw::platefile;
using namespace vw;

// A dummy method for passing to the RPC calls below.
static void null_closure() {}

// -----------------------------------------------------------------------------
//                                  MAIN
// -----------------------------------------------------------------------------

int main(int argc, char** argv) {

  po::options_description general_options("AMQP Performance Test Program");
  general_options.add_options()
    ("help", "Display this help message");

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(general_options).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << "\n\n";
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  std::string queue_name = "index_perftest_queue";

  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection());
  boost::shared_ptr<AmqpRpcClient> rpc_controller( new AmqpRpcClient(conn, DEV_INDEX,
                                                                     queue_name, "index") );
  boost::shared_ptr<IndexService> index_service( new IndexService::Stub(rpc_controller.get() ) );
  rpc_controller->bind_service(index_service, queue_name);
  
  int32 i = 0;
  while (1) {
    IndexTestRequest request;
    request.set_value(i);

    IndexTestReply response;
    index_service->TestRequest(rpc_controller.get(), &request, &response, 
                               google::protobuf::NewCallback(&null_closure));

    if (i != response.value())
      std::cout << "Error: IndexTestMessage failed!\n";

    ++i;
  }
}

