// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/Index.h>
#include <vw/Plate/Amqp.h>
#include <vw/Plate/ProtoBuffers.pb.h>
#include <vw/Plate/common.h>
#include <vw/Plate/RpcServices.h>
#include <vw/Plate/IndexManagerService.h>

#include <google/protobuf/descriptor.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace vw::platefile;
using namespace vw;

// -----------------------------------------------------------------------------
//                                  MAIN
// -----------------------------------------------------------------------------

int main(int argc, char** argv) {
  std::string queue_name, root_directory;

  po::options_description general_options("Handles queries that do index server discovery and inquery");
  general_options.add_options()
    ("help", "Display this help message");

  po::options_description options("Allowed Options");
  options.add(general_options);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << std::endl << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  AmqpRpcServer server(INDEX_MGR_EXCHANGE, "index_mgr");
  boost::shared_ptr<google::protobuf::Service> service( new IndexManagerServiceImpl() );
  server.export_with_routing_key(service, "index_mgr");

  std::cout << "\t--> Listening for messages.\n";
  server.run();

  return 0;
}

