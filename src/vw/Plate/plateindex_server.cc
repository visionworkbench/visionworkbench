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
#include <vw/Plate/IndexService.h>

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

  po::options_description general_options("Runs a mosaicking daemon that listens for mosaicking requests coming in over the AMQP bus..\n\nGeneral Options:");
  general_options.add_options()
    ("queue_name,q", po::value<std::string>(&queue_name)->default_value("index"), 
     "Specify the name of the AMQP queue to create and listen to for mosaicking requests. (Defaults to the \"ngt_mosaic_worker\" queue.")
    ("help", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("root-directory", po::value<std::string>(&root_directory));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("root-directory", -1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [-q <queue name>] root_directory" <<std::endl << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  if( vm.count("root-directory") != 1 ) {
    std::cerr << "Error: must specify a root directory that contains plate files!" 
              << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }
  
  AmqpRpcServer server(INDEX_EXCHANGE, queue_name);
  boost::shared_ptr<google::protobuf::Service> service( new IndexServiceImpl(root_directory) );
  server.export_with_routing_key(service, "index");
  
  std::cout << "\t--> Listening for messages.\n";
  server.run();
  
  return 0;
}

