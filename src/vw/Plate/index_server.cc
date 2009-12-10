// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Core/Stopwatch.h>
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
//                                  MAIN LOOP
// -----------------------------------------------------------------------------

class ServerTask {
  boost::shared_ptr<AmqpRpcServer> m_server;

public:

  ServerTask(boost::shared_ptr<AmqpRpcServer> server) : m_server(server) {}

  void operator()() {
    std::cout << "\t--> Listening for messages.\n";
    m_server->run();
  }

  void kill() { m_server->shutdown(); }
};

// -----------------------------------------------------------------------------
//                                  MAIN
// -----------------------------------------------------------------------------

int main(int argc, char** argv) {
  std::string queue_name, root_directory;
  std::string hostname;
  int port;

  po::options_description general_options("Runs a mosaicking daemon that listens for mosaicking requests coming in over the AMQP bus..\n\nGeneral Options:");
  general_options.add_options()
    ("queue_name,q", po::value<std::string>(&queue_name)->default_value("index"),
     "Specify the name of the AMQP queue to create and listen to for mosaicking requests. (Defaults to the \"ngt_mosaic_worker\" queue.")
    ("hostname,h", po::value<std::string>(&hostname)->default_value("localhost"),
     "Specify the hostname of the AMQP server to use for remote procedure calls (RPCs).")
    ("port,p", po::value<int>(&port)->default_value(5672),
     "Specify the port of the AMQP server to use for remote procedure calls (RPCs).")
    ("debug", "Output debug messages.")
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

  boost::shared_ptr<AmqpConnection> connection( new AmqpConnection(hostname, port) );
  boost::shared_ptr<AmqpRpcServer> server( new AmqpRpcServer(connection, INDEX_EXCHANGE, 
                                                             queue_name, vm.count("debug")) );
  boost::shared_ptr<google::protobuf::Service> service( new IndexServiceImpl(root_directory) );
  server->bind_service(service, "index");

  // Start the server task in another thread
  boost::shared_ptr<ServerTask> server_task( new ServerTask(server) );
  Thread server_thread( server_task );

  long long t0 = Stopwatch::microtime();

  std::cout << "\n\n";
  while(1) {
    int queries = server->queries_processed();
    size_t bytes = server->bytes_processed();
    int n_outstanding_messages = server->incoming_message_queue_size();
    server->reset_stats();

    float dt = float(Stopwatch::microtime() - t0) / 1e6;
    t0 = Stopwatch::microtime();

    std::cout << "[index_server] : "
              << float(queries/dt) << " qps    "
              << float(bytes/dt)/1000.0 << " kB/sec    " 
              << n_outstanding_messages << " outstanding messages                          \r" 
              << std::flush;
    sleep(1.0);
  }

  return 0;
}

