// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Core/Stopwatch.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/AmqpConnection.h>
#include <vw/Plate/ProtoBuffers.pb.h>
#include <vw/Plate/common.h>
#include <vw/Plate/RpcServices.h>
#include <vw/Plate/WMSService.h>
#include <signal.h>

#include <google/protobuf/descriptor.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

using namespace vw::platefile;
using namespace vw;

#define EXCHANGE PLATE_EXCHANGE_NAMESPACE ".wms"
VW_DEFINE_EXCEPTION(Usage, Exception);

// Global variable.  (Makes signal handling a lot easier...)
boost::shared_ptr<WMSServiceImpl> g_service;

// ------------------------------ SIGNAL HANDLER -------------------------------

volatile bool process_messages = true;

void sig_unexpected_shutdown(int sig_num) {
  signal(sig_num, SIG_IGN);
  process_messages = false;
  signal(sig_num, sig_unexpected_shutdown);
}

// ------------------------------    MAIN LOOP   -------------------------------

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

// ------------------------------      MAIN      -------------------------------

struct Options {
  std::string exchange_name;
  std::string platefile_root;
  std::string cache_root;
  std::string rabbit_hostname;
  int rabbit_port;
  bool debug;
};

void handle_arguments(int argc, char *argv[], Options& opt) {
  po::options_description general_options("Runs a daemon that listens for wms requests coming in over the AMQP bus.\n\nGeneral Options:");
  general_options.add_options()
    ("exchange,e", po::value(&opt.exchange_name)->default_value(EXCHANGE),
     "Specify the name of the AMQP exchange to use for the wms service.")
    ("hostname,h", po::value(&opt.rabbit_hostname)->default_value("localhost"),
     "Specify the hostname of the AMQP server to use for remote procedure calls (RPCs).")
    ("port,p", po::value(&opt.rabbit_port)->default_value(5672),
     "Specify the port of the AMQP server to use for remote procedure calls (RPCs).")
    ("cache-root", po::value(&opt.cache_root)->default_value("/tmp/platefile_cache"),
     "Specify the cache directory. Both the server and the client need to be able to see this.")
    ("debug", "Output debug messages.")
    ("help", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("root-directory", po::value(&opt.platefile_root));

  opt.platefile_root = fs::system_complete(opt.platefile_root).string();
  opt.cache_root     = fs::system_complete(opt.cache_root).string();

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("root-directory", -1);

  po::variables_map vm;

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " platefile_root" << std::endl
        << std::endl << general_options << std::endl;

  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error &e) {
    vw_throw(Usage() << "Error parsing input:\n\t" << e.what() << "\n" << usage.str());
  }

  if( vm.count("help") )
    vw_throw(Usage() << usage.str());

  if( vm.count("root-directory") != 1 ) {
    vw_throw(Usage() << "Error: must specify a root directory that contains plate files!\n\n"
                     << usage.str());
  }

  opt.debug = vm.count("debug");
}

void run(const Options& opt) {
  std::string queue_name = AmqpRpcClient::UniqueQueueName("wms_server");

  boost::shared_ptr<AmqpConnection> connection( new AmqpConnection(opt.rabbit_hostname, opt.rabbit_port) );
  boost::shared_ptr<AmqpRpcServer> server( new AmqpRpcServer(connection, opt.exchange_name,
                                                             queue_name, opt.debug) );
  g_service.reset( new WMSServiceImpl(opt.platefile_root, opt.cache_root) );
  server->bind_service(g_service, "wms");

  // Start the server task in another thread
  boost::shared_ptr<ServerTask> server_task( new ServerTask(server) );
  Thread server_thread( server_task );

  // Install Unix Signal Handlers.
  signal(SIGINT, sig_unexpected_shutdown);

  std::cout << "\n\n";
  long long t0 = Stopwatch::microtime();

  while(process_messages) {

    int queries = server->queries_processed();
    size_t bytes = server->bytes_processed();
    int n_outstanding_messages = server->incoming_message_queue_size();
    server->reset_stats();

    float dt = float(Stopwatch::microtime() - t0) / 1e6;
    t0 = Stopwatch::microtime();

    std::cout << "[wms_server] : "
              << float(queries/dt) << " qps    "
              << float(bytes/dt)/1000.0 << " kB/sec    "
              << n_outstanding_messages << " outstanding messages                          \r"
              << std::flush;
    sleep(1.0);
  }

  std::cout << "\nShutting down!" << std::endl;
}

int main(int argc, char** argv) {
  Options opt;
  try {
    handle_arguments(argc, argv, opt);
    run(opt);
  } catch (const Usage& e) {
    std::cout << e.what() << std::endl;
    return 1;
  } catch (const Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}


