// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/Amqp.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Core/ThreadQueue.h>
#include <csignal>

using namespace vw;

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace vw::platefile;
using namespace vw;

void run_client(std::string exchange, std::string client_queue) {
  std::cout << "Running client...\n";

  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection());
  AmqpChannel chan(conn);

  chan.exchange_declare(exchange, "direct", false, false);  // not durable, no autodelete
  chan.queue_declare(client_queue, false, true, true);  // not durable, exclusive, auto-delete
  chan.queue_bind(client_queue, exchange, client_queue);

  int msgs = 0;
  long long t0 = Stopwatch::microtime();


  ThreadQueue<SharedByteArray> q;
  boost::shared_ptr<AmqpConsumer> c = chan.basic_consume(client_queue, boost::bind(&ThreadQueue<SharedByteArray>::push, boost::ref(q), _1));

  SharedByteArray result;
  while (1) {
    if (!q.timed_wait_pop(result, 3000)) {
      vw_out(0) << "No messages for 5 seconds" << std::endl;
      break;
    }

    ++msgs;

    float diff = (Stopwatch::microtime() - t0) / 1e6;
    if (diff > 1.0) {
      std::cout << "messages / second : " << (msgs/diff) << "\n";
      t0 = Stopwatch::microtime();
      msgs = 0;
    }
  }
}

bool go = true;

void sighandle(int sig) {
  signal(sig, SIG_IGN);
  go = false;
  signal(sig, sighandle);
}

void run_server(const std::string exchange, const std::string client_queue, const std::string& msg_text) {
  std::cout << "Running server...\n";
  signal(SIGINT, sighandle);

  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection());
  AmqpChannel chan(conn);

  // Create a queue and bind it to the index server exchange.
  chan.exchange_declare(exchange, "direct", false, false);

  // Set up the message
  ByteArray msg(msg_text.size());
  std::copy(msg_text.begin(), msg_text.end(), msg.begin());

  while(go) {
    chan.basic_publish(msg, exchange, client_queue);
  }
}


// -----------------------------------------------------------------------------
//                                  MAIN
// -----------------------------------------------------------------------------

int main(int argc, char** argv) {

  po::options_description general_options("AMQP Performance Test Program");
  general_options.add_options()
    ("client", "Act as client.")
    ("server", "Act as server.")
    ("get", "Use basic_get() instead of basic_consume().")
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

  if( vm.count("server") )
    run_server("ptest_exchange", "ptest_queue", "test");

  if( vm.count("client") )
    run_client("ptest_exchange", "ptest_queue");

  return 0;
}

