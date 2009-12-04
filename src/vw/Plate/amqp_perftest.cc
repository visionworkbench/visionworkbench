// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/Amqp.h>
#include <vw/Core/Stopwatch.h>
using namespace vw;

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace vw::platefile;
using namespace vw;

void run_client(std::string exchange, std::string client_queue, bool use_get) {
  std::cout << "Running client...\n";

  AmqpConnection m_conn;

  // Create a queue and bind it to the index server exchange.
  m_conn.exchange_declare(exchange, "direct", false, false);  // not durable, no autodelete
  m_conn.queue_declare(client_queue, false, true, true);  // not durable, exclusive, auto-delete
  m_conn.queue_bind(client_queue, exchange, client_queue);
  
  int msgs = 0;
  long long t0 = Stopwatch::microtime();

  while (1) {
    boost::shared_array<uint8> result;
    if (use_get) {
      result = m_conn.basic_get(client_queue, true, -1);
    } else { 
      std::string routing_key;
      result = m_conn.basic_consume(client_queue, routing_key, true);
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


void run_server(std::string exchange, std::string client_queue, int msg_size) {
  std::cout << "Running server...\n";

  AmqpConnection m_conn;

  // Create a queue and bind it to the index server exchange.
  m_conn.exchange_declare(exchange, "direct", false, false);

  // Set up the message
  boost::shared_array<uint8> msg ( new uint8[msg_size] );

  while(1) {
    m_conn.basic_publish(msg, msg_size, exchange, client_queue);
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
    run_server("ptest_exchange", "ptest_queue", 4);

  if( vm.count("client") ) 
    run_client("ptest_exchange", "ptest_queue", vm.count("get"));

  return 0;
}

