// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/RpcServices.h>
#include <vw/Plate/Amqp.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Core/ThreadQueue.h>
#include <csignal>

using namespace vw;

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace vw::platefile;
using namespace vw;

// --------------------------------------------------------------
//                           CLIENT
// --------------------------------------------------------------

void run_client(std::string exchange, std::string client_queue,
                const std::string server_queue, unsigned message_size, bool rtt, unsigned exchanges) {
  std::cerr << "Running client...\n";

  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection());
  AmqpRpcDumbClient client(conn, exchange, client_queue, exchanges);
  client.bind_service(client_queue);

  int msgs = 0;
  long long t0 = Stopwatch::microtime();

  uint8 frame_seq = 0;
  SharedByteArray result;
  while (1) {
    try {
      client.get_bytes(result, 3000);
    } catch (const AMQPTimeout&) {
      vw_out(0) << "No messages for 3 seconds" << std::endl;
      break;
    }

    if (result->size() != message_size) {
      std::cerr << "Error -- unexpected message size : " << result->size() << "\n";
    }

    for (unsigned i = 1; i < result->size(); ++i) {
      bool bad = false;
      if ( (*result)[i] != i % 256 ) {
        std::cerr << "      Bad byte: " << int((*result)[i]) << "\n";
        bad = true;
      }

      if (bad)
        std::cerr << "Error -- corrupt message payload.\n";

    }

    if ((*result)[0] != frame_seq) {
      std::cerr << "Out of order message [expected " << uint32(frame_seq) << " got " << uint32((*result)[0]) << "]" << std::endl;
      frame_seq = (*result)[0];
    } else {
      frame_seq++;
    }

    // Echo the message back to the server
    if (rtt)
      client.send_bytes(*(result.get()), server_queue);

    ++msgs;

    float diff = (Stopwatch::microtime() - t0) / 1e6;
    if (diff > 1.0) {
      std::cout << "messages / second : " << (msgs/diff) << "\n";
      t0 = Stopwatch::microtime();
      msgs = 0;
    }
  }
}

volatile bool go = true;

void sighandle(int sig) {
  signal(sig, SIG_IGN);
  go = false;
  signal(sig, sighandle);
}

// --------------------------------------------------------------
//                           SERVER
// --------------------------------------------------------------

void run_server(const std::string exchange, const std::string client_queue,
                const std::string server_queue, unsigned message_size, bool rtt, unsigned batch, unsigned exchanges) {
  std::cerr << "Running server...\n";

  signal(SIGINT, sighandle);

  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection());
  AmqpRpcDumbClient server(conn, exchange, server_queue, exchanges);
  server.bind_service(server_queue);

  // Set up the message
  ByteArray msg(message_size);
  uint8 frame_seq = 0;
  for (unsigned i = 1; i < message_size; ++i) {
    msg[i] = i % 256;
  }

  SharedByteArray result;
  long long t0 = Stopwatch::microtime();
  long long msgs = 0;

  while(go) {
    for (unsigned i = 0; i < batch; ++i) {
      msg[0] = frame_seq++;
      server.send_bytes(msg, client_queue);
      msgs++;
    }
    float diff = (Stopwatch::microtime() - t0) / 1e6;
    if (diff > 1.0) {
      std::cout << "sent messages / second : " << (msgs/diff) << "\n";
      t0 = Stopwatch::microtime();
      msgs = 0;
    }

    if (rtt) {
      try {
        server.get_bytes(result, 3000);
      } catch (const AMQPTimeout&) {
        vw_out(0) << "No ACKs for 3 seconds" << std::endl;
        break;
      }

      // Verify message
      if (result->size() != msg.size()) {
        std::cerr << "Error -- unexpected echo message size : " << result->size() << "\n";
      }

      for (unsigned i = 1; i < result->size(); ++i) {
        bool bad = false;
        if ( (*result)[i] != msg[i] ) {
          std::cerr << "      Bad byte: " << int((*result)[i]) << "\n";
          bad = true;
        }

        if (bad)
          std::cerr << "Error -- corrupt message payload.\n";
      }
    }
  }
}

// -----------------------------------------------------------------------------
//                                  MAIN
// -----------------------------------------------------------------------------

int main(int argc, char** argv) {

  unsigned message_size;
  unsigned batch;
  unsigned exchanges;

  po::options_description general_options("AMQP Performance Test Program");
  general_options.add_options()
    ("client", "Act as client.")
    ("server", "Act as server.")
    ("batch",     po::value(&batch)->default_value(1), "How many messages to publish before checking for responses (server-only)")
    ("exchanges", po::value(&exchanges)->default_value(1), "How many exchanges should we stripe across?")
    ("rtt", "Measure round trip messages/sec instead of one way messages/sec.")
    ("message-size", po::value(&message_size)->default_value(5), "Message size in bytes")
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
    run_server("ptest_exchange", "ptest_client_queue", "ptest_server_queue",
               message_size, vm.count("rtt"), batch, exchanges);

  if( vm.count("client") )
    run_client("ptest_exchange", "ptest_client_queue", "ptest_server_queue",
               message_size, vm.count("rtt"), exchanges);

  return 0;
}

