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

// --------------------------------------------------------------
//                           CLIENT
// --------------------------------------------------------------

void run_client(std::string exchange, std::string client_queue,
                const std::string server_queue, int message_size, bool rtt) {
  std::cerr << "Running client...\n";

  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection());
  AmqpChannel chan(conn);

  chan.exchange_declare(exchange, "direct", false, false);  // not durable, no autodelete
  chan.queue_declare(client_queue, false, true, true);  // not durable, exclusive, auto-delete
  chan.queue_bind(client_queue, exchange, client_queue);

  int msgs = 0;
  long long t0 = Stopwatch::microtime();

  ThreadQueue<SharedByteArray> q;
  boost::shared_ptr<AmqpConsumer> c = chan.basic_consume(client_queue,
                                                         boost::bind(&ThreadQueue<SharedByteArray>::push,
                                                                     boost::ref(q), _1));

  uint8 frame_seq = 0;
  SharedByteArray result;
  while (1) {
    if (!q.timed_wait_pop(result, 3000)) {
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
      chan.basic_publish(*(result.get()), exchange, server_queue);

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
                const std::string server_queue, int message_size, bool rtt, unsigned batch) {
  std::cerr << "Running server...\n";

  signal(SIGINT, sighandle);

  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection());
  AmqpChannel chan(conn);

  // Create a queue and bind it to the index server exchange.
  chan.exchange_declare(exchange, "direct", false, false);
  chan.queue_declare(server_queue, false, true, true);  // not durable, exclusive, auto-delete
  chan.queue_bind(server_queue, exchange, server_queue);

  ThreadQueue<SharedByteArray> q;
  boost::shared_ptr<AmqpConsumer> c = chan.basic_consume(server_queue,
                                                         boost::bind(&ThreadQueue<SharedByteArray>::push,
                                                                     boost::ref(q), _1));

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
    for (int i = 0; i < batch; ++i) {
      msg[0] = frame_seq++;
      chan.basic_publish(msg, exchange, client_queue);
      msgs++;
    }
    float diff = (Stopwatch::microtime() - t0) / 1e6;
    if (diff > 1.0) {
      std::cout << "sent messages / second : " << (msgs/diff) << "\n";
      t0 = Stopwatch::microtime();
      msgs = 0;
    }

    if (rtt) {
      if (!q.timed_wait_pop(result, 3000)) {
        vw_out(0) << "No ACKs for 3 seconds" << std::endl;
        break;
      }

      // Verify message
      if (result->size() != msg.size()) {
        std::cerr << "Error -- unexpected echo message size : " << result->size() << "\n";
      }

      for (int i = 1; i < result->size(); ++i) {
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

  int message_size;
  unsigned batch;

  po::options_description general_options("AMQP Performance Test Program");
  general_options.add_options()
    ("client", "Act as client.")
    ("server", "Act as server.")
    ("batch",  po::value(&batch)->default_value(1), "How many messages to publish before checking for responses (server-only)")
    ("rtt", "Measure round trip messages/sec instead of one way messages/sec.")
    ("message-size", po::value<int>(&message_size)->default_value(5), "Message size in bytes")
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
               message_size, vm.count("rtt"), batch);

  if( vm.count("client") )
    run_client("ptest_exchange", "ptest_client_queue", "ptest_server_queue",
               message_size, vm.count("rtt"));

  return 0;
}

