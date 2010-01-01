#include <iostream>
#include <vw/Plate/RpcServices.h>
#include <vw/Plate/Amqp.h>
#include <vw/Core/ProgressCallback.h>

using namespace vw;
using namespace vw::platefile;

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

#if 0
  SharedByteArray result;
      client.get_bytes(result, 3000);
      server.send_bytes(msg, client_queue);
    return 0;
}
#endif

struct Options {
  std::string filename;
  uint32 clients;
  uint32 id;
  uint32 timeout;
  uint32 block_size;
  uint32 block_count;
  bool verify;
};

VW_DEFINE_EXCEPTION(Usage, Exception);

void handle_arguments(int argc, char *argv[], Options& opt)
{
  bool help = false;

  po::options_description options("Options");
  options.add_options()
    ("file,f",        po::value(&opt.filename), "File to append to")
    ("num-clients,n", po::value(&opt.clients),  "Number of clients.")
    ("id,i",          po::value(&opt.id),       "This client's id.")
    ("timeout,t",     po::value(&opt.timeout)->default_value(5000), "timeout in ms.")
    ("verify,v",      "Check the file instead of writing it")
    ("block-size,b",  po::value(&opt.block_size)->default_value(4),  "How much should each client write?")
    ("block-count,c", po::value(&opt.block_count)->default_value(std::numeric_limits<uint32>().max()), "How many times should each client write?")
    ("help,h",        "Display this help message.");

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).run(), vm );
  po::notify( vm );

  std::ostringstream out;

#define REQUIRE(name) do {\
  if (!vm.count(#name)) {\
    out << "Error: missing required parameter \"" #name "\"" << std::endl;\
    help = true;\
  }\
} while (0)

  if (vm.count("verify"))
    opt.verify = 1;
  if (vm.count("help"))
    help = true;

  if (!opt.verify)
    REQUIRE(id);

  REQUIRE(num-clients);
  REQUIRE(file);

#undef REQUIRE

  if( help )
    vw_throw(Usage() << out.str() << "Usage: " << argv[0] << "\n" << options << "\n");

  VW_ASSERT( opt.clients > 0,      Usage() << "Error: # of clients must be > 0."     );
  VW_ASSERT( opt.id < opt.clients, Usage() << "Error: # of clients must be > id."    );
  VW_ASSERT( opt.clients < 10000,  Usage() << "Error: # of clients must be < 10000." );
}

std::string queue_name(const Options& opt) {
  return std::string("client_") + vw::stringify(opt.id);
}
std::string next_queue_name(const Options& opt) {
  return std::string("client_") + vw::stringify((opt.id+1) % opt.clients );
}

void run(const Options& opt) {
  // On first pass with client id 0, there's no message to read.
  bool first_node = (opt.id == 0);

  boost::shared_ptr<AmqpConnection> conn(new AmqpConnection());
  AmqpRpcDumbClient node(conn, "lustre_torture", queue_name(opt));
  node.bind_service(queue_name(opt));

  std::string next = next_queue_name(opt);
  ByteArray send_msg(next.size() + 2);
  send_msg[0] = 'G';
  send_msg[1] = 'O';
  std::copy(next.begin(), next.end(), send_msg.begin()+2);

  TerminalProgressCallback pc;

  vw_out(InfoMessage) << "Ready to create "
                      << opt.clients * opt.block_size * opt.block_count / 1024. / 1024.
                      << "MB file! Press enter." << std::endl;
  getchar();

  uint32 blocks_written = 0;
  uint32 tick = opt.block_count / 100;
  while (blocks_written++ < opt.block_count) {
    if (blocks_written % tick == 0) {
      pc.report_incremental_progress(.01);
      pc.print_progress();
    }
    if (first_node) {
      first_node = false;
      // truncate file (and close it implicitly)
      std::ofstream file(opt.filename.c_str(), std::ios::binary|std::ios::trunc);
    } else {
      SharedByteArray msg;
      node.get_bytes(msg, opt.timeout);
      VW_ASSERT(msg->begin()[0] == 'G' && msg->begin()[1] == 'O', LogicErr() << "Buh?");
    }

    std::ofstream file(opt.filename.c_str(), std::ios::binary|std::ios::app);
    VW_ASSERT(file.is_open(), LogicErr() << "Could not open file: " << opt.filename);

    std::string s_id = boost::lexical_cast<std::string>(opt.id);
    s_id.insert(0, opt.block_size-s_id.size(), '0');
    VW_ASSERT(s_id.size() == opt.block_size, LogicErr() << "Must pad s_id to block size");

    file.write(s_id.c_str(), s_id.size());
    file.close();
    // to exacerbate the problem, don't put anything, put the close and the
    // unlock-next-client as close as possible
    node.send_bytes(send_msg, next);
  }
  pc.report_finished();
}

void verify(const Options& opt) {
    std::ifstream file(opt.filename.c_str(), std::ios::binary);
    VW_ASSERT(file.is_open(), LogicErr() << "Could not open file: " << opt.filename);

    uint32 id = 0, next;
    std::vector<char> bytes(opt.block_size+1, 0);
    uint32 failed = 0;

    TerminalProgressCallback pc;
    uint32 tick = (opt.clients * opt.block_count) / 100;

    uint32 blocks_read = 0;
    while (file.read(&bytes[0], opt.block_size)) {
      if (blocks_read++ % tick == 0) {
        pc.report_incremental_progress(.01);
        pc.print_progress();
      }
      try {
        next = boost::lexical_cast<uint32>(&bytes[0]);
      } catch (const boost::bad_lexical_cast&) {
        next = -1;
      }
      if (next != id)
        failed++;
      id = (id + 1) % opt.clients;
    }
    pc.report_finished();

    VW_ASSERT(failed == 0, LogicErr() << "Verification failed on " << failed << " blocks");
}

int main(int argc, char *argv[]) {
  Options opt;
  try {
    handle_arguments(argc, argv, opt);
    if (opt.verify)
      verify(opt);
    else
      run(opt);
  } catch (const Usage& e) {
    std::cout << e.what() << std::endl;
    return 1;
  } catch (const AMQPTimeout&) {
    std::cout << "Timeout!" << std::endl;
    return 0;
  } catch (const Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
