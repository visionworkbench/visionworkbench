
#include <vw/Plate/Index.h>
#include <vw/Plate/Amqp.h>
#include <vw/Plate/IndexRpc.pb.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace vw::platefile;


#define EXCHANGE "ngt_platefile_index"
#define QUEUE "index"

int main(int argc, char** argv) {
  std::string queue_name, plate_file_name;

  po::options_description general_options("Runs a mosaicking daemon that listens for mosaicking requests coming in over the AMQP bus..\n\nGeneral Options:");
  general_options.add_options()
    ("queue_name,q", po::value<std::string>(&queue_name)->default_value(""), "Specify the name of the AMQP queue to create and listen to for mosaicking requests. (Defaults to the \"ngt_mosaic_worker\" queue.")
    ("help", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("plate-file", po::value<std::string>(&plate_file_name));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("plate-file", -1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [-q <queue name>]" <<std::endl << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  if( vm.count("plate-file") != 1 ) {
    std::cerr << "Error: must specify an input platefile!" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  std::cout << "Starting Index Server\n";
  try {

    std::cout << "\t--> Opening plate file: " << plate_file_name << "\n";
    Index idx(plate_file_name);
    
    std::cout << "\t--> Opening connection to AMQP broker\n";
    AmqpConnection conn;
    
    // Create a "topic" exchange on the AMQP broker.  This, in
    // conjunction with the queue and binding below, will match any
    // messages that are sent using the 'indexserver.*' routing key.
    // The routing key will be used later to determine what message type
    // was sent.
    std::cout << "\t--> Declaring exchange, queue, and binding.\n";
    conn.exchange_declare(EXCHANGE, "topic", true, false);
    conn.queue_declare(QUEUE, true, true, false);
    conn.queue_bind(QUEUE, EXCHANGE, "index.#");

    while(1) {
      std::string routing_key;
      std::string message = conn.basic_consume(QUEUE, routing_key, false);

      if (routing_key == "index.read_request") {
        std::cout << "Processing \"index.read_request\"\n";

        IndexReadRequest r;
        r.ParseFromString(message);
        std::cout << r.DebugString() << "\n";

        // Access the data in the index
        try {
          IndexRecord record = idx.read_request(r.col(), r.row(), r.depth());
          IndexReadResponse read_response;
          *(read_response.mutable_index_record()) = record;
          conn.basic_publish_protobuf(read_response, EXCHANGE, r.requestor() + ".read_response");
        } catch (TileNotFoundErr &e) {
          IndexError err; 
          err.set_message(e.what());
          conn.basic_publish_protobuf(err, EXCHANGE, r.requestor() + ".index_error");
        }

      } else if (routing_key == "index.write_request") {
        std::cout << "Processing \"index.write_request\"\n";


      } else if (routing_key == "index.write_complete") {
        std::cout << "Processing \"index.complete_request\"\n";


      } else {
        std::cout << "WARNING: unhandled message with routing key \"" << routing_key << "\".\n";
      }
    }

  } catch (vw::IOErr &e) {
    std::cout << "A server error occurred: " << e.what() << "\nExiting.\n\n";
    exit(0);
  } 
  return 0;
}

