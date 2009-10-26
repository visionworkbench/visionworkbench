// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/Index.h>
#include <vw/Plate/Amqp.h>
#include <vw/Plate/ProtoBuffers.pb.h>
#include <vw/Plate/common.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace vw::platefile;
using namespace vw;

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
    
    // Create a "topic" exchange on the AMQP broker to handle index
    // server requests.  This, in conjunction with the queue and
    // binding below, will match any messages that are sent using the
    // 'index.#' routing key.  The '.*' portion of the routing
    // key will be used later to determine what message type was sent.
    std::cout << "\t--> Declaring \'"<< INDEX_EXCHANGE <<"\'exchange binding it to the \'"
              << INDEX_QUEUE <<"' queue.\n";
    conn.exchange_declare(INDEX_EXCHANGE, "topic", true, false);
    conn.queue_declare(INDEX_QUEUE, true, true, false);
    conn.queue_bind(INDEX_QUEUE, INDEX_EXCHANGE, "index.#");

    while(1) {
      std::string routing_key;
      std::string message = conn.basic_consume(INDEX_QUEUE, routing_key, false);

      if (routing_key == "index.read_request") {
        IndexReadRequest r;
        r.ParseFromString(message);
        vw_out(InfoMessage, "plate") << "Processing index.read_request: "
                                     << r.DebugString() << "\n";

        // Access the data in the index
        try {
          IndexRecord record = idx.read_request(r.col(), r.row(), r.depth(), r.transaction_id());
          IndexReadReply read_response;
          *(read_response.mutable_index_record()) = record;
          conn.basic_publish_protobuf(read_response, INDEX_EXCHANGE,
                                      r.requestor() + ".read_reply");
        } catch (TileNotFoundErr &e) {
          IndexError err; 
          err.set_message(e.what());
          conn.basic_publish_protobuf(err, INDEX_EXCHANGE, r.requestor() + ".index_error");
        }

      } else if (routing_key == "index.write_request") {

        IndexWriteRequest r;
        r.ParseFromString(message);
        vw_out(InfoMessage, "plate") << "Processing index.write_request: "
                                     << r.DebugString() << "\n";

        int32 blob_id = idx.write_request(r.size());
        IndexWriteReply response;
        response.set_blob_id(blob_id);
        conn.basic_publish_protobuf(response, INDEX_EXCHANGE,
                                    r.requestor() + ".write_reply");

      } else if (routing_key == "index.write_complete") {

        IndexWriteComplete r;
        r.ParseFromString(message);
        vw_out(InfoMessage, "plate") << "Processing index.write_complete: "
                                     << r.DebugString() << "\n";

        // Access the data in the index
        idx.write_complete(r.header(), r.record());
        IndexSuccess response;
        conn.basic_publish_protobuf(response, INDEX_EXCHANGE,
                                    r.requestor() + ".write_complete");

      } else {
        std::cout << "WARNING: unknown routing key \"" << routing_key 
                  << "\".  Ignoring request.\n";
      }
    }

  } catch (vw::IOErr &e) {
    std::cout << "A server error occurred: " << e.what() << "\nExiting.\n\n";
    exit(0);
  } 
  return 0;
}

