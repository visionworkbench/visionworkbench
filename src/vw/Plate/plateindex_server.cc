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
#include <boost/regex.hpp>
namespace po = boost::program_options;
using namespace vw::platefile;
using namespace vw;

std::vector<std::string> plate_filenames(std::string const& root_directory)  {

  std::vector<std::string> result;

  if ( !fs::exists( root_directory ) ) 
    vw_throw(IOErr() << "Could not access the root directory: \"" << root_directory << "\"");

  // Create a regular expression for matching the pattern 'plate_<number>.blob"
  boost::regex re;
  re.assign(".*\\.plate", boost::regex_constants::icase);

  fs::directory_iterator end_itr; // default construction yields past-the-end

  // Iterate through the files in the platefile directory and return
  // any that match the regex above.
  for ( fs::directory_iterator itr( root_directory ); itr != end_itr; ++itr ) {
    if (boost::regex_match(itr->leaf(), re))
      result.push_back(itr->leaf());
  }

  return result;
}

int main(int argc, char** argv) {
  std::string queue_name, root_directory;

  po::options_description general_options("Runs a mosaicking daemon that listens for mosaicking requests coming in over the AMQP bus..\n\nGeneral Options:");
  general_options.add_options()
    ("queue_name,q", po::value<std::string>(&queue_name)->default_value(""), "Specify the name of the AMQP queue to create and listen to for mosaicking requests. (Defaults to the \"ngt_mosaic_worker\" queue.")
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

  std::cout << "Starting Index Server\n";
  std::vector<std::string> platefiles = plate_filenames(root_directory);
  if (platefiles.size() < 1) {
    std::cout << "Error: could not find any platefiles in the root directory.\nExiting.\n\n";
    exit(1);
  }
  
  try {
    for (int i = 0 ; i<platefiles.size(); ++i) {
      std::cout << "\t--> Loading platefile: " << root_directory << "/" << platefiles[i] << "\n";
    }    

    Index idx(platefiles[0]);

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

      if (routing_key == "index.index_open_request") {
        IndexOpenRequest r;
        r.ParseFromString(message);
        vw_out(InfoMessage, "plate") << "Processing index.index_open_request: "
                                     << r.DebugString() << "\n";
        
        // For now we return a dummy response, since we already have
        // the index open.
        IndexOpenReply reply;
        reply.set_platefile_id(0);
        reply.set_secret(0);
        conn.basic_publish_protobuf(reply, INDEX_EXCHANGE, 
                                    r.requestor()+".index_open_reply");


      } else if (routing_key == "index.read_request") {
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

  //     } else if (routing_key == "index.write_request") {

  //       IndexWriteRequest r;
  //       r.ParseFromString(message);
  //       vw_out(InfoMessage, "plate") << "Processing index.write_request: "
  //                                    << r.DebugString() << "\n";

  //       int32 blob_id = idx.write_request(r.size());
  //       IndexWriteReply response;
  //       response.set_blob_id(blob_id);
  //       conn.basic_publish_protobuf(response, INDEX_EXCHANGE,
  //                                   r.requestor() + ".write_reply");

  //     } else if (routing_key == "index.write_complete") {

  //       IndexWriteComplete r;
  //       r.ParseFromString(message);
  //       vw_out(InfoMessage, "plate") << "Processing index.write_complete: "
  //                                    << r.DebugString() << "\n";

  //       // Access the data in the index
  //       idx.write_complete(r.header(), r.record());
  //       IndexSuccess response;
  //       conn.basic_publish_protobuf(response, INDEX_EXCHANGE,
  //                                   r.requestor() + ".write_complete");

      } else {
        std::cout << "WARNING: unknown routing key \"" << routing_key 
                  << "\".  Ignoring request.\n";
      }
    }
  } catch (vw::Exception &e) {
    std::cout << "A server error occurred: " << e.what() << "\nExiting.\n\n";
    exit(0);
  } 
  return 0;
}

