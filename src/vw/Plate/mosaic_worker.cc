// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Plate/Index.h>
#include <vw/Plate/Amqp.h>
#include <vw/Plate/common.h>

// Protocols
#include <vw/Plate/IndexRpc.pb.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;
using namespace vw::platefile;

class RemoteIndex : public IndexBase {

  std::string m_platefile;
  std::string m_requestor;
  AmqpConnection m_conn;

public:
  /// Constructor
  RemoteIndex(std::string const& platefile, std::string const& requestor) :
    m_platefile(platefile), m_requestor(requestor) {}
  
  /// Destructor
  virtual ~RemoteIndex() {}
  
  /// Attempt to access a tile in the index.  Throws an
  /// TileNotFoundErr if the tile cannot be found.
  virtual IndexRecord read_request(int col, int row, int depth) {
    IndexReadRequest req;
    req.set_col(col);
    req.set_row(col);
    req.set_depth(col);
    req.set_requestor(m_requestor);
    req.set_platefile(m_platefile);
    m_conn.basic_publish_protobuf(req, EXCHANGE, m_requestor);
    
    std::string routing_key;
    std::string response = m_conn.basic_consume(QUEUE, routing_key, false);
    if (routing_key != m_requestor + ".read_response")
      vw_throw(LogicErr() << "Was expcting read response, but received " << routing_key << "instead.");

    IndexReadResponse r;
    r.ParseFromString(response);
    std::cout << r.DebugString() << "\n";
    return r.index_record();    
  }
  
  // Writing, pt. 1: Locks a blob and returns the blob id that can
  // be used to write a block.
  virtual int write_request(int size) = 0;
  
  // Writing, pt. 2: Supply information to update the index and
  // unlock the blob id.
  virtual void write_complete(int col, int row, int depth, IndexRecord record) = 0;
  
  
  virtual void save(std::string const& filename) = 0;
  
  
  virtual int32 version() const = 0;
  virtual int32 default_block_size() const = 0;
  virtual std::string default_block_filetype() const = 0;
  virtual int32 max_depth() const = 0;
};



// -----------------------------------------------------------------------------

int main(int argc, char** argv) {
  std::string queue_name;

  po::options_description general_options("Runs a mosaicking daemon that listens for mosaicking requests coming in over the AMQP bus..\n\nGeneral Options:");
  general_options.add_options()
    ("queue_name,q", po::value<std::string>(&queue_name)->default_value(""), "Specify the name of the AMQP queue to create and listen to for mosaicking requests. (Defaults to the \"ngt_mosaic_worker\" queue.")
    ("help", "Display this help message");

  po::options_description options("Allowed Options");
  options.add(general_options);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [-q <queue name>]" <<std::endl << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  std::cout << "Opening AMQP connection.\n";
  AmqpConnection conn;
  
  conn.exchange_declare("ngt", "direct", true);
  conn.basic_publish("hello there", "ngt", "jobs");

  conn.queue_declare("test", true, true, true);
  conn.queue_bind("test", "ngt", "jobs");
  //  while(1) 
    //    conn.basic_consume("test", "testtag", false);

  return 0;
}
