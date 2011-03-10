// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/Blob.h>

using namespace vw;
using namespace vw::platefile;

#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
namespace po = boost::program_options;

int main( int argc, char *argv[] ) {

  std::string blob_name;
  int transaction_id;

  po::options_description general_options("Writes out all inner images of a blob.\n\nGeneral Options");
  general_options.add_options()
    ("transaction_id,t", po::value(&transaction_id),
     "Transaction id to pull from blob.")
    ("help,h", "Display this help message.");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("blob", po::value(&blob_name));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("blob", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <blob name>..." <<std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    std::cout << "An error occured while parsing command line arguments.\n\n";
    std::cout << usage.str();
    return 0;
  }

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }

  if( blob_name.empty() ) {
    std::cerr << "Error: must specify a platefile url.\n";
    std::cout << usage.str();
    return 0;
  }

  // -------------------

  Blob blob( blob_name, true );

  size_t dot = blob_name.rfind('.');
  std::string blob_prefix = blob_name.substr(0,dot);

  std::cout << "Started!\n";

  for ( Blob::iterator it = blob.begin(); it != blob.end(); it++ ) {
    if ( transaction_id < 0 ||
         it->transaction_id() == uint64(transaction_id) ) {
      std::ostringstream ostr;
      ostr << blob_prefix << "_" << it.current_base_offset() << "_"
           << it->col() << "_" << it->row() << "_" << it->level()
           << "." << it->filetype();
      blob.read_to_file( ostr.str(), it.current_base_offset() );
    }
  }

  std::cout << "Finished!\n";
}
