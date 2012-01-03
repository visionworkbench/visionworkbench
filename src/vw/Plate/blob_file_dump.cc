// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/Blob.h>
#include <vw/Plate/FundamentalTypes.h>

using namespace vw;
using namespace vw::platefile;

#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
namespace po = boost::program_options;

int main( int argc, char *argv[] ) {

  std::string blob_name;
  TransactionOrNeg transaction_id;

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
  } catch (const po::error& e) {
    std::cout << "An error occured while parsing command line arguments.\n\n";
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 1;
  }

  if( blob_name.empty() ) {
    std::cerr << "Error: must specify a platefile url.\n";
    std::cout << usage.str();
    return 1;
  }

  // -------------------

  ReadBlob blob(blob_name);

  size_t dot = blob_name.rfind('.');
  std::string blob_prefix = blob_name.substr(0,dot);

  std::cout << "Started!\n";

  BOOST_FOREACH(const BlobTileRecord& rec, blob) {
    if (transaction_id >= 0 && rec.hdr.transaction_id() != transaction_id)
      continue;
    std::ostringstream ostr;
    ostr << blob_prefix << "_" << rec.hdr.level() << "_" << rec.hdr.row() << "_" << rec.hdr.col() << "." << rec.hdr.filetype();

    std::ofstream ofile(ostr.str().c_str(), std::ios::binary);
    VW_ASSERT(ofile.is_open(), IOErr() << "could not open dst file for writing (" << ostr.str() << ")");
    ofile.write(reinterpret_cast<char*>(&rec.data->operator[](0)), rec.data->size());
    VW_ASSERT(!ofile.fail(), IOErr() << ": failed to write to " << ostr.str());
    ofile.close();
  }

  std::cout << "Finished!\n";

  return 0;
}
