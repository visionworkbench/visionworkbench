// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/PlateFile.h>
#include <vw/Plate/LocalIndex.h>
using namespace vw;
using namespace vw::platefile;

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

// --------------------------------------------------------------------------
//                                    MAIN
// --------------------------------------------------------------------------

int main( int argc, char *argv[] ) {

  std::string filename;

  po::options_description general_options("\nRebuild a platefile index.\n");
  general_options.add_options()
    ("help,h", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::string>(&filename), "");

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <plate_filename>\n";
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    std::cout << "An error occured while parsing command line arguments.\n\n";
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 1;
  }

  if (filename.empty()) {
    std::cout << usage.str();
    return 1;
  }

  try {

    //--------------------------- OPEN THE PLATE FILE -----------------------------

    // Check to see whether the index exists:
    std::string index_str = filename + "/index";
    if (!fs::exists(filename)) {
      std::cout << "Error: could not open platefile: \"" << filename << "\".\n";
      exit(1);
    }

    if (fs::exists(index_str)) {
      std::cout << "Are you sure you want to delete & rebuild \"" << index_str << "\"? ";
      std::string user_input;
      std::cin >> user_input;
      if (user_input == "y")
        fs::remove_all(index_str);
      else
        exit(0);
    }

    LocalIndex index(filename);
    index.rebuild_index();

 }  catch (vw::Exception &e) {
    std::cout << "An error occured: " << e.what() << "\nExiting.\n\n";
    exit(1);
  }

}
