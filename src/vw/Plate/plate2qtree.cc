// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

/// deplate.cc
///
/// Converts a plate file to a nested set of directories and tiles on
/// disk.

#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Plate/PlateFile.h>

using namespace vw;
using namespace vw::platefile;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

// Erases a file suffix if one exists and returns the base string
static std::string prefix_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1) 
    result.erase(index, result.size());
  return result;
}

class SaveTileFunc : public TreeMapFunc {
  boost::shared_ptr<PlateFile> m_platefile;
  std::string m_output_filename;

public:
  SaveTileFunc(boost::shared_ptr<PlateFile> platefile, std::string output_filename) : 
    m_platefile(platefile), m_output_filename(output_filename) {}
  virtual ~SaveTileFunc() {}

  virtual void operator() (int32 col, int32 row, int32 level) {

    try {
      
      // Create the level directory (if it doesn't exist)
      std::ostringstream ostr;
      ostr << m_output_filename << "/" << level;
      if ( !fs::exists(ostr.str()) )
        fs::create_directory(ostr.str());
      
      // Create the column directory (if it doesn't exist)
      ostr << "/" << col;
      if ( !fs::exists(ostr.str()) )
        fs::create_directory(ostr.str());
      
      // Create the file (with the row as the filename)
      ostr << "/" << row;

      std::string output_filename = m_platefile->read_to_file(ostr.str(), col, row, level);
      std::cout << "\t--> [ " << col << " " << row << " " << level << "] : Writing " 
                << output_filename << "\n";
      
    } catch (TileNotFoundErr &e) { 
      std::cout << "\t--> [ " << col << " " << row << " " << level << "] : Missing tile\n"; }
  }
};

int main( int argc, char *argv[] ) {
 
  std::string output_file_name;
  std::string plate_file_name;

  po::options_description general_options("Turns georeferenced image(s) into a TOAST quadtree.\n\nGeneral Options");
  general_options.add_options()
    ("output-name,o", po::value<std::string>(&output_file_name), "Specify the base output directory")
    ("help", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("plate-file", po::value<std::string>(&plate_file_name));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("plate-file", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <filename>..." <<std::endl << std::endl;
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

  if( vm.count("plate-file") != 1 ) {
    std::cerr << "Error: must specify an input platefile!" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  if( output_file_name == "" )
    output_file_name = prefix_from_filename(plate_file_name);

  // Open the plate file
  boost::shared_ptr<PlateFile> platefile(new PlateFile(plate_file_name));
  std::cout << "Writing " << platefile->depth() << " levels of tiles to " 
            << output_file_name << "\n";

  // Create the output directory
  if ( !fs::exists(output_file_name) )
    fs::create_directory(output_file_name);

  // Spider across the platefile, saving all valid tiles.
  boost::shared_ptr<TreeMapFunc> func(new SaveTileFunc(platefile, output_file_name));
  platefile->map(func);
}
