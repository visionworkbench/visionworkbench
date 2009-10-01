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

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <filename>..." <<std::endl << std::endl;
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

  if( output_file_name == "" )
    output_file_name = prefix_from_filename(plate_file_name) + ".toast";

  // Open the plate file
  PlateFile platefile(plate_file_name);
  std::cout << "Writing " << platefile.depth() << " levels of tiles to " 
            << output_file_name << "\n";


  // For debugging:
  //  platefile.print();

  // Create the output directory
  if ( !fs::exists(output_file_name) )
    fs::create_directory(output_file_name);

  for (int n = 0; n <= platefile.depth(); ++n) {
    int block_cols = pow(2,n);
    int block_rows = pow(2,n);
    
    for (int j = 0; j < block_rows; ++j) {
      for (int i = 0; i < block_cols; ++i) {
        
        ImageView<PixelRGB<uint8> > tile;
        try {
          IndexRecord rec = platefile.read(tile, i, j, n);
          if (!rec.valid())
            vw_throw(TileNotFoundErr() << "\tTile was found, but was marked invalid.");
          
          // Create the level directory (if it doesn't exist)
          std::ostringstream ostr;
          ostr << output_file_name << "/" << n;
          if ( !fs::exists(ostr.str()) )
            fs::create_directory(ostr.str());
          
          // Create the column directory (if it doesn't exist)
          ostr << "/" << i;
          if ( !fs::exists(ostr.str()) )
            fs::create_directory(ostr.str());
          
          // Create the file (with the row as the filename)
          ostr << "/" << j << "." << rec.block_filetype();
          std::cout << "\t--> Writing " << ostr.str() << "\n";
          write_image(ostr.str(), tile);

        } catch (TileNotFoundErr &e) {}
      }
    }    

    platefile.save();

  }
}
