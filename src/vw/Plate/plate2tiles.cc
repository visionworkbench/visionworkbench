// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// Converts a plate file to a nested set of directories and tiles on
/// disk.

#include <sstream>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/TileManipulation.h>

using namespace vw;
using namespace vw::platefile;

#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

// Erases a file suffix if one exists and returns the base string
static std::string prefix_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1)
    result.erase(index, result.size());
  return result;
}
// ----------------------------------------------------------------------------
//                          save_toast_tile()
// ----------------------------------------------------------------------------

void save_toast_tile(std::string base_output_name, boost::shared_ptr<PlateFile> platefile,
                     int32 col, int32 row, int32 level, int32 transaction_id) {

  try {

    // Create the level directory (if it doesn't exist)
    std::ostringstream ostr;
    ostr << base_output_name << "/" << level;
    if ( !fs::exists(ostr.str()) )
      fs::create_directory(ostr.str());

    // Create the column directory (if it doesn't exist)
    ostr << "/" << col;
    if ( !fs::exists(ostr.str()) )
      fs::create_directory(ostr.str());

    // Create the file (with the row as the filename)
    ostr << "/" << row;

    // transaction_id = -1 returns the latest tile available
    std::string output_filename = platefile->read_to_file(ostr.str(), col, row, level, transaction_id);
    // std::cout << "\t--> [ " << col << " " << row << " " << level << "] : Writing "
    //           << output_filename << "\n";

  } catch (TileNotFoundErr &e) {
    //    std::cout << "\t--> [ " << col << " " << row << " " << level << "] : Missing tile\n";
  }
}

// ----------------------------------------------------------------------------
//                          save_gigapan_tile()
// ----------------------------------------------------------------------------

void save_gigapan_tile(std::string base_output_name, boost::shared_ptr<PlateFile> platefile,
                       int32 col, int32 row, int32 level, int32 transaction_id) {

  try {
    std::stringstream filename_stream;
    std::stringstream directory_stream;

    directory_stream << base_output_name << '/';

    filename_stream << 'r';

    for (int32 l = level - 1; l >= 0; l--) {
      uint32 bit = 1 << l;
      int index = 0;
      if ( col & bit )
        index = 1;
      if ( row & bit )
        index += 2;

      filename_stream << index;
      }

    std::string filename = filename_stream.str();

    int size = filename.size();
    while ( (size > 3) && filename.size() >= 3 ) {
      directory_stream << filename.substr(0, 3);
      filename.erase(0, 3);

      if ( !fs::exists(directory_stream.str()) )
        fs::create_directory( directory_stream.str() );

      directory_stream << '/';
    }

    filename = directory_stream.str() + filename_stream.str();

    // transaction_id = -1 returns the latest tile available
    std::string output_filename = platefile->read_to_file(filename, col, row, level, transaction_id);
    //    std::cout << "\t--> [ " << col << " " << row << " " << level << "] : Writing "
    //              << output_filename << "\n";

  } catch (TileNotFoundErr &e) {
    //    std::cout << "\t--> [ " << col << " " << row << " " << level << "] : Missing tile\n";
  }
}

// ----------------------------------------------------------------------------
//                                do_level()
// ----------------------------------------------------------------------------

void do_level(int level, BBox2i tile_region, boost::shared_ptr<PlateFile> platefile,
              std::string output_name, std::string output_format, int transaction_id) {

  //  std::cout << "\t--> Exporting tiles -- " << tile_region << " @ level " << level << ".\n";

  // Subdivide the bbox into smaller workunits if necessary.  This
  // helps to keep operations efficient.
  std::list<BBox2i> tile_workunits = bbox_tiles(tile_region, 1024, 1024);
  for ( std::list<BBox2i>::iterator region_iter = tile_workunits.begin();
        region_iter != tile_workunits.end(); ++region_iter) {

    // Fetch the list of valid tiles in this particular workunit.
    std::list<TileHeader> tile_records = platefile->search_by_region(level, *region_iter,
                                                                     transaction_id,
                                                                     transaction_id, 1);

    if (!tile_records.empty()) {
      vw_out() << "\t--> Exporting " << tile_records.size() << " tiles in " << *region_iter
               << " @ level " << level << "\n";
    }

    for ( std::list<TileHeader>::iterator header_iter = tile_records.begin();
          header_iter != tile_records.end(); ++header_iter) {
      if (output_format == "toast") {
        save_toast_tile(output_name, platefile,
                        header_iter->col(), header_iter->row(),
                        header_iter->level(), header_iter->transaction_id());
      } else if (output_format == "gigapan") {
        save_gigapan_tile(output_name, platefile,
                          header_iter->col(), header_iter->row(),
                          header_iter->level(), header_iter->transaction_id());
      } else {
        vw_out() << "Error -- unknown output format: " << output_format << "\n";
        exit(1);
      }

    }
  }
}

// ----------------------------------------------------------------------------
//                                do_all_levels()
// ----------------------------------------------------------------------------

void do_all_levels(std::string platefile_name, std::string output_name,
                   std::string output_format, int transaction_id) {

  // Open the plate file
  boost::shared_ptr<PlateFile> platefile(new PlateFile(platefile_name));
  std::cout << "Exporting " << platefile->num_levels() << " levels of tiles to " << output_name << "\n";

  // Create the output directory
  if ( !fs::exists(output_name) )
    fs::create_directory(output_name);

  // Iterate over the levels
  for (int32 level = 0; level < platefile->num_levels(); ++level) {

    // The entire region contains 2^level tiles.
    int32 region_size = 1 << level;
    int32 subdivided_region_size = region_size / 16;
    if (subdivided_region_size < 1024) subdivided_region_size = 1024;
    BBox2i full_region(0,0,region_size,region_size);
    std::list<BBox2i> workunits = bbox_tiles(full_region,
                                             subdivided_region_size,
                                             subdivided_region_size);
    for ( std::list<BBox2i>::iterator region_iter = workunits.begin();
          region_iter != workunits.end(); ++region_iter) {
      do_level(level, *region_iter, platefile, output_name, output_format, transaction_id);
    }
  }
}

// ----------------------------------------------------------------------------
//                                    main()
// ----------------------------------------------------------------------------

int main( int argc, char *argv[] ) {

  std::string output_file_name;
  std::string output_format;
  std::string plate_file_name;
  int transaction_id;

  po::options_description general_options("Turns georeferenced image(s) into a TOAST quadtree.\n\nGeneral Options");
  general_options.add_options()
    ("output-format,f", po::value<std::string>(&output_format)->default_value("toast"),
     "Output tree format, one of [ toast, toast_dem, gigapan ] ")
    ("output-name,o", po::value<std::string>(&output_file_name),
     "Specify the base output directory")
    ("transaction-id,t", po::value<int>(&transaction_id)->default_value(-1),
     "Specify the transaction id to save.")
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

  // Spider across the platefile, saving all valid tiles.
  do_all_levels(plate_file_name, output_file_name, output_format, transaction_id);
}
