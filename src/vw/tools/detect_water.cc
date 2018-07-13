// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

/**
  The central program for using the water detection tools added to Vision Workbench.
*/

#include <vw/Image/ImageIO.h>
#include <vw/tools/modis_utilities.h>
#include <vw/tools/modis_water_detection.h>
#include <vw/tools/radar.h>
#include <vw/tools/landsat.h>
#include <vw/tools/multispectral.h>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;


int main(int argc, char **argv) {

  std::vector<std::string> input_file_names;
  std::string output_path, mode, dem_path;
  int num_threads    = 0;
  int tile_size      = 512;
  double sensitivity = 1.0;
  cartography::GdalWriteOptions write_options;

  // Parse the command line options.
  
  po::options_description general_options("Runs VW's water detection tools.\n\nGeneral Options");
  general_options.add_options()
    ("output-path,o",    po::value<std::string>(&output_path), "The output file path")
    ("dem-path,d",       po::value<std::string>(&dem_path), "Path to a DEM file to be used with processing.")
    ("num-threads",      po::value<int>(&num_threads)->default_value(0), 
                         "Number of threads to use for writing")
    ("tile-size",        po::value<int>(&tile_size)->default_value(512), 
                         "Tile size used for parallel processing")
    ("sensitivity",      po::value<double>(&sensitivity)->default_value(1.0), 
                         "Lower this to make the algorithm detect more water, increase for less water.")
    ("mode,m",           po::value<std::string>(&mode), 
        "The processing mode. Required.  Options: [sentinel1, landsat, worldview]")
    ("debug",  "Record debugging information.")
    ("help,h", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::vector<std::string> >(&input_file_names));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  std::ostringstream usage;
  usage << "Usage: detect_water [options] <image files>..." <<std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 0;
  }
  
  if (!vm.count("mode")) {
    std::cout << "Error: The 'mode' option must be specified.\n";
    return 0;
  }
  bool debug = vm.count("debug");

  // Check the input mode
  // - SPOT mode is available as a hidden option since there is an 
  //   implementation but it does not work well.
  boost::algorithm::to_lower(mode);
  bool radar_mode     = (mode.find("sentinel1") != std::string::npos);
  bool landsat_mode   = (mode.find("landsat"  ) != std::string::npos);
  bool worldview_mode = (mode.find("worldview") != std::string::npos);
  bool spot_mode      = (mode.find("spot"     ) != std::string::npos);
  if (!radar_mode && !landsat_mode && !worldview_mode && !spot_mode) {
    std::cout << "Error: Unrecognized 'mode' option!  Choices are " 
              << "'sentinel1' , 'landsat', and 'worldview'.\n";
    return -1;
  }
  
  if ((!dem_path.empty()) && (!radar_mode))
    std::cout << "Warning: DEM file is only used when mode is 'sentinel1'\n";
  
  // If a DEM file was provided, make sure it actually exists before we start doing any processing.
  if ((!dem_path.empty()) && !boost::filesystem::exists(dem_path)) {
    std::cout << "Error: Provided DEM file not found!\n";
    return -2;
  }
  
  // Load settings into VW internals.
  write_options.raster_tile_size = Vector2i(tile_size, tile_size);
  vw_settings().set_default_tile_size(tile_size);
  if (num_threads > 0) {
    write_options.num_threads = num_threads;
    vw_settings().set_default_num_threads(write_options.num_threads);
  }

  if( vm.count("output-path") < 1 ) {
    std::cerr << "Error: must specify the output file!" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  // Execute the requested flood detection mode.
  
  if (radar_mode) {
    std::cout << "Processing sentinel-1 image!\n";
    radar::sar_martinis(input_file_names[0], output_path, write_options, dem_path, debug, 
                        tile_size, sensitivity);
    return 0;
  }
  
  if (landsat_mode) {
    std::cout << "Processing Landsat image!\n";
    landsat::detect_water(input_file_names, output_path, write_options, sensitivity, debug);
    return 0;
  }

  if (worldview_mode) {
    std::cout << "Processing WorldView image!\n";
    multispectral::detect_water_worldview23(input_file_names, output_path, write_options, sensitivity, debug);
    return 0;
  }

  if (spot_mode) {
    std::cout << "Processing SPOT image!\n";
    multispectral::detect_water_spot67(input_file_names, output_path, write_options, sensitivity, debug);
    return 0;
  }
  
  // Basic MODIS flood detection was written but never fully developed.
}
