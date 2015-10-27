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


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Cartography/GeoReference.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>

using namespace vw;

/// \file georef_query.cc Tool for making queries about GeoRef information on an image file.


int main( int argc, char *argv[] )
{
  // Output formats list
  const int OUTPUT_LAT_LON   = 0;
  const int OUTPUT_PROJECTED = 1;
  const int OUTPUT_PIXELS    = 2;

  std::string inputImagePath;
  bool        printBounds;
  double      lat, lon, row, col;
  int         outputFormat;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h", "Display this help message")
    ("print-bounds", po::bool_switch(&printBounds)->default_value(false),    "Prints the lat/lon bounds")
    ("row",           po::value<double>(&row),  "Enter the row    of a pixel")
    ("col",           po::value<double>(&col),  "Enter the column of a pixel")
    ("lat",           po::value<double>(&lat),  "Enter a latitude  value")
    ("lon",           po::value<double>(&lon),  "Enter a longitude value")
    ("output-format", po::value<int>(&outputFormat)->default_value(OUTPUT_LAT_LON),
                           "Specifies the output format: 0 = Lat/Lon, 1 = Projected, 2 = Pixels");

  po::options_description positional("");
  positional.add_options()
    ("inputImage",   po::value(&inputImagePath), "Path to input geotiff file");

  po::positional_options_description positional_desc;
  positional_desc.add("inputImage", 1);
  
  std::string usage("[options] <input image>\n");
  po::variables_map vm;
  try {
    po::options_description all_options;
    all_options.add(general_options).add(positional);

    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).style( po::command_line_style::unix_style ).run(), vm );

    po::notify( vm );
  } catch (po::error const& e) {
    vw::vw_throw( vw::ArgumentErr() << "Error parsing input:\n"
                  << e.what() << "\n" << usage << general_options );
  }

  // Check for required inputs
  if ( !vm.count("inputImage") )
    vw_throw( vw::ArgumentErr() << "Requires <input image> in order to proceed.\n\n"
              << usage << general_options );
  if ( !vm.count("print-bounds") && !(vm.count("row") && vm.count("col")) && !(vm.count("lat") && vm.count("lon")) )
    vw_throw( vw::ArgumentErr() << "Requires either print-bounds, row/col, or lat/lon in order to proceed.\n\n"
              << usage << general_options );



  // Load the input DEM georeference
  cartography::GeoReference georef; 
  boost::scoped_ptr<DiskImageResource> inputImage(DiskImageResource::open(inputImagePath));
  if (!read_georeference(georef, *inputImage)) {
    //vw_out() << "Failed to read input image!\n";
    std::cout << "Failed to read input image georeference!\n";
    return false;
  }

  // Set output display text
  std::string yString, xString;
  switch (outputFormat) {
    case OUTPUT_LAT_LON:
      yString = "latitude ";
      xString = "longitude";
      break;
    case OUTPUT_PROJECTED:
      yString = "y";
      xString = "x";
      break;
    case OUTPUT_PIXELS:
      yString = "row";
      xString = "col";
      break;
    default: // Error
          vw_throw( vw::ArgumentErr() << "Invalid value given for outputFormat!.\n\n"
                        << usage << general_options );
  }


  if (printBounds) { // Print the bounds the user requested
    // Get the bounding box of the entire image
    vw::BBox2 fullImageBox(Vector2(0,0), Vector2(inputImage->cols(), inputImage->rows()));

    // Get the requested output bounding box
    vw::BBox2 outputBox;
    switch (outputFormat) {
      case OUTPUT_LAT_LON: // Geodetic coordinates
        outputBox = georef.pixel_to_lonlat_bbox(fullImageBox);
        break;
      case OUTPUT_PROJECTED:
        outputBox = georef.pixel_to_point_bbox(fullImageBox);
        break;
      case OUTPUT_PIXELS:
        outputBox = fullImageBox; // Not much to do here
        break;
    }
    std::cout <<   "Min " << yString << " = " << outputBox.min()[1]
              << "\nMax " << yString << " = " << outputBox.max()[1]
              << "\nMin " << xString << " = " << outputBox.min()[0]
              << "\nMax " << xString << " = " << outputBox.max()[0] << "\n";
    return 0;
  } // end printBounds case

  // Set up the input coordinate
  vw::Vector2 inputLoc, outputLoc;
  if (vm.count("row") && vm.count("col")) {
    inputLoc[0] = col;
    inputLoc[1] = row;

    switch (outputFormat){
      case OUTPUT_LAT_LON: // Geodetic coordinates
        outputLoc = georef.pixel_to_lonlat(inputLoc);
        break;
      case OUTPUT_PROJECTED:
        outputLoc = georef.pixel_to_point(inputLoc);
        break;
      case OUTPUT_PIXELS:
        outputLoc = inputLoc; // Not much to do here
        break;
    }
  } // End pixel input case

  if (vm.count("lat") && vm.count("lon")) {
    inputLoc[0] = lon;
    inputLoc[1] = lat;

    switch (outputFormat){
      case OUTPUT_LAT_LON: // Geodetic coordinates
        outputLoc = inputLoc; // Not much to do here
        break;
      case OUTPUT_PROJECTED:
        outputLoc = georef.lonlat_to_point(inputLoc);
        break;
      case OUTPUT_PIXELS:
        outputLoc = georef.lonlat_to_pixel(inputLoc);
        break;
    }
  } // End lat/lon input case
  std::cout << xString << " = " << outputLoc[0] << "\n"
            << yString << " = " << outputLoc[1] << "\n";


  return 0;
}








