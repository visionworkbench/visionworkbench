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
  Program to exercise VW's water detection tools.
*/

#include <vw/Image/ImageIO.h>
#include <vw/tools/modis_utilities.h>
#include <vw/tools/modis_water_detection.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;


int main(int argc, char **argv) {

  std::vector<std::string> input_file_names;
  std::string output_path;
  float manual_threshold;


  po::options_description general_options("Runs VW's water detection tools.\n\nGeneral Options");
  general_options.add_options()
    ("output-path,o",    po::value<std::string>(&output_path    ), "The output file path")
    ("manual-threshold", po::value<float>(&manual_threshold)->default_value(0), "Manually specify a threshold to use")
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

  // Handle user threshold
  float threshold = 100;
  if( vm.count("manual-threshold") )
    threshold = manual_threshold;

  if( vm.count("output-path") < 1 ) {
    std::cerr << "Error: must specify the output file!" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  std::cout << "Loading MODIS image...\n";
  modis::ModisImage modis_image;
  modis::load_modis_image(modis_image, input_file_names);
/*
  write_image("b1.tif", select_channel(modis_image, modis::B1));
  write_image("b2.tif", select_channel(modis_image, modis::B2));
  write_image("b3.tif", select_channel(modis_image, modis::B3));
  write_image("b4.tif", select_channel(modis_image, modis::B4));
*/
  cartography::GeoReference georef;
  modis::load_modis_georef(input_file_names, georef);
  std::cout << "Read georef: \n" << georef << std::endl;
  
  std::cout << "Forming MODIS products...\n";
  modis::ModisProductImage product_image;
  modis::form_modis_product_image(modis_image, product_image);
  
  // DEBUG
  write_image("ndvi.tif", select_channel(product_image, modis::NDVI));
  write_image("ndwi.tif", select_channel(product_image, modis::NDWI));
  write_image("evi.tif" , select_channel(product_image, modis::EVI ));
  write_image("lswi.tif", select_channel(product_image, modis::LSWI));
  

  std::cout << "Classifying MODIS image...\n";
  modis::WaterDetection water_result;
  modis::for_each_pixel(modis_image, product_image, water_result, modis::detect_water_evi_functor());
  write_image("evi.tif", water_result);
  modis::for_each_pixel(modis_image, product_image, water_result, modis::detect_water_xiao_functor());
  write_image("xiao.tif", water_result);
  modis::for_each_pixel(modis_image, product_image, water_result, modis::detect_water_diff_functor(threshold));
  write_image("diff.tif", water_result);
  modis::for_each_pixel(modis_image, product_image, water_result, modis::detect_water_dartmouth_functor(threshold));
  write_image("dartmouth.tif", water_result);
  modis::for_each_pixel(modis_image, product_image, water_result, modis::detect_water_mod_ndwi_functor(threshold));
  write_image("ndwi.tif", water_result);
  modis::for_each_pixel(modis_image, product_image, water_result, modis::detect_water_fai_functor(threshold));
  write_image("fai.tif", water_result);
  
  //write_image(output_path, water_result);
  
  
  
  std::cout << "Finished!\n";


  return 0;
}
