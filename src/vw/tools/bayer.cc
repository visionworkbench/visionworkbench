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

#include <vw/Core/Exception.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Camera/BayerFilter.h>

#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;

int main( int argc, char *argv[] ) {

  std::string input_file_name, output_file_name;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Display this help message")
    ("input-file", po::value<std::string>(&input_file_name), "Explicitly specify the input file")
    ("output-file,o", po::value<std::string>(&output_file_name)->default_value("output.png"), "Specify the output file");
  po::positional_options_description p;
  p.add("input-file", 1);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << desc;
    return 1;
  }

  if( vm.count("help") ) {
    std::cout << desc << std::endl;
    return 1;
  }

  if( vm.count("input-file") != 1 ) {
    std::cout << "Error: Must specify exactly one input file!" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  try {
    ImageView<PixelGray<float> > image;
    read_image( image, input_file_name );

    ImageView<PixelRGB<float> > result = camera::inverse_bayer_filter(image);
    write_image( output_file_name, result );
  }
  catch (const Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}
