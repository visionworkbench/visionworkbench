// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Image/ImageView.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Camera/BayerFilter.h>

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
  } catch ( po::error &e ) {
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
  catch( Exception& e ) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}
