// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file vwv.cc
///
/// The Vision Workbench image viewer.
///

// Qt
#include <QApplication>
#include <QWidget>

// Boost
#include <boost/program_options.hpp>
using namespace boost;
namespace po = boost::program_options;

// VW
#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/Math.h>
#include <vw/FileIO.h>
#include <vw/Stereo.h>
using namespace vw;

#include <vw/gui/MainWindow.h>

// Allows FileIO to correctly read/write unusual pixel types
namespace vw {
  template<> struct PixelFormatID<Vector3>   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_3_CHANNEL; };
}

void print_usage(po::options_description const& visible_options) {
  vw_out() << "\nUsage: vwv [options] <image file> \n";
  vw_out() << visible_options << std::endl;
}

int main(int argc, char *argv[]) {

  unsigned cache_size;
  std::string image_filename;
  float nodata_value;

  // Set up command line options
  po::options_description visible_options("Options");
  visible_options.add_options()
    ("help,h", "Display this help message")
    ("normalize,n", "Attempt to normalize the image before display.")
    ("nodata-value", po::value<float>(&nodata_value), "Choose a \"nodata\" value in the image to treat as transparent.")
    ("cache", po::value<unsigned>(&cache_size)->default_value(1000), "Cache size, in megabytes");

  po::options_description positional_options("Positional Options");
  positional_options.add_options()
    ("image", po::value<std::string>(&image_filename), "Input Image Filename");
  po::positional_options_description positional_options_desc;
  positional_options_desc.add("image", 1);

  po::options_description all_options;
  all_options.add(visible_options).add(positional_options);

  // Parse and store options
  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_options_desc).run(), vm );
  po::notify( vm );

  // If the command line wasn't properly formed or the user requested
  // help, we print an usage message.
  if( vm.count("help") ) {
    print_usage(visible_options);
    exit(0);
  }

  if( vm.count("image") != 1 ) {
    vw_out() << "Error: Must specify exactly one input file!" << std::endl;
    print_usage(visible_options);
    return 1;
  }

  // Set the Vision Workbench cache size
  vw_settings().set_system_cache_size( cache_size*1024*1024 );

  // Start up the Qt GUI
  QApplication app(argc, argv);
  vw::gui::MainWindow main_window(image_filename, nodata_value, vm.count("normalize"), vm);
  main_window.show();
  try {
    app.exec();
  } catch (vw::Exception &e) {
    vw_out() << "An unexpected error occurred: " << e.what() << "\nExiting\n\n";
  }

  return 0;
}

