// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

/// \file hdr_merge.cc
/// 
/// This is a simple example program that stitches a HDR stack into an
/// HDR image, performs tone-mapping, and saves several versions with
/// different post-processing applied for comparison.  Usually the
/// best image is produced by re-applying the camera response curves
/// and then gamma correcting.
///
/// You can specify the exposure spacing in the stack, however if EXIF
/// information is available the exposure ratios will be inferred
/// automatically.

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/HDR.h>

#include <vector>
#include <string>

using namespace std;
using namespace vw;
using namespace vw::hdr;

int main( int argc, char *argv[] ) {
  try {
    std::vector<std::string> input_filenames;
    std::string output_filename, curve_file;
    float exposure_ratio = 0;
    
    po::options_description desc("Options");
    desc.add_options()
      ("help", "Display this help message")
      ("input-filenames", po::value<std::vector<std::string> >(&input_filenames), "Specify the input files")
      ("output-filename,o", po::value<std::string>(&output_filename)->default_value("hdr.exr"), "Specify the output filename")
      ("exposure-ratio,e", po::value<float>(&exposure_ratio), "Manually specified exposure ratio for the images (in units of f-stops).")
      ("save-curves,c", po::value<std::string>(&curve_file), "Write the curve lookup tables to a file on disk.");

    po::positional_options_description p;
    p.add("input-filenames", -1);
    
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
    po::notify( vm );
    
    if( vm.count("help") ) {
      std::cout << desc << std::endl;
      return 1;
    }

    if( vm.count("input-filenames") != 1 ) {
      std::cout << "Usage: " << argv[0] << " <input-filenames> [args]\n";
      std::cout << desc << std::endl;
      return 1;
    }

    // Process HDR stack using Exif tags.
    ImageView<PixelRGB<float> > hdr;

    // In the absense of EXIF data or if the user has provided an
    // explicit exposure ratio, we go with that value here.
       std::vector<double> brightness_values;
    if( vm.count("exposure-ratio") != 0 ) 
      brightness_values = brightness_values_from_exposure_ratio(exposure_ratio, input_filenames.size());
    else
      brightness_values = brightness_values_from_exif(input_filenames);

    // For debugging: (print out the brightness values
    //     for (int i = 0; i < brightness_values.size(); ++i) 
    //       vw_out(0) << "BV: " << brightness_values[i] << "\n";

    vector<ImageViewRef<PixelRGB<float> > > images(input_filenames.size());
    for ( unsigned i=0; i < input_filenames.size(); ++i )
      images[i] = DiskImageView<PixelRGB<float> >(input_filenames[i]);

    // Compute the camera curves
    CameraCurveFn curves = camera_curves(images, brightness_values);

    // Write out the curves to disk as a tabulated file
    if ( vm.count("save-curves") != 0 ) {
      write_curves(curve_file, curves);
    }

    // Create the HDR images and write the results to the file
    write_image(output_filename, HighDynamicRangeView<PixelRGB<float> > (images, curves, brightness_values) );

  } catch (vw::Exception &e) {
    std::cout << argv[0] << ": a Vision Workbench error occurred: \n\t" << e.what() << "\nExiting.\n\n";  
    return 1;
  }  
  return 0;
}
