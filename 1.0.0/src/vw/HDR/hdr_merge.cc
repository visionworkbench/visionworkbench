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

#include <vw/Image/ImageView.h>
#include <vw/FileIO.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageMath.h>
#include <vw/HDR/LDRtoHDR.h>
#include <vw/HDR/GlobalToneMap.h>

#include <vector>

#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;
using namespace vw;
using namespace vw::hdr;

int main( int argc, char *argv[] ) {
  try {
    std::vector<std::string> filenames;
    std::string output_filename, generated_curves_file, curves_input_file;
    float ratio = 1.414; // default ratio is a factor of two in exposure
    
    po::options_description desc("Options");
    desc.add_options()
      ("help", "Display this help message")
      ("filenames", po::value<std::vector<std::string> >(&filenames), "Specify the input files")
      ("output-filename,o", po::value<std::string>(&output_filename), "Specify the output filename")
      ("exposure-ratio,e", po::value<float>(&ratio), "Manually specified exposure ratio for the images (in units of f-stops).")
      ("generate-curves,g", po::value<std::string>(&generated_curves_file), "Write the curves file to disk using this filename.")
      ("use-curves,c", po::value<std::string>(&curves_input_file), "Read the curves file from disk using this filename.")
      ("verbose,v", "Print out status messages.");

    po::positional_options_description p;
    p.add("filenames", -1);
    
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
    po::notify( vm );
    
    if( vm.count("help") ) {
      std::cout << desc << std::endl;
      return 1;
    }

    if( vm.count("filenames") != 1 ) {
      std::cout << "Usage: " << argv[0] << " <filenames>\n";
      std::cout << desc << std::endl;
      return 1;
    }

    bool verbose = true;
    if ( vm.count("verbose") == 0 ) {
      verbose = false;
    }

    if (vm.count("output-filename") != 1) {
      // Replace this with a more sophisticated algorithm that finds
      // the longest substring shared by all of the input filenames.
      output_filename = "merged-hdr.exr";
    }

    vector<Vector<double> > curves;
    // Process HDR stack using Exif tags.
    if (verbose) { std::cout << "Merging images into a single HDR luminance map.\n"; }
    ImageView<PixelRGB<float> > hdr;

    // User has provided an explicit exposure ratio
    if( vm.count("exposure-ratio") != 0 ) {
      // In the absence of EXIF data, we must rely on the
      // user-provided exposure ratio and read in the images
      // ourselvelves.
      if (verbose) { cout << "One or more files have no EXIF metadata.  Falling back on command line exposure ratio.\n"; }
      if (ratio == 0) {
        cout << "Error: the exposure ratio was not supplied at the command line.  Exiting.\n\n";
        return 1;
      }
      vector<ImageView<PixelRGB<float> > > images(filenames.size());
      for ( unsigned i=0; i < filenames.size(); ++i )
        read_image(images[i], filenames[i]);

      // If the user has supplied a camera curve file, read it in and
      // use it.  Otherwise, ottempt to deduce the camera curves based
      // on the image data itself.
      if ( vm.count("use-curves") != 0 ) 
        read_curves_file(curves, curves_input_file);
      else
        curves = camera_curves(images, ratio);
      
      // Using the images and the camera curves, assemble the HDR image.
      hdr = process_ldr_images(images, curves, ratio);


    // No exposure ratio supplied.  
    } else { 
      try {
        hdr = process_ldr_images_exif(filenames, curves);
      } catch (vw::camera::ExifErr &e) {
        std::cout << "Error: could not parse EXIF data from images and no exposure ratio was supplied.  \nExiting.\n\n";
        exit(1);
      }
    }

    if( vm.count("generate-curves") != 0 ) {
      write_curves_file(generated_curves_file, curves);
    }

    // Write out the result file
    write_image(output_filename, hdr);

  } catch (vw::Exception &e) {
    std::cout << argv[0] << ": a Vision Workbench error occurred: \n\t" << e.what() << "\nExiting.\n\n";  
    return 1;
  }  
  return 0;
}
