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

/// \file tonemap.cc
/// 
/// This is a simple application that reads in an high dynamic range
/// image on disk and tone maps it using one of the operators included
/// in this library.  

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
#include <vw/Image/ImageViewRef.h>

#include <vw/HDR/CameraCurve.h>
#include <vw/HDR/PolynomialCameraCurve.h>
#include <vw/HDR/GlobalToneMap.h>
#include <vw/HDR/LocalToneMap.h>

#include <iostream>
#include <vector>

using namespace std;
using namespace vw;
using namespace vw::hdr;

int main( int argc, char *argv[] ) {
  try {
    std::string input_filename, output_filename, curves_input_file;
    float bias, gamma;
    
    po::options_description desc("Options");
    desc.add_options()
      ("help", "Display this help message")
      ("input-filename", po::value<std::string>(&input_filename), "Specify the input filename.")
      ("output-filename,o", po::value<std::string>(&output_filename)->default_value("tonemapped.png"), "Specify the output filename.")
      ("bias,b", po::value<float>(&bias)->default_value(vw::hdr::DRAGO_DEFAULT_BIAS), "Drago Tonemapping Parameter.  (The default of 0.85 works well for most images)")
      ("gamma,g", po::value<float>(&gamma)->default_value(2.2), "Apply a gamma correction to the tonemapped image.");

    
    po::positional_options_description p;
    p.add("input-filename", 1);
    p.add("output-filename", 2);
    
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
    po::notify( vm );
    
    if( vm.count("help") ) {
      std::cout << desc << std::endl;
      return 1;
    }
    
    if( vm.count("input-filename") != 1 ) {
      std::cout << "Usage: " << argv[0] << " <input filename>\n";
      std::cout << desc << std::endl;
      return 1;
    }
    
    
    // Read in the HDR Image
    DiskImageView<PixelRGB<float> > hdr(input_filename);
    
    // Apply Drago tone-mapping operator.
    ImageViewRef<PixelRGB<double> > tone_mapped = drago_tone_map(hdr, bias);
    
    // Apply gamma correction.
    tone_mapped = pow(tone_mapped, gamma);
    
    // Write out the result to disk.
    write_image(output_filename, tone_mapped);    
    
  } catch (vw::Exception &e) {
    std::cout << argv[0] << ": a Vision Workbench error occurred: \n\t" << e.what() << "\nExiting.\n\n";  
    return 1;
  }
}
