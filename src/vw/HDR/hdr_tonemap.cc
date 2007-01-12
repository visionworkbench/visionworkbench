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

#include <vw/HDR/CameraCurve.h>
#include <vw/HDR/GlobalToneMap.h>
#include <vw/HDR/LocalToneMap.h>

#include <iostream>
#include <vector>

using namespace std;
using namespace vw;
using namespace vw::hdr;

const double DRAGO_DEFAULT_BIAS = 0.85;
const double DRAGO_DEFAULT_EXPOSURE_FACTOR = 1.0;
const double DRAGO_DEFAULT_LDMAX = 100;

int main( int argc, char *argv[] ) {
  try {
    std::string input_filename, output_filename, curves_input_file;
    float bias, exposure_factor, Ld_max, gamma;
    
    po::options_description desc("Options");
    desc.add_options()
      ("help", "Display this help message")
      ("input-filename", po::value<std::string>(&input_filename), "Specify the input filename.")
      ("output-filename", po::value<std::string>(&output_filename), "Specify the output filename.")
      ("bias,b", po::value<float>(&bias)->default_value(DRAGO_DEFAULT_BIAS), "Drago Tonemapping Parameter: Bias")
      ("exposure-factor,e", po::value<float>(&exposure_factor)->default_value(DRAGO_DEFAULT_EXPOSURE_FACTOR), "Drago Tonemapping Parameter: Exposure Factor")
      ("ld-max,l", po::value<float>(&Ld_max)->default_value(DRAGO_DEFAULT_LDMAX), "Drago Tonemapping Parameter: Ld_max")
      ("use-curves,c", po::value<std::string>(&curves_input_file), "Apply the inverse of a curves file on disk to the tonemapped image.")
      ("gamma,g", po::value<float>(&gamma), "Apply a gamma correction to the tonemapped image.")
      ("verbose,v", "Print out status messages.");
    
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
      std::cout << "Usage: " << argv[0] << " <input filename> <output filename>\n";
      std::cout << desc << std::endl;
      return 1;
    }
    
    if( vm.count("output-filename") != 1 ) {
      std::cout << "Usage: " << argv[0] << " <input filename> <output filename>\n";
      std::cout << desc << std::endl;
      return 1;
    }
    
    bool verbose = true;
    if ( vm.count("verbose") == 0 ) {
      verbose = false;
    }
    
    // Read in the HDR Image
    ImageView<PixelRGB<float> > hdr;
    read_image(hdr, input_filename);
    
    // Apply Drago tone-mapping operator.
    if (verbose) { std::cout << "Applying the Drago tone mapping operator.\n"; }

    ImageView<PixelRGB<float> > tone_mapped = drago_tone_map(hdr, bias, exposure_factor, Ld_max);
    
    // Apply gamma correction.
    if (vm.count("gamma") != 0) {
      if (verbose) { std::cout << "Applying a gamma correction of " << gamma << ".\n"; }
      tone_mapped = pow(copy(tone_mapped), 1/gamma);
    }
    
    // Re-apply camera response curves.
    // First must invert curves calculated earlier.
    if (vm.count("use-curves") != 0) {
      if (verbose) { std::cout << "Applying camera curves from \"" << curves_input_file << "\"\n"; }

      vector<Vector<double> > curves;
      read_curves_file(curves, curves_input_file);
      if (curves.size() != tone_mapped.channels()) {
        cout << "Error: The number of curves does not match the number of channels in the tone mapped image.  Exiting. \n\n";
        exit(1);
      }
      
      // Invert the curves
      vector<Vector<double> > inverse_curves(curves.size());
      for ( unsigned i = 0; i < curves.size(); ++i ) {
        invert_curve(curves[i], inverse_curves[i], VW_HDR_RESPONSE_POLYNOMIAL_ORDER);
      }
      
      // Apply the curves
      psi(tone_mapped, curves);
    }
    
    // Write out the result to disk.
    write_image(output_filename, tone_mapped);    
    
  } catch (vw::Exception &e) {
    std::cout << argv[0] << ": a Vision Workbench error occurred: \n\t" << e.what() << "\nExiting.\n\n";  
    return 1;
  }
}
