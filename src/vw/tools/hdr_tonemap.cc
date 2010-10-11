// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Filter.h>
#include <vw/FileIO.h>

#include <vw/HDR/CameraCurve.h>
#include <vw/HDR/GlobalToneMap.h>
#include <vw/HDR/LocalToneMap.h>

#include <iostream>

using std::cout;
using std::endl;
using std::string;
using namespace vw;
using namespace vw::hdr;

int main( int argc, char *argv[] ) {
  try {
    string input_filename, output_filename, curve_file;
    float bias, gamma;
    int bit_depth;

    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "Display this help message")
      ("input-filename", po::value<string>(&input_filename), "Specify the input filename.")
      ("output-filename,o", po::value<string>(&output_filename)->default_value("tonemap.png"), "Specify the output filename.")
      ("bias,b", po::value<float>(&bias)->default_value(vw::hdr::DRAGO_DEFAULT_BIAS), "Drago Tonemapping Parameter.  (The default of 0.85 works well for most images)")
      ("gamma,g", po::value<float>(&gamma)->default_value(2.2), "Apply a gamma correction to the tonemapped image.")
      ("bit-depth", po::value<int>(&bit_depth)->default_value(8), "Select a bit depth (8, 16, 32, or 64  [selecting 32 or 64 bit saves as IEE float if supported]")
      ("curves-for-raw-image,c", po::value<string>(&curve_file), "Read the curve lookup tables to a file on disk.  Only use this option if you are processing a raw image from a camera and you have a curves file to match.  You need not specify this option if you are processing a HDR luminance image.");

    po::positional_options_description p;
    p.add("input-filename", 1);
    p.add("output-filename", 2);

    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
    po::notify( vm );

    if( vm.count("help") ) {
      cout << desc << endl;
      return 1;
    }

    if( vm.count("input-filename") != 1 ) {
      cout << "Usage: " << argv[0] << " <input filename>\n";
      cout << desc << endl;
      return 1;
    }


    // Read in the HDR Image
    ImageViewRef<PixelRGB<double> > tone_mapped = DiskImageView<PixelRGB<double> >(input_filename);

    // If the user has specified curves *and* is processing a raw
    // camera image, we apply them to the image here *before* applying
    // the tane mapping operator, since the tone mapping operator is
    // designed to operate on luminance values (though scale here is not important).
    if ( vm.count("curves-for-raw-image") != 0 ) {
      cout << "Reading Camera Curves" << endl;
      CameraCurveFn curves = read_curves(curve_file);
      tone_mapped = luminance_image(tone_mapped, curves, 1.0);
    }

    cout << "Applying Drago tone-mapping" << endl;
    // Apply Drago tone-mapping operator.
    tone_mapped = pow(drago_tone_map(tone_mapped, bias), gamma);

    TerminalProgressCallback tpc( "tools.hdr_tonemap", "Processing");

    // Write out the result to disk.
    if (bit_depth == 8)
      write_image(output_filename, channel_cast_rescale<uint8>(clamp(tone_mapped)), tpc);
    else if (bit_depth == 16)
      write_image(output_filename, channel_cast_rescale<uint16>(clamp(tone_mapped)), tpc);
    else if (bit_depth == 32)
      write_image(output_filename, channel_cast_rescale<float32>(clamp(tone_mapped)), tpc);
    else if (bit_depth == 64)
      write_image(output_filename, tone_mapped, tpc);
    else {
      vw_out() << "Unknown bit depth specified by user.  Please choose from the following options: [ 8, 16, 32, 64 ]\n";
      exit(1);
    }

  } catch (vw::Exception &e) {
    vw_out() << argv[0] << ": a Vision Workbench error occurred: \n\t" << e.what() << "\nExiting.\n\n";
    return 1;
  }
}
