// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

using std::cout;
using std::endl;
using std::vector;
using std::string;

using namespace vw;
using namespace vw::hdr;

int main( int argc, char *argv[] ) {
  try {
    vector<string> input_filenames;
    string output_filename, curve_file;
    float exposure_ratio = 0;

    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "Display this help message")
      ("input-filenames,i", po::value<vector<string> >(&input_filenames), "Specify the input files")
      ("output-filename,o", po::value<string>(&output_filename)->default_value("merged_hdr_image.exr"), "Specify the output filename")
      ("exposure-ratio,e", po::value<float>(&exposure_ratio), "Manually specified exposure ratio for the images (e.g. for increasing power of 2 expoures, you would use a exposure-ratio of 2.0).")
      ("save-curves,c", po::value<string>(&curve_file), "Write the curve lookup tables to a file on disk.")
      ("use-curves,u", po::value<string>(&curve_file), "Read the curve lookup tables to a file on disk.  These curves will be used instead of computing new curves.");

    po::positional_options_description p;
    p.add("input-filenames", -1);

    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).options(desc).positional(p).run(), vm );
    po::notify( vm );

    if( vm.count("help") ) {
      cout << desc << endl;
      return 1;
    }

    if( vm.count("input-filenames") != 1 ) {
      cout << "Usage: " << argv[0] << " <input-filenames> [args]\n";
      cout << desc << endl;
      return 1;
    }

    // Process HDR stack using Exif tags.
    ImageView<PixelRGB<float> > hdr;

    cout << "Getting Brightness Values" << endl;
    // In the absense of EXIF data or if the user has provided an
    // explicit exposure ratio, we go with that value here.
       vector<double> brightness_values;
    if( vm.count("exposure-ratio") != 0 )
      brightness_values = brightness_values_from_exposure_ratio(exposure_ratio, input_filenames.size());
    else
      brightness_values = brightness_values_from_exif(input_filenames);

    // For debugging: (print out the brightness values
    //     for (int i = 0; i < brightness_values.size(); ++i)
    //       vw_out() << "BV: " << brightness_values[i] << "\n";

    vector<ImageViewRef<PixelRGB<float> > > images(input_filenames.size());
    for ( unsigned i=0; i < input_filenames.size(); ++i )
      images[i] = DiskImageView<PixelRGB<float> >(input_filenames[i]);

    cout << "Getting Camera Curves" << endl;
    // Compute the camera curves
    CameraCurveFn curves;
    if ( vm.count("use-curves") != 0 )
      curves = read_curves(curve_file);
    else {
      curves = camera_curves(images, brightness_values);

      cout << "Saving Camera Curves" << endl;
      // Write out the curves to disk as a tabulated file
      if ( vm.count("save-curves") != 0 ) {
        write_curves(curve_file, curves);
      }
    }

    TerminalProgressCallback tpc( "tools.hdr_merge", "Processing");
    // Create the HDR images and write the results to the file
    write_image(output_filename, HighDynamicRangeView<PixelRGB<float> > (images, curves, brightness_values), tpc);

  } catch (vw::Exception &e) {
    vw_out() << argv[0] << ": a Vision Workbench error occurred: \n\t"
             << e.what() << "\nExiting.\n\n";
    return 1;
  }
  return 0;
}
