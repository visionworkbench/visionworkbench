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

/// Tool to generate a new pinhole camera model that approximates an existing model.
/// - The only current usage is to swap from one model with slow point-to-pixel performance
///   to a TsaiLensDistortion or BrownConradyDistortion with fast point-to-pixel performance.


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/Core/Exception.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Interpolation.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/CameraUtilities.h>
#include <vw/FileIO/GdalWriteOptions.h>
#include <vw/FileIO/FileUtils.h>

#include <iostream>

// TODO(oalexan1): Add option --min-rpc-coeff-val. Check if this results
// in RPC distortion optimizing very small coefficients.

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;
using namespace vw::camera;

// TODO(oalexan1): Make all options below use the Options structure
struct Options: vw::GdalWriteOptions {};

int main(int argc, char *argv[]) {

  std::string image_file_name, output_file_name, output_model_type, camera_file_name;
  int sample_spacing, rpc_degree;
  double camera_to_ground_dist;
  Vector2 image_size;
  std::string image_size_str;

  Options opt;
  po::options_description general_options
  ("Usage: convert_pinhole_model [options] <input image> <camera model> ""\n\nOptions");
  general_options.add_options()
    ("help,h",
     "Display this help message.")
    ("input-file", po::value<std::string>(&image_file_name), 
     "Explicitly specify the input file.")
    ("camera-file", po::value<std::string>(&camera_file_name), 
     "Explicitly specify the camera file.")
    ("sample-spacing", po::value(&sample_spacing)->default_value(0),    
     "Pick one out of this many consecutive pixels to sample. if not specified, it will "
     "be auto-computed.")
    ("output-type", po::value<std::string>(&output_model_type)->default_value("TsaiLensDistortion"), 
     "The output model type. Options: TsaiLensDistortion, BrownConradyDistortion, RPC.")
    ("rpc-degree", po::value(&rpc_degree)->default_value(3),
     "The degree of the polynomials, if the output distortion model is RPC.")
    ("camera-to-ground-dist", po::value(&camera_to_ground_dist)->default_value(0),    
     "The distance from the camera to the ground, in meters. This is necessary to convert an optical bar model to pinhole.")
    ("output-file,o", 
     po::value<std::string>(&output_file_name)->default_value("output.tsai"), 
     "Specify the output file. It is expected to have the .tsai extension.")
    ("image-size",
     po::value(&image_size_str)->default_value(""),
     "Image width and height, specified as two numbers in quotes and separated by a "
     "space, unless the input image file is provided.");

  general_options.add(vw::GdalWriteOptionsDescription(opt));
  
  po::positional_options_description p;
  p.add("input-file",  1);
  p.add("camera-file", 1);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(general_options).positional(p).run(), vm);
    po::notify(vm);
  } catch (const po::error& e) {
    vw_out() << "An error occurred while parsing command line arguments.\n";
    vw_out() << "\t" << e.what() << "\n\n";
    vw_out() << general_options;
    return 1;
  }

  if (vm.count("help")) {
    vw_out() << general_options << std::endl;
    return 1;
  }

  // This must happen after the options are parsed
  opt.setVwSettingsFromOpt();
  create_out_dir(output_file_name);

  // Parse the image size from a string
  if (image_size_str != "") {
    std::istringstream iss(image_size_str);
    if (! (iss >> image_size[0] >> image_size[1])) {
      vw_out() << "Could not parse correctly the image size. "
               << "The image dimensions must be in quotes.\n";
      return 1;
    }
  }
  
  // Check if the image dimensions were specified
  if (image_size[0] && image_size[1] > 0) {

    // If the image size is specified, just the camera must be specified
    if (camera_file_name != "") {
      vw_out() << "If the image size is provided, the input image file must not be specified.\n";
      return 1;
    }

    // If only the camera is specified, it will populate the image file field, as that's the
    // first one.
    if (image_file_name == "") {
      vw_out() << "The input camera file is not specified.\n";
      return 1;
    } else {
      camera_file_name = image_file_name;
      image_file_name = "";
    }
    
  } else {
    // Must specify both the image and camera
    if ((vm.count("input-file") != 1) || (vm.count("camera-file") != 1)) {
      vw_out() << "Error: Must specify exactly one image file and one camera file. "
               << "At least provide the image size via the provided option "
               << "if not the image itself." << std::endl;
      vw_out() << general_options << std::endl;
      return 1;
    }
  }
  
  if (rpc_degree <= 0) {
    vw_out() << "Error: The RPC degree must be positive." << std::endl;
    vw_out() << general_options << std::endl;
    return 1;
  }

  try {
    // Get the size of the input image, unless the dimensions were already specified
    if (image_size[0] <= 0 || image_size[1] <= 0) {
      boost::shared_ptr<vw::DiskImageResource> image_in
        (vw::DiskImageResource::open(image_file_name));
      image_size = vw::Vector2(image_in->format().cols, image_in->format().rows);
    }

    vw_out() << "Image width and height: " << image_size[0] << ' ' << image_size[1] << "\n";
    if (sample_spacing <= 0) {
      sample_spacing = auto_compute_sample_spacing(image_size);
      vw_out() << "Sample the image by picking one in every " 
               << sample_spacing << " pixels.\n";
    }
     
    // Here we will accept an optical bar model too, then in_model
    // will be a pointer to either that or to pinhole.
    vw_out() << "Loading camera model file: " << camera_file_name.c_str() << "\n";
    CameraModel * in_model;
    OpticalBarModel opb;
    PinholeModel    pin;
    bool success = false;
    try{
      pin.read(camera_file_name);
      in_model = &pin;
      success = true;
      vw_out() << "Read a Pinhole camera model.\n";
    }catch(...){}

    if (!success) {
      try {
        opb.read(camera_file_name);
        in_model = &opb;
        success = true;
        vw_out() << "Read an OpticalBarModel camera.\n";
        if (camera_to_ground_dist <= 0) {
          vw_out() << "Must set the camera to ground distance "
                   << "if the input is an optical bar model.\n";
          return 1;
        }
      } catch(...){}
    }
    
    if (!success) 
      vw_throw(ArgumentErr() << "Could not read the camera model.\n");
    
    bool force_conversion = true;
    PinholeModel out_model 
      = fitPinholeModel(in_model, image_size, output_model_type,
                        force_conversion,
                        sample_spacing, rpc_degree, camera_to_ground_dist);

    vw_out() << "Writing output model: " << output_file_name.c_str() << "\n";
    out_model.write(output_file_name);
  }
  catch (const Exception& e) {
    vw_out() << "Error: " << e.what() << std::endl;
  }

  return 0;
}
