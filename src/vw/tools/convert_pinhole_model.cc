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

#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;
using namespace vw::camera;

int main( int argc, char *argv[] ) {

  std::string image_file_name, output_file_name, output_model_type, camera_file_name;
  int sample_spacing, rpc_degree;
  
  po::options_description desc("Usage: convert_pinhole_model [options] <input image> <camera model> \n\nOptions");
  desc.add_options()
    ("help,h",        "Display this help message.")
    ("input-file",    po::value<std::string>(&image_file_name), 
                      "Explicitly specify the input file.")
    ("camera-file",    po::value<std::string>(&camera_file_name), 
                      "Explicitly specify the camera file.")
    ("sample-spacing", po::value(&sample_spacing)->default_value(0),    
                       "Pick one out of this many consecutive pixels to sample. if not specified, it will be auto-computed.")
    ("output-type", po::value<std::string>(&output_model_type)->default_value("TsaiLensDistortion"), 
     "The output model type. Options: TsaiLensDistortion, BrownConradyDistortion, RPC.")
    ("rpc-degree", po::value(&rpc_degree)->default_value(3),    
                       "The degree of the polynomials, if the output distortion model is RPC.")
    ("output-file,o", po::value<std::string>(&output_file_name)->default_value("output.tsai"), 
     "Specify the output file. It is expected to have the .tsai extension.");
  
  po::positional_options_description p;
  p.add("input-file",  1);
  p.add("camera-file", 1);

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
  if( (vm.count("input-file") != 1) || (vm.count("camera-file") != 1) ) {
    std::cout << "Error: Must specify exactly one image file and one camera file!" << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  if (rpc_degree <= 0) {
    std::cout << "Error: The RPC degree must be positive." << std::endl;
    std::cout << desc << std::endl;

  }

  try {
    // Get the size of the input image
    boost::shared_ptr<vw::DiskImageResource> image_in(vw::DiskImageResource::open(image_file_name));
    Vector2i image_size(image_in->format().cols, image_in->format().rows);
    
    if (sample_spacing <= 0) {
      sample_spacing = auto_compute_sample_spacing(image_size);
      vw_out() << "Sample the image by picking one in every " << sample_spacing << " pixels.\n";
    }

    vw_out() << "Loading camera model file: " << camera_file_name.c_str() << "\n";

    // Here we will accept an optical bar model too, then in_model
    // will be a pointer to either that or to pinhole.
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
        vw_out() << "Read an OpticalBarModel camera model.\n";
      }catch(...){}
    }
    
    if (!success) 
      vw_throw(ArgumentErr() << "Could not read the camera model.\n");
    
    PinholeModel out_model;
    bool force_conversion = true;
    
    //double error;
    if (output_model_type == "TsaiLensDistortion") {
      //error =
      create_approx_pinhole_model<TsaiLensDistortion>
	(in_model, out_model, image_size, sample_spacing, force_conversion, rpc_degree);
    }else if (output_model_type == "BrownConradyDistortion") {
      //error =
      create_approx_pinhole_model<BrownConradyDistortion>
        (in_model, out_model, image_size, sample_spacing, force_conversion, rpc_degree);
    } else if (output_model_type == RPCLensDistortion::class_name()) {
      //error =
      create_approx_pinhole_model<RPCLensDistortion>
	(in_model, out_model, image_size, sample_spacing, force_conversion, rpc_degree);
    }else{
      vw_out() << "Unsupported output model type: " << output_model_type << "\n";
      return 1;
    }
    
    //std::cout << "Approximation error = " << error << std::endl;
    printf("Writing output model: %s\n", output_file_name.c_str());
    out_model.write(output_file_name);
  }
  catch (const Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  printf("Finished!\n");

  return 0;
}
