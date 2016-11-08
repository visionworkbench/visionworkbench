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

/// Tool to undistort a pinhole camera image given the camera model file.


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
#include <vw/Camera/LensDistortion.h>
#include <vw/FileIO/FileUtils.h>
#include <vw/tools/Common.h>

#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace vw;
using vw::camera::PinholeModel;
using vw::camera::LensDistortion;

// Global variables, to make it easier to invoke the function do_work
// with many channels and channel types.
std::string input_file_name, output_file_name, camera_file_name;

/// Generates an undistorted view of an input image.
/// - This should be moved somewhere else.
/// - Could easily add the option to keep the image the same size.
template <class PixelT>
void do_work() {

  // Load the input image
  printf("Loading input image: %s\n", input_file_name.c_str());
  ImageView<PixelT> image_in;
  read_image(image_in, input_file_name);

  // Interpolate the input image
  InterpolationView<EdgeExtensionView<ImageView<PixelT>, ZeroEdgeExtension>, BilinearInterpolation> interp_image_in = interpolate(image_in, BilinearInterpolation(), ZeroEdgeExtension());
  
  // Output image
  ImageView<PixelT> image_out;

  // Load the camera
  printf("Loading camera model file: %s\n", camera_file_name.c_str());
  PinholeModel camera_model(camera_file_name);
  
  printf("Undistorting the image...\n");

  const int width_in  = interp_image_in.cols();
  const int height_in = interp_image_in.rows();
  const LensDistortion* lens_ptr = camera_model.lens_distortion();
  const double pitch = camera_model.pixel_pitch();

  //std::cout << camera_model << std::endl;

  // Figure out the size of the undistorted image
  // - Iterate along each side of the image and record where the output pixels go
  Vector2 out_loc, lens_loc;
  BBox2   output_area;
  for (int r=0; r<height_in; ++r) {
    lens_loc = elem_prod(Vector2(0, r), pitch);
    out_loc  = lens_ptr->undistorted_coordinates(camera_model, lens_loc);
    output_area.grow(elem_quot(out_loc, pitch));
    
    lens_loc = elem_prod(Vector2(width_in-1, r), pitch);
    out_loc  = lens_ptr->undistorted_coordinates(camera_model, lens_loc);
    output_area.grow(elem_quot(out_loc, pitch));
  }
  for (int c=0; c<width_in; ++c) {
    lens_loc = elem_prod(Vector2(c, 0), pitch);
    out_loc  = lens_ptr->undistorted_coordinates(camera_model, lens_loc);
    output_area.grow(elem_quot(out_loc, pitch));
    
    lens_loc = elem_prod(Vector2(c, height_in-1), pitch);
    out_loc  = lens_ptr->undistorted_coordinates(camera_model, lens_loc);
    output_area.grow(elem_quot(out_loc, pitch));
    //std::cout << lens_loc << " --> " << out_loc << " --> " << elem_quot(out_loc, pitch) << std::endl;
  }
  //std::cout << camera_model << std::endl;
  image_out.set_size(floor(output_area.width()), floor(output_area.height()));
  Vector2 offset = output_area.min();

  // Fill in the undistorted image
  for (int r=0; r<image_out.rows(); ++r) {
    for (int c=0; c<image_out.cols(); ++c) {
      // For each output pixel, find the associated input pixel.
      // - Need to add in an offset so that the pixel coords we pass in are in the coordinate
      //   frame of the input image rather than of the output image.
      lens_loc = elem_prod(Vector2(c,r)+offset, pitch);
      //std::cout << lens_loc << std::endl;
      out_loc  = lens_ptr->distorted_coordinates(camera_model, lens_loc);
      Vector2 in_loc = elem_quot(out_loc, pitch);
      image_out(c,r) = interp_image_in(in_loc[0], in_loc[1]);
    }
  } // End double loop through output image

  printf("Writing output image: %s\n", output_file_name.c_str());
  write_image(output_file_name, image_out);
  
  printf("Finished!\n");
  
} // End function do_work

// Magic to make do_work() work with many number of channels and channel types
#define DO_WORK(PIXELTYPEMACRO, CHANNELTYPE) DO_WORK_(PIXELTYPEMACRO, CHANNELTYPE)
#define DO_WORK_(PIXELTYPEMACRO, CHANNELTYPE) DO_WORK__(PIXELTYPEMACRO, CHANNELTYPE, PIXELTYPEMACRO ## _ ## CHANNELTYPE)
#define DO_WORK__(PIXELTYPEMACRO, CHANNELTYPE, FUNCSUFFIX) \
  void do_work_##FUNCSUFFIX(void) {                                    \
    do_work<PIXELTYPEMACRO<CHANNELTYPE > >();                          \
  }

// The channel type can be integer or float
#define DO_WORK_ALL_CHANNELS( PIXELTYPE )  \
  DO_WORK( PIXELTYPE, uint8 );             \
  DO_WORK( PIXELTYPE, int8 );              \
  DO_WORK( PIXELTYPE, uint16 );            \
  DO_WORK( PIXELTYPE, int16 );             \
  DO_WORK( PIXELTYPE, float32 );           \
  DO_WORK( PIXELTYPE, float64 ); 

// Support 1 or 3 channel images, perhaps with an alpha channel
DO_WORK_ALL_CHANNELS(PixelGray)
DO_WORK_ALL_CHANNELS(PixelGrayA)
DO_WORK_ALL_CHANNELS(PixelRGB)
DO_WORK_ALL_CHANNELS(PixelRGBA)

#define SWITCH_ON_CHANNEL_TYPE( PIXELTYPE )                            \
  switch (fmt.channel_type) {                                          \
  case VW_CHANNEL_UINT8:   do_work_##PIXELTYPE##_uint8();   break;     \
  case VW_CHANNEL_INT8:    do_work_##PIXELTYPE##_int8();    break;     \
  case VW_CHANNEL_UINT16:  do_work_##PIXELTYPE##_uint16();  break;     \
  case VW_CHANNEL_INT16:   do_work_##PIXELTYPE##_int16();   break;     \
  case VW_CHANNEL_FLOAT32: do_work_##PIXELTYPE##_float32(); break;     \
  default:                 do_work_##PIXELTYPE##_float64(); break;     \
  }

void do_work_all_channels(std::string const& input_file_name){

  ImageFormat fmt = vw::image_format(input_file_name);
  try {
    switch (fmt.pixel_format) {
    case VW_PIXEL_GRAY:  SWITCH_ON_CHANNEL_TYPE(PixelGray);  break; // 1 channel
    case VW_PIXEL_GRAYA: SWITCH_ON_CHANNEL_TYPE(PixelGrayA); break; // 1 channel  + alpha channel
    case VW_PIXEL_RGB:   SWITCH_ON_CHANNEL_TYPE(PixelRGB);   break; // 3 channels 
    default:             SWITCH_ON_CHANNEL_TYPE(PixelRGBA);  break; // 3 channels + alpha channel
    }
  }catch (const Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}
    
int main( int argc, char *argv[] ) {

  po::options_description desc("Usage: undistort_image [options] <input image> <camera model> \n\nOptions");
  desc.add_options()
    ("help,h",        "Display this help message")
    ("input-file",    po::value<std::string>(&input_file_name), 
                      "Explicitly specify the input file")
    ("camera-file",    po::value<std::string>(&camera_file_name), 
                      "Explicitly specify the camera file")
    ("output-file,o", po::value<std::string>(&output_file_name)->default_value("output.png"), 
                      "Specify the output file");
  po::positional_options_description p;
  p.add("input-file", 1);
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

  vw::create_out_dir(output_file_name);
  do_work_all_channels(input_file_name);
  
  return 0;
}
