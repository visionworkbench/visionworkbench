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


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/tools/geoblend.h>

// Global Variables that are defined by the command line
std::vector<std::string> image_files;
std::string mosaic_name;
std::string output_file_type;
std::string channel_type_str;
bool draft;
unsigned int tilesize;
bool tile_output = false;
unsigned int patch_size, patch_overlap;
float nodata_value;
bool has_nodata_value = false;

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

#define SWITCH_ON_CHANNEL_TYPE( PIXELTYPE ) \
  switch (fmt.channel_type) {                                           \
  case VW_CHANNEL_UINT8:   do_blend_##PIXELTYPE##_uint8();   break;     \
  case VW_CHANNEL_INT16:   do_blend_##PIXELTYPE##_int16();   break;     \
  case VW_CHANNEL_UINT16:  do_blend_##PIXELTYPE##_uint16();  break;     \
  default:                 do_blend_##PIXELTYPE##_float32(); break;     \
  }

int main( int argc, char *argv[] ) {
  try {

    po::options_description general_options("Options");
    general_options.add_options()
      ("mosaic-name,o", po::value(&mosaic_name)->default_value("mosaic"), "Specify base output directory")
      ("output-file-type,t", po::value(&output_file_type)->default_value("tif"), "Output file type")
      ("tile-output", "Output the leaf tiles of a quadtree, instead of a single blended image.")
      ("tiled-tiff", po::value(&tilesize)->default_value(0), "Output a tiled TIFF image, with given tile size (0 disables, TIFF only)")
      ("patch-size", po::value(&patch_size)->default_value(256), "Patch size for tiled output, in pixels")
      ("patch-overlap", po::value(&patch_overlap)->default_value(0), "Patch overlap for tiled output, in pixels")
      ("draft", "Draft mode (no blending)")
      ("ignore-alpha", "Ignore the alpha channel of the input images, and don't write an alpha channel in output.")
      ("nodata-value", po::value(&nodata_value), "Pixel value to use for nodata in input and output (when there's no alpha channel)")
      ("channel-type", po::value(&channel_type_str), "Images' channel type. One of [uint8, uint16, int16, float].")
      ("help,h", "Display this help message");

    po::options_description hidden_options("");
    hidden_options.add_options()
      ("input-files", po::value(&image_files));

    po::options_description options("Allowed Options");
    options.add(general_options).add(hidden_options);

    po::positional_options_description p;
    p.add("input-files", -1);

    std::ostringstream usage;
    usage << "Description: merges several DEMs" << std::endl << std::endl;
    usage << "Usage: geoblend [options] <filename1> <filename2> ..." << std::endl << std::endl;
    usage << general_options << std::endl;

    po::variables_map vm;
    try {
      po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
      po::notify( vm );
    } catch (const po::error& e) {
      std::cout << "An error occured while parsing command line arguments.\n";
      std::cout << "\t" << e.what() << "\n\n";
      std::cout << usage.str();
      return 1;
    }

    if( vm.count("help") ) {
      std::cerr << usage.str() << std::endl;
      return 1;
    }

    draft = false;
    if( vm.count("draft") ) {
      draft = true;
    }

    if( vm.count("input-files") < 1 ) {
      std::cerr << "Error: Must specify at least one input file!" << std::endl << std::endl;
      std::cerr << usage.str();
      return 1;
    }

    if(vm.count("tile-output")) tile_output = true;

    if( patch_size <= 0 ) {
      std::cerr << "Error: The patch size must be a positive number!  (You specified " << patch_size << ".)" << std::endl;
      std::cerr << usage << std::endl;
      return 1;
    }

    if( patch_overlap>=patch_size || patch_overlap%2==1 ) {
      std::cerr << "Error: The patch overlap must be an even number nonnegative number" << std::endl;
      std::cerr << "smaller than the patch size!  (You specified " << patch_overlap << ".)" << std::endl;
      std::cerr << usage << std::endl;
      return 1;
    }

    if(tilesize > 0 && vm.count("tile-output")) {
        std::cerr << "Error: Cannot output both a unified tiled TIFF (single image output) and individual tiles (multiple image output)." << std::endl;
        std::cerr << "\tThe two options --tiled-tiff (> 0) and --tile-output conflict." << std::endl;
        return 1;
    }

    if(output_file_type != "tif") {
        std::cerr << "Error: Cannot output a unified tiled TIFF if output file type is not set to tif." << std::endl;
        std::cerr << usage << std::endl;
        return 1;
    }

    if(vm.count("nodata-value")) has_nodata_value = true;

    ImageFormat fmt = vw::image_format(image_files[0]);

    if (vm.count("channel-type")) {
      fmt.channel_type = channel_name_to_enum(channel_type_str);
      switch (fmt.channel_type) {
        case VW_CHANNEL_UINT8:  case VW_CHANNEL_INT16:
        case VW_CHANNEL_UINT16: case VW_CHANNEL_FLOAT32:
          break;
        default:
          std::cerr << "Error: Bad channel type specified." << std::endl;
          std::cerr << usage.str() << std::endl;
          exit(1);
      }
    }

    if (vm.count("ignore-alpha")) {
      if (fmt.pixel_format == VW_PIXEL_RGBA)  fmt.pixel_format = VW_PIXEL_RGB;
      if (fmt.pixel_format == VW_PIXEL_GRAYA) fmt.pixel_format = VW_PIXEL_GRAY;
    }

    vw_out(vw::VerboseDebugMessage) << "Using pixel type " << pixel_format_name(fmt.pixel_format)
                                    << ":" << channel_type_name(fmt.channel_type) << std::endl;

    switch (fmt.pixel_format) {
    case VW_PIXEL_GRAY:
      SWITCH_ON_CHANNEL_TYPE(PixelGray); break;
    case VW_PIXEL_GRAYA:
      SWITCH_ON_CHANNEL_TYPE(PixelGrayA); break;
    case VW_PIXEL_RGB:
      SWITCH_ON_CHANNEL_TYPE(PixelRGB); break;
    default:
      SWITCH_ON_CHANNEL_TYPE(PixelRGBA); break;
    }

  }
  catch (const std::exception& err) {
    vw::vw_out(vw::ErrorMessage) << "Error: " << err.what() << std::endl << "Aborting!" << std::endl;
    return 1;
  }
  return 0;
}
