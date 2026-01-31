// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

// TODO(oalexan1): Move all logic having to do with colors only, to
// src/vw/Image/Colormap.cc.

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/Core/Functors.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Colormap.h>
#include <vw/Image/Filter.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/tools/Common.h>
#include <vw/FileIO/FileUtils.h>
#include <vw/FileIO/GdalWriteOptions.h>
#include <vw/FileIO/GdalWriteOptionsDesc.h>
#include <vw/Cartography/Hillshade.h>

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>

#include <cstdlib>
#include <vw/Image/Manipulation.h>

namespace fs = boost::filesystem;
namespace po = boost::program_options;
using namespace vw;
using namespace vw::cm;

struct Options: vw::GdalWriteOptions {
  
  // Input
  std::string input_file_name;
  std::string shaded_relief_file_name;

  // Settings
  std::string output_file_name, colormap_style;
  float       nodata_value, min_val, max_val;
  bool        draw_legend;
  
  std::map<float, Vector3u> lut_map;
  
  // For hillshade
  bool hillshade;
  double azimuth, elevation, scale; 
  double hillshade_nodata_value;
  double blur_sigma;
  bool   align_to_georef;

  Options(): draw_legend(false), hillshade(false), azimuth(300), elevation(20),
  scale(0.0), hillshade_nodata_value(std::numeric_limits<double>::quiet_NaN()), 
  blur_sigma(std::numeric_limits<double>::quiet_NaN()), align_to_georef(false) {}
};

template <class PixelT>
ImageViewRef<PixelMask<PixelRGB<uint8>>>
apply_colormap(ImageViewRef<PixelT> input_image,
               float min_val, float max_val,
               bool has_alpha, bool has_nodata,
               float nodata_value,
               std::map<float, Vector3u> & lut_map) {

  ImageViewRef<PixelGray<float>> adj_image =
    pixel_cast<PixelGray<float>>(select_channel(input_image, 0));
  
  // Compute min/max of input image values
  if (min_val == 0 && max_val == 0) {
    min_max_channel_values(create_mask(adj_image, nodata_value),
                            min_val, max_val);
    vw_out() << "\t--> Image color map range: ["
             << min_val << "  " << max_val << "]\n";
  } else {
    vw_out() << "\t--> Using user-specified color map range: ["
             << min_val << "  " << max_val << "]\n";
  }

  // Mask input
  ImageViewRef<PixelMask<PixelGray<float>>> img;
  if (has_alpha)
    img = alpha_to_mask(channel_cast<float>(input_image));
  else if (nodata_value != std::numeric_limits<float>::max()) {
    img = channel_cast<float>(create_mask(adj_image, nodata_value));
  }else if (has_nodata) {
    img = create_mask(adj_image, nodata_value);
  } else {
    img = pixel_cast<PixelMask<PixelGray<float>>>(adj_image);
  }
  
  // Apply colormap
  ImageViewRef<PixelMask<PixelRGB<uint8>>> colorized_image
    = vw::per_pixel_filter(normalize(img, min_val, max_val, 0, 1.0), Colormap(lut_map));

  return colorized_image;
}

template <class PixelT>
ImageViewRef<PixelMask<PixelRGB<uint8>>> do_colormap(Options& opt) {

  // Attempt to extract nodata value
  boost::shared_ptr<vw::DiskImageResource>
    disk_img_rsrc(vw::DiskImageResourcePtr(opt.input_file_name));

  bool has_nodata = disk_img_rsrc->has_nodata_read();
  bool has_alpha = PixelHasAlpha<PixelT>::value;
  
  if (opt.nodata_value != std::numeric_limits<float>::max()) {
    vw_out() << "\t--> Using user-supplied nodata value: " << opt.nodata_value << ".\n";
  } else if (has_nodata) {
    opt.nodata_value = disk_img_rsrc->nodata_read();
    vw_out() << "\t--> Extracted nodata value from file: " << opt.nodata_value << ".\n";
  }

  DiskImageView<PixelT> input_image(opt.input_file_name);
  ImageViewRef<PixelMask<PixelRGB<uint8>>> colorized_image
    = apply_colormap<PixelT>(input_image, opt.min_val, opt.max_val, has_alpha,
                             has_nodata, opt.nodata_value, opt.lut_map);

  return colorized_image;
}

void save_colormap(Options const& opt, 
                   ImageViewRef<PixelMask<PixelRGB<uint8>>> colorized_image) {

  cartography::GeoReference georef;
  bool has_georef = cartography::read_georeference(georef, opt.input_file_name);

  bool has_nodata = false;
  double nodata_val = std::numeric_limits<double>::quiet_NaN();

  if (!opt.shaded_relief_file_name.empty()) { // Using a hillshade file
    vw_out() << "\t--> Incorporating hillshading from: "
             << opt.shaded_relief_file_name << ".\n";
    DiskImageView<PixelGray<float>>
      shaded_relief_image(opt.shaded_relief_file_name);
    ImageViewRef<PixelMask<PixelRGB<uint8>>> shaded_image =
      copy_mask(channel_cast<uint8>(colorized_image*pixel_cast<float>(shaded_relief_image)),
                colorized_image);
    vw_out() << "Writing color-mapped image: " << opt.output_file_name << "\n";
    vw::cartography::block_write_gdal_image(opt.output_file_name, shaded_image,
                                            has_georef, georef, has_nodata, nodata_val, opt,
                                            TerminalProgressCallback("tools.colormap", "Writing:"));

  } else { // Not using a hillshade file
    vw_out() << "Writing color-mapped image: " << opt.output_file_name << "\n";
    vw::cartography::block_write_gdal_image(opt.output_file_name, colorized_image,
                                            has_georef, georef, has_nodata, nodata_val, opt,
                                            TerminalProgressCallback("tools.colormap", "Writing:"));
  }
}

void save_legend(Options const& opt) {
  ImageView<PixelGray<float>> img(100, 500);
  for (int j = 0; j < img.rows(); ++j) {
    float val = float(j) / img.rows();
    for (int i = 0; i < img.cols(); ++i) {
      img(i,j) = val;
    }
  }

  ImageViewRef<PixelMask<PixelRGB<uint8>>> colorized_image
    = vw::per_pixel_filter(img, Colormap(opt.lut_map));
  write_image("legend.png", channel_cast_rescale<uint8>(apply_mask(colorized_image)));
}

void handle_arguments(int argc, char *argv[], Options& opt) {
  po::options_description general_options("");
  general_options.add_options()
    ("output-file,o", po::value(&opt.output_file_name),
                      "Specify the output file.")
    ("colormap-style",po::value(&opt.colormap_style)->default_value("binary-red-blue"),
     "Specify the colormap style. Options: binary-red-blue (default), "
     "jet, black-body, viridis, kindlmann, cubehelix, plasma, inferno, rainbow, turbo. "
     "Or specify the name of a file having the colormap, on each line of "
     "which there must be a normalized or percentage intensity and "
     "the three integer rgb values it maps to.")
    ("nodata-value",  po::value(&opt.nodata_value)->default_value
     (std::numeric_limits<float>::max()),
                      "Remap the nodata default value to the min altitude value.")
    ("min",           po::value(&opt.min_val)->default_value(0),
                      "Minimum height of the color map.")
    ("max",           po::value(&opt.max_val)->default_value(0),
                      "Maximum height of the color map.")
    ("moon",          "Set the min and max values to [-8499 10208] meters, which is "
     "suitable for covering elevations on the Moon.")
    ("mars",          "Set the min and max values to [-8208 21249] meters, which is "
     "suitable for covering elevations on Mars.")
    ("shaded-relief-file,s", po::value(&opt.shaded_relief_file_name),
      "Specify a shaded relief image (grayscale) to apply to the colorized image. "
      "For example, this can be a hillshade image.")
    ("hillshade",    "Create a hillshaded image first, then incorporate it in the "
     "colormap. This is equivalent to using an external file with the "
     "--shaded-relief-file option.")
    ("azimuth,a",       po::value(&opt.azimuth)->default_value(300), 
     "Sets the direction the light source is coming from (in degrees). "
     "Zero degrees is to the right, with positive degree counter-clockwise. "
     "To be used with the --hillshade option.")
    ("elevation,e",     po::value(&opt.elevation)->default_value(20), 
     "Set the elevation of the light source (in degrees). To be used with "
     "the --hillshade option.")
    ("legend",        "Generate the colormap legend.  This image is saved (without "
     "labels) as 'legend.png'.");

  general_options.add(vw::GdalWriteOptionsDescription(opt));
  
  po::options_description positional("");
  positional.add_options()
    ("input-file", po::value(&opt.input_file_name), "Input disparity map.");

  po::positional_options_description positional_desc;
  positional_desc.add("input-file", 1);

  po::options_description all_options;
  all_options.add(general_options).add(positional);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(all_options)
              .positional(positional_desc).run(), vm);
    po::notify(vm);
  } catch (const po::error& e) {
    vw_throw(ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << general_options);
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <input img> \n";

  if (vm.count("help"))
    vw_throw(ArgumentErr() << usage.str() << general_options);
  if (opt.input_file_name.empty())
    vw_throw(ArgumentErr() << "Missing input file!\n"
              << usage.str() << general_options);
  if (vm.count("moon")) {
    opt.min_val = -8499;
    opt.max_val = 10208;
  }
  if (vm.count("mars")) {
    opt.min_val = -8208;
    opt.max_val = 21249;
  }
  if (opt.output_file_name.empty())
    opt.output_file_name =
      fs::path(opt.input_file_name).replace_extension().string() + "_CMAP.tif";
  opt.draw_legend = vm.count("legend");
  opt.hillshade = vm.count("hillshade");

  if (opt.shaded_relief_file_name != "" && opt.hillshade) 
    vw_throw(ArgumentErr() << "Cannot use both --shaded-relief-file and --hillshade.\n");

  opt.setVwSettingsFromOpt();
  create_out_dir(opt.output_file_name);
}

int main(int argc, char *argv[]) {

  Options opt;
  try {
    handle_arguments(argc, argv, opt);

    parse_color_style(opt.colormap_style, opt.lut_map);
    
    // Get the right pixel/channel type.
    ImageFormat fmt = vw::image_format(opt.input_file_name);
    ImageViewRef<PixelMask<PixelRGB<uint8>>> out;

    switch(fmt.pixel_format) {
    case VW_PIXEL_GRAY:
      switch(fmt.channel_type) {
      case VW_CHANNEL_UINT8:  out = do_colormap<PixelGray<uint8>>  (opt); break;
      case VW_CHANNEL_INT16:  out = do_colormap<PixelGray<int16>>  (opt); break;
      case VW_CHANNEL_UINT16: out = do_colormap<PixelGray<uint16>> (opt); break;
      default:                out = do_colormap<PixelGray<float32>>(opt); break;
      }
      break;
    case VW_PIXEL_GRAYA:
      switch(fmt.channel_type) {
      case VW_CHANNEL_UINT8:  out = do_colormap<PixelGrayA<uint8>>  (opt); break;
      case VW_CHANNEL_INT16:  out = do_colormap<PixelGrayA<int16>>  (opt); break;
      case VW_CHANNEL_UINT16: out = do_colormap<PixelGrayA<uint16>> (opt); break;
      default:                out = do_colormap<PixelGrayA<float32>>(opt); break;
      }
      break;
    default:
      vw_throw(ArgumentErr() << "Unsupported pixel format. The image must have only one channel.");
    }
    
    if (opt.hillshade) {
      // Do the hillshade first, then will use it when colorizing.
      // Create opt.shaded_relief_file_name by removing any extension first
      opt.shaded_relief_file_name =
        fs::path(opt.output_file_name).replace_extension().string();
      // Remove _CMAP
      int len = opt.shaded_relief_file_name.size();
      if (len >= 5 && opt.shaded_relief_file_name.substr(len - 5) == "_CMAP")
        opt.shaded_relief_file_name = opt.shaded_relief_file_name.substr(0, len - 5);
      // In either case, append _HILLSHADE.tif
      opt.shaded_relief_file_name = opt.shaded_relief_file_name + "_HILLSHADE.tif";
      // Create the hillshade and save it to disk
      vw::cartography::do_multitype_hillshade(opt.input_file_name,
                                              opt.shaded_relief_file_name,
                                              opt.azimuth, opt.elevation, opt.scale,
                                              opt.nodata_value, opt.blur_sigma,
                                              opt.align_to_georef, opt);
    }

    save_colormap(opt, out);
    
    // Draw legend
    if (opt.draw_legend)
      save_legend(opt);

  } catch (const ArgumentErr& e) {
    vw_out() << e.what() << std::endl;
    return 1;
  } catch (const Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
