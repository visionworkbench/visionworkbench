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

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <vw/Cartography/Hillshade.h>
#include <vw/FileIO/FileUtils.h>
#include <vw/FileIO/GdalWriteOptions.h>
#include <vw/FileIO/GdalWriteOptionsDesc.h>
#include <vw/tools/Common.h>

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace vw;

struct Options : vw::GdalWriteOptions {
  Options() : nodata_value(std::numeric_limits<double>::quiet_NaN()), blur_sigma(std::numeric_limits<double>::quiet_NaN()) {}
  // Input
  std::string input_file_name;

  // Settings
  std::string output_file_name;
  double azimuth, elevation, scale;
  double nodata_value;
  double blur_sigma;
  bool   align_to_georef;
};

void handle_arguments(int argc, char *argv[], Options& opt) {
  po::options_description general_options("Description: Outputs image of a DEM lighted as specified\n\nUsage: hillshade [options] <input file> \n\nOptions");
  general_options.add_options()
    ("input-file", po::value(&opt.input_file_name),
     "Explicitly specify the input file")
    ("output-file,o", po::value(&opt.output_file_name),
     "Specify the output file")
    ("azimuth,a", po::value(&opt.azimuth)->default_value(300),
     "Sets the direction the light source is coming from (in degrees). Zero "
     "degrees is to the right, with positive degree counter-clockwise.")
    ("elevation,e", po::value(&opt.elevation)->default_value(20),
     "Set the elevation of the light source (in degrees).")
    ("scale,s", po::value(&opt.scale)->default_value(0),
     "Set the scale of a pixel (in the same units as the DTM height values).")
    ("nodata-value", po::value(&opt.nodata_value),
     "Remap the DEM default value to the min altitude value.")
    ("blur", po::value(&opt.blur_sigma),
     "Pre-blur the DEM with the specified sigma.")
    ("align-to-georef", po::bool_switch(&opt.align_to_georef),
     "The azimuth is relative to East instead of +x in the image.");

  general_options.add(vw::GdalWriteOptionsDescription(opt));

  po::positional_options_description p;
  p.add("input-file", 1);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(general_options).positional(p).run(), vm);
    po::notify(vm);
  } catch (const po::error& e) {
    vw_throw(ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << general_options);
  }

  if(vm.count("help"))
    vw_throw(ArgumentErr() << general_options);

  if (vm.count("input-file") != 1)
    vw_throw(ArgumentErr() << "Error: Must specify exactly one input file!\n\n" << general_options);

  if (opt.output_file_name.empty())
    opt.output_file_name =
      fs::path(opt.input_file_name).replace_extension().string() + "_HILLSHADE.tif";

  opt.setVwSettingsFromOpt();
  
  create_out_dir(opt.output_file_name);
}

int main(int argc, char *argv[]) {

  Options opt;
  try {
    handle_arguments(argc, argv, opt);
    vw::cartography::do_multitype_hillshade(opt.input_file_name,
                                            opt.output_file_name,
                                            opt.azimuth, opt.elevation, opt.scale,
                                            opt.nodata_value, opt.blur_sigma,
                                            opt.align_to_georef, opt);

  } catch (const ArgumentErr& e) {
    vw_out() << e.what() << std::endl;
    return 1;
  } catch (const Exception& e) {
    vw_out() << e.what() << std::endl;
    return 1;
  }

  return 0;
}
