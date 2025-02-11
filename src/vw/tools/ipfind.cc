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

/// \file ipfind.cc
///
/// Finds the interest points in an image and outputs them an Binary
/// (default) or ASCII format.  The ASCII format is compatible with
/// the popular Lowe-SIFT toolchain.
///

#include <vw/Core/System.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>
#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/InterestPoint/Descriptor.h>
#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/IntegralDescriptor.h>
#include <vw/InterestPoint/IntegralDetector.h>
#include <vw/InterestPoint/IntegralInterestOperator.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/InterestPointUtils.h>
#include <vw/InterestPoint/MatcherIO.h>
#include <vw/FileIO/FileUtils.h>
#include <vw/FileIO/GdalWriteOptions.h>

using namespace vw;
using namespace vw::ip;

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>
namespace po = boost::program_options;

struct Options: public vw::GdalWriteOptions {
  std::vector<std::string> input_file_names;
  std::string output_folder, interest_operator, descriptor_generator;
  float ip_gain;
  int ip_per_image, ip_per_tile;
  int nodata_radius, print_num_ip, debug_image;
  bool no_orientation;
};
  
int main(int argc, char** argv) {

  const float IDEAL_LOG_THRESHOLD    = .03;
  const float IDEAL_OBALOG_THRESHOLD = .07;
  const float IDEAL_HARRIS_THRESHOLD = 1.2e-5;

  Options opt;

  po::options_description general_options("Options");
  general_options.add_options()
    ("interest-operator",    po::value(&opt.interest_operator)->default_value("sift"), 
     "Choose an interest point detector from: sift (default), orb, OBALoG, LoG, Harris, IAGD.")
    ("descriptor-generator", po::value(&opt.descriptor_generator)->default_value("sift"), 
     "Choose a descriptor generator from: sift (default), orb, sgrad, sgrad2, patch, pca. Some descriptors work only with certain interest point operators (for example, for OBALoG use sgrad, sgrad2, patch, and pca).")
    ("ip-per-image",           po::value(&opt.ip_per_image)->default_value(0), 
     "Set the maximum number of IP to find in the whole image. If not specified, use instead --ip-per-tile.")
    ("ip-per-tile",           po::value(&opt.ip_per_tile)->default_value(250), 
     "Set the maximum number of IP to find in each tile. The tile size is set with --tile-size.")
    ("gain,g",               po::value(&opt.ip_gain)->default_value(1.0), 
     "Increasing this number will increase the gain at which interest points are detected. Default: 1.")
    ("single-scale", "Turn off scale-invariant interest point detection. This option only searches for interest points in the first octave of the scale space. Harris and LoG only.")
    ("no-orientation",       po::bool_switch(&opt.no_orientation), "Turn off rotational invariance.")
    ("normalize",     "Normalize the input. Use for images that have non-standard values such as ISIS cube files.")
    ("per-tile-normalize",     "Individually normalize each processing tile (only with OpenCV).")
    ("nodata-radius", po::value(&opt.nodata_radius)->default_value(1), 
     "Don't detect IP within this many pixels of image borders or nodata. Default: 1.")
    ("output-folder",  po::value(&opt.output_folder)->default_value(""), 
     "Write output files to this location.")
    // Debug options
    ("debug-image,d", po::value(&opt.debug_image)->default_value(0), 
     "Write out a low-resolution or full-resolution debug image with interest points on it if the value of this flag is respectively 1 or 2. Default: 0 (do nothing).")
    ("print-ip",             po::value(&opt.print_num_ip)->default_value(0), 
     "Print information for this many detected IP. Default: 0.")
    ("lowe",          "Save the interest points in an ASCII data format that is compatible with the Lowe-SIFT toolchain.");

  general_options.add(vw::GdalWriteOptionsDescription(opt));

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-files", po::value(&opt.input_file_names));

  po::options_description options("Allowed options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-files", -1);

  std::ostringstream usage;
  usage << "Usage: ipfind [options] <images>\n\n";
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(options).positional(p).run(), vm);
    po::notify(vm);
  } catch (const po::error& e) {
    vw_out() << "An error occurred while parsing command line arguments.\n";
    vw_out() << "\t" << e.what() << "\n\n";
    vw_out() << usage.str();
    return 1;
  }

  opt.setVwSettingsFromOpt();

  if (vm.count("help")) {
    vw_out() << usage.str();
    return 1;
  }

  if (opt.input_file_names.size() < 1) {
    vw_out() << "Error: Must specify at least one input file!" << std::endl << std::endl;
    vw_out() << usage.str();
    return 1;
  }


  if (opt.nodata_radius < 1) {
    vw_out() << "Error: nodata-radius must be >= 1!" << std::endl << std::endl;
    vw_out() << usage.str();
    return 1;
  }

  bool opencv_normalize = false;
  if (vm.count("per-tile-normalize"))
    opencv_normalize = true;

  // Checking strings
  boost::to_lower(opt.interest_operator);
  boost::to_lower(opt.descriptor_generator);

  // The OpenCV detector types are handled differently than the custom detector types.  
  const bool detector_is_opencv = ((opt.interest_operator == "brisk") ||
                                    (opt.interest_operator == "surf") ||
                                    (opt.interest_operator == "orb" ) ||
                                    (opt.interest_operator == "sift")  );
  const bool descriptor_is_opencv = ((opt.descriptor_generator == "brisk") ||
                                      (opt.descriptor_generator == "surf") ||
                                      (opt.descriptor_generator == "orb" ) ||
                                      (opt.descriptor_generator == "sift")  );

  if (opencv_normalize && !detector_is_opencv) {
    vw_out() << "Cannot use per-tile normalize with a non-OpenCV detector!\n";
    exit(1);
  }

  // Determine if opt.interest_operator is legitimate
  if (!(opt.interest_operator == "iagd"   ||
          opt.interest_operator == "harris" ||
          opt.interest_operator == "log"    ||
          opt.interest_operator == "obalog" || 
          //opt.interest_operator == "brisk"  ||
          opt.interest_operator == "orb"    ||
          //opt.interest_operator == "surf"   ||
          opt.interest_operator == "sift"    )) {
    vw_out() << "Unknown interest operator: " << opt.interest_operator
             << ". Options are: sift, orb, OBALoG, LoG, Harris, IAGD.\n";
    exit(1);
  }

  // Determine if opt.descriptor_generator is legitimate
  if (!(opt.descriptor_generator == "patch"  ||
          opt.descriptor_generator == "pca"    ||
          opt.descriptor_generator == "sgrad"  ||
          opt.descriptor_generator == "sgrad2" ||
          //opt.descriptor_generator == "brisk"  ||
          opt.descriptor_generator == "orb"    ||
          //opt.descriptor_generator == "surf"   ||
          opt.descriptor_generator == "sift"   )) {
    vw_out() << "Unkown descriptor generator: " << opt.descriptor_generator
             << ". Options are: sift, orb, sgrad, sgrad2, patch, pca.\n";
    exit(1);
  }

  if ((detector_is_opencv && !descriptor_is_opencv) &&
      (opt.descriptor_generator != opt.interest_operator)) {
    vw_out() <<"The value of --descriptor-generator is '" << opt.descriptor_generator << "'. ";
    opt.descriptor_generator = opt.interest_operator;
    vw_out() << "Switching it to '" << opt.descriptor_generator << "' to match --interest-operator.\n";
  }

  if (!detector_is_opencv && descriptor_is_opencv) {
    vw_out() <<"The value of --descriptor-generator is '" << opt.descriptor_generator << "'. ";
    opt.descriptor_generator = "sgrad";
    vw_out() << "Switching it to '" << opt.descriptor_generator << "' to match --interest-operator.\n";
    vw_out() << "Can use any of: 'sgrad', 'sgrad2', 'patch', 'pca'.\n";
  }

  // Iterate over the input files and find interest points in each.
  for (size_t i = 0; i < opt.input_file_names.size(); i++) {
    vw_out() << "Finding interest points in \"" << opt.input_file_names[i] << "\".\n";
    if (opt.ip_per_image > 0) {
      vw_out() << "Computing " << opt.ip_per_image << " ip per image.\n";
      int tile_size = vw_settings().default_tile_size();
      Vector2 image_size = vw::file_image_size(opt.input_file_names[i]);
      double num_tiles = image_size[0]*image_size[1]/(double(tile_size*tile_size));
      opt.ip_per_tile = (int)ceil(opt.ip_per_image / num_tiles);
      opt.ip_per_tile = std::max(1, opt.ip_per_tile);
    }
   vw_out() << "Computing " << opt.ip_per_tile << " ip per tile.\n";
    
    std::string file_prefix = boost::filesystem::path(opt.input_file_names[i]).replace_extension().string();
    if (opt.output_folder != "") {
      // Write the output files to the specified output folder.
      boost::filesystem::path of(opt.output_folder);
      of /= boost::filesystem::path(file_prefix).filename();
      file_prefix = of.string();

      vw::create_out_dir(file_prefix);
    }

    // Get reference to the file and convert to floating point
    boost::shared_ptr<vw::DiskImageResource> image_rsrc(vw::DiskImageResourcePtr(opt.input_file_names[i]));
    DiskImageView<PixelGray<float>> raw_image(*image_rsrc);
    ImageViewRef<PixelGray<float>> image;
    ImageViewRef<PixelMask<PixelGray<float>>> masked_image;
    // Load nodata value if present.
    double nodata = 0;
    bool   has_nodata = image_rsrc->has_nodata_read();
    if (has_nodata) {
      nodata = image_rsrc->nodata_read();
      vw_out() << "Read in nodata value: " << nodata << std::endl;
    }

    if (vm.count("normalize")) {
      vw_out() << "Normalizing the input image...\n";
      if (has_nodata) {
        masked_image = normalize(create_mask(raw_image, nodata));
        image        = apply_mask(masked_image); // A backup option for detectors which don't handle nodata.
      }
      else
        image = normalize(raw_image);
    } else { // Not normalizing
      if (has_nodata) {
        masked_image = create_mask(raw_image, nodata);
        image        = apply_mask(masked_image); // A backup option for detectors which don't handle nodata.
      }
      else
        image = raw_image;
    }

    // Potentially mask image on a no data value
    if (has_nodata)
      vw_out(DebugMessage,"interest_point") << "Image has a nodata value: " << nodata << "\n";

    const bool describeInDetect = true;
  
    // Detecting Interest Points
    // - Note that only the OpenCV detectors handle masks.
    InterestPointList ip;
    if (opt.interest_operator == "harris") {
      // Harris threshold is inversely proportional to gain.
      HarrisInterestOperator interest_operator(IDEAL_HARRIS_THRESHOLD/opt.ip_gain);
      if (!vm.count("single-scale")) {
        ScaledInterestPointDetector<HarrisInterestOperator> detector(interest_operator, opt.ip_per_tile);
        ip = detect_interest_points(image, detector, opt.ip_per_tile);
      } else {
        InterestPointDetector<HarrisInterestOperator> detector(interest_operator, opt.ip_per_tile);
        ip = detect_interest_points(image, detector, opt.ip_per_tile);
      }
    } else if (opt.interest_operator == "log") {
      // Use a scale-space Laplacian of Gaussian feature detector. The
      // associated threshold is abs(interest) > interest_threshold.
      // LoG threshold is inversely proportional to gain..
      LogInterestOperator interest_operator(IDEAL_LOG_THRESHOLD/opt.ip_gain);
      if (!vm.count("single-scale")) {
        ScaledInterestPointDetector<LogInterestOperator> detector(interest_operator, opt.ip_per_tile);
        ip = detect_interest_points(image, detector, opt.ip_per_tile);
      } else {
        InterestPointDetector<LogInterestOperator> detector(interest_operator, opt.ip_per_tile);
        ip = detect_interest_points(image, detector, opt.ip_per_tile);
      }
    } else if (opt.interest_operator == "obalog") {
      // OBALoG threshold is inversely proportional to gain ..
      OBALoGInterestOperator interest_operator(IDEAL_OBALOG_THRESHOLD/opt.ip_gain);
      IntegralInterestPointDetector<OBALoGInterestOperator> detector(interest_operator, opt.ip_per_tile);
      ip = detect_interest_points(image, detector, opt.ip_per_tile);
    } else if (opt.interest_operator == "iagd") {
      // This is the default ASP implementation
      IntegralAutoGainDetector detector(opt.ip_per_tile);
      ip = detect_interest_points(image, detector, opt.ip_per_tile);
#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1
    } else if (detector_is_opencv) {

      OpenCvIpDetectorType ocv_type = OPENCV_IP_DETECTOR_TYPE_SIFT;
      if (opt.interest_operator == "brisk") {
        ocv_type = OPENCV_IP_DETECTOR_TYPE_BRISK;
      } else if (opt.interest_operator == "orb") {
        ocv_type = OPENCV_IP_DETECTOR_TYPE_ORB;
      } else if (opt.interest_operator == "sift") {
        ocv_type = OPENCV_IP_DETECTOR_TYPE_SIFT;
      } else if (opt.interest_operator == "surf") {
        ocv_type = OPENCV_IP_DETECTOR_TYPE_SURF;
      }
      OpenCvInterestPointDetector detector(ocv_type, opencv_normalize, describeInDetect, opt.ip_per_tile);
      if (has_nodata)
        ip = detect_interest_points(masked_image, detector, opt.ip_per_tile);
      else
        ip = detect_interest_points(image, detector, opt.ip_per_tile);
    }
#else // End OpenCV section
    } else {
      vw_out() << "Error: Cannot use ORB, BRISK, SIFT, or SURF if not compiled with OpenCV!" << std::endl;
      exit(1); 
    }
#endif

    vw_out() << "Detected " << ip.size() << " raw keypoints.\n";
    if (ip.size() == 0) {
      vw_out() << "No IP, quitting this image. Perhaps the --normalize and/or --ip-per-tile options could be used.\n";
      continue;
    }

    // Removing Interest Points on nodata or within 1/px
    if (opt.nodata_radius > 0) {
      size_t before_size = ip.size();
      remove_ip_near_nodata(raw_image, nodata, ip, opt.nodata_radius);
      
      vw_out(InfoMessage,"interest_point") << "Removed " << before_size-ip.size() << " points close to nodata.\n";
    } // End clearing IP's near nodata

    vw_out() << "\t Found " << ip.size() << " points.\n";

    // Delete the orientation value if requested
    if (opt.no_orientation) {
      BOOST_FOREACH(InterestPoint& i, ip) {
        i.orientation = 0;
      }
    }

    // Generate descriptors for interest points.
    vw_out(InfoMessage) << "\tRunning " << opt.descriptor_generator << " descriptor generator.\n";
    if (opt.descriptor_generator == "patch") {
      PatchDescriptorGenerator descriptor;
      describe_interest_points(image, descriptor, ip);
    } else if (opt.descriptor_generator == "pca") {
      PCASIFTDescriptorGenerator descriptor("pca_basis.exr", "pca_avg.exr");
      describe_interest_points(image, descriptor, ip);
    } else if (opt.descriptor_generator == "sgrad") {
      SGradDescriptorGenerator descriptor;
      describe_interest_points(image, descriptor, ip);
    } else if (opt.descriptor_generator == "sgrad2") {
      SGrad2DescriptorGenerator descriptor;
      describe_interest_points(image, descriptor, ip);
#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1
    // Nothing happens here because we compute the descriptors at the same time we detect IPs.
    }
#else
    } else {
      vw_out() << "Error: Cannot use SIFT, SURF, ORB or BRISK if not compiled with OpenCV!\n";
      exit(1); 
    }
#endif

    // Print raw IP info if requested
    int limit = opt.print_num_ip;
    for (InterestPointList::const_iterator iter=ip.begin(); iter!=ip.end(); ++iter) {
      if (limit <= 0)
        break;
      vw_out() << iter->to_string() << "\n";  
      --limit;
    }

    // If ASCII output was requested, write it out.  Otherwise stick with binary output.
    if (vm.count("lowe")) {
      vw_out() << "Writing output file " << file_prefix + ".key" << std::endl;
      write_lowe_ascii_ip_file(file_prefix + ".key", ip);
    } else {
      vw_out() << "Writing output file " << file_prefix + ".vwip" << std::endl;
      write_binary_ip_file(file_prefix + ".vwip", ip);
    }

    // Write Debug image
    if (opt.debug_image > 0) {
      const int WRITE_FULL_RES_DEBUG = 2;
      write_ip_debug_image(file_prefix + "_debug.png",
                           raw_image, ip, has_nodata, nodata,
                           (opt.debug_image==WRITE_FULL_RES_DEBUG),
                           // Orb has large scale values, so shrink them
                           (opt.interest_operator == "orb")); 
    }
                         
  } // End loop through input files

  return 0;
}
