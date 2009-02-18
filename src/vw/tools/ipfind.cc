// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file ipfind.cc
///
/// Finds the interest points in an image and outputs them an Binary
/// (default) or ASCII format.  The ASCII format is compatible with
/// the popular Lowe-SIFT toolchain.
///
#include <vw/Core.h>
#include <vw/InterestPoint.h>
#include <vw/Image.h>
#include <vw/FileIO.h>

using namespace vw;
using namespace vw::ip;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

static std::string prefix_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1) 
    result.erase(index, result.size());
  return result;
}

int main(int argc, char** argv) {
  std::vector<std::string> input_file_names;
  std::string interest_operator, descriptor_generator;
  float harris_threshold, log_threshold, surf_threshold;
  int max_points;
  int tile_size;
  int num_threads;
  ImageView<double> integral;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help", "Display this help message")
    ("num-threads", po::value<int>(&num_threads)->default_value(0), "Set the number of threads for interest point detection.  Setting the num_threads to zero causes ipfind to use the visionworkbench default number of threads.")
    ("tile-size,t", po::value<int>(&tile_size)->default_value(2048), "Specify the tile size for processing interest points. (Useful when working with large images)")
    ("lowe,l", "Save the interest points in an ASCII data format that is compatible with the Lowe-SIFT toolchain.")
    
    // Interest point detector options
    ("interest-operator", po::value<std::string>(&interest_operator)->default_value("LoG"), "Choose an interest point metric from [LoG, Harris]")
    ("log-threshold", po::value<float>(&log_threshold)->default_value(0.03), "Sets the threshold for the Laplacian of Gaussian interest operator")
    ("harris-threshold", po::value<float>(&harris_threshold)->default_value(1e-5), "Sets the threshold for the Harris interest operator")
    ("max-points", po::value<int>(&max_points)->default_value(0), "Set the maximum number of interest points you want returned.  The most \"interesting\" points are selected.")
    ("single-scale", "Turn off scale-invariant interest point detection.  This option only searches for interest points in the first octave of the scale space.")

    // Descriptor generator options
    ("descriptor-generator", po::value<std::string>(&descriptor_generator)->default_value("patch"), "Choose a descriptor generator from [patch,pca]");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-files", po::value<std::vector<std::string> >(&input_file_names));
  
  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-files", -1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <filenames>..." << std::endl << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    vw_out(0) << usage.str();
    return 1;
  }

  if( input_file_names.size() < 1 ) {
    vw_out(0) << "Error: Must specify at least one input file!" << std::endl << std::endl;
    vw_out(0) << usage.str();
    return 1;
  }

  if (num_threads == 0) {
    num_threads = vw_settings().default_num_threads();
  }

  // Checking strings
  boost::to_lower( interest_operator );
  boost::to_lower( descriptor_generator );
  // Determine if interest_operator is legitimate 
  if ( !( interest_operator == "harris" ||
	  interest_operator == "log" ) ) {
    vw_out(0) << "Unknown interest operator: " << interest_operator 
	      << ". Options are : [ Harris, LoG]\n";
    exit(0);
  }
  // Determine if descriptor_generator is legitimate
  if ( !( descriptor_generator == "patch" ||
	  descriptor_generator == "pca" ) ) {
    vw_out(0) << "Unkown descriptor generator: " << descriptor_generator
	      << ". Options are : [ Patch, PCA ]\n";
    exit(0);
  }

  // Iterate over the input files and find interest points in each.
  for (unsigned i = 0; i < input_file_names.size(); ++i) {

    vw_out(0) << "Finding interest points in \"" << input_file_names[i] << "\".\n";
    std::string file_prefix = prefix_from_filename(input_file_names[i]);
    DiskImageView<PixelRGB<float> > image(input_file_names[i]);

    // Detecting Interest Points
    InterestPointList ip;
    if ( interest_operator == "harris" ) {
      HarrisInterestOperator interest_operator(harris_threshold);
      if (!vm.count("single-scale")) {
        ScaledInterestPointDetector<HarrisInterestOperator> detector(interest_operator, max_points);
        ip = detect_interest_points(image, detector, num_threads);
      } else {
        InterestPointDetector<HarrisInterestOperator> detector(interest_operator, max_points);
        ip = detect_interest_points(image, detector, num_threads);
      }
    } else if ( interest_operator == "log") {
      // Use a scale-space Laplacian of Gaussian feature detector. The
      // associated threshold is abs(interest) > interest_threshold.
      LogInterestOperator interest_operator(log_threshold);
      if (!vm.count("single-scale")) {
        ScaledInterestPointDetector<LogInterestOperator> detector(interest_operator, max_points);
        ip = detect_interest_points(image, detector, num_threads);
      } else {
        InterestPointDetector<LogInterestOperator> detector(interest_operator, max_points);
        ip = detect_interest_points(image, detector, num_threads);
      }
    }
  
    vw_out(0) << "\t Found " << ip.size() << " points.\n";

    // Generate descriptors for interest points.
    vw_out(InfoMessage) << "\tRunning " << descriptor_generator << " descriptor generator.\n";
    if (descriptor_generator == "patch") {
      PatchDescriptorGenerator descriptor;
      descriptor(image, ip);
    } else if (descriptor_generator == "pca") {
      PCASIFTDescriptorGenerator descriptor("pca_basis.exr", "pca_avg.exr");
      descriptor(image, ip);
    }

    // If ASCII output was requested, write it out.  Otherwise stick
    // with binary output.
    if (vm.count("lowe"))
      write_lowe_ascii_ip_file(file_prefix + ".key", ip);
    else 
      write_binary_ip_file(file_prefix + ".vwip", ip);
  }
}
