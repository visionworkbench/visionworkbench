// __BEGIN_LICENSE__
//  Copyright (c) 2006-2012, United States Government as represented by the
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


/// \file ipalign.cc
///
/// Example program demonstrating how to align two images using the
/// Interest Point module.

#include <vw/InterestPoint.h>
#include <vw/Math.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/FileIO.h>

using namespace vw;
using namespace vw::ip;

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem/path.hpp>
namespace fs = boost::filesystem;
#include <boost/foreach.hpp>

#define MAX_POINTS_TO_DRAW 1000

struct Options {
  // Input
  std::vector<std::string> input_filenames;

  // Settings
  std::string output_prefix, interest_operator, descriptor_generator;
  float matcher_threshold, detect_gain, tile_size;
  float inlier_threshold;
  int ransac_iterations;
  bool single_scale, homography, debug_images;
  bool save_intermediate;
};

// ------------------------------------------------------------------

// Draw the interest points and write as an image.
template <class ViewT>
void write_point_image(std::string out_file_name,
                       ImageViewBase<ViewT> const& src,
                       InterestPointList const& points) {

  ImageView<PixelRGB<uint8> > viz = pixel_cast<PixelRGB<uint8> >(channel_cast_rescale<uint8>(src));

  // Draw points into color planes
  int n = 0;
  for (InterestPointList::const_iterator pt = points.begin();
       pt != points.end() && n < MAX_POINTS_TO_DRAW; ++pt, ++n) {
    // Draw a red line from the point outward along the orientation
    for (int r=0; r<(int)(8*(*pt).scale); ++r){
      int i = (int)(0.5 + (*pt).x + r*cos((*pt).orientation));
      int j = (int)(0.5 + (*pt).y + r*sin((*pt).orientation));
      // red, green, blue
      viz(i,j) = PixelRGB<uint8>(255, 0, 0);
    }
    // Draw a green 3x3 filled square at the point to indicate center
    int i0 = (int)(0.5 + (*pt).x);
    int j0 = (int)(0.5 + (*pt).y);
    for (int j=j0-1; j<=j0+1; ++j){
      for (int i=i0-1; i<=i0+1; ++i){
        // red, green, blue
        viz(i,j) = PixelRGB<uint8>(0, 255, 0);
      }
    }
  }

  write_image(out_file_name, viz);
}

// Produces Interest Points from files
void find_interest_points( std::string const& image_name,
                           InterestPointList& ip,
                           Options& opt ) {
  const float IDEAL_LOG_THRESHOLD = .03;
  const float IDEAL_OBALOG_THRESHOLD = .07;
  const float IDEAL_HARRIS_THRESHOLD = 1.2e-5;

  DiskImageView<PixelRGB<uint8> > input_image( image_name );
  ip.clear();
  // Image Alignment
  //
  // Images are aligned by computing interest points, matching
  // them using a standard 2-Norm nearest-neighor metric, and then
  // rejecting outliers by fitting a similarity between the
  // putative matches using RANSAC.
  vw_out(InfoMessage) << "\nInterest Point Detection: "
                      << image_name << "\n";
  if (opt.interest_operator == "harris") {
    HarrisInterestOperator
      interest_operator(IDEAL_HARRIS_THRESHOLD/opt.detect_gain);
    if (!opt.single_scale) {
      ScaledInterestPointDetector<HarrisInterestOperator> detector(interest_operator);
      ip = detect_interest_points(input_image, detector);
    } else {
      InterestPointDetector<HarrisInterestOperator> detector(interest_operator);
      ip = detect_interest_points(input_image, detector);
    }
  } else if (opt.interest_operator == "log" ) {
    // Use a scale-space Laplacian of Gaussian feature detector. The
    // associated threshold is abs(interest) > interest_threshold.
    LogInterestOperator interest_operator(IDEAL_LOG_THRESHOLD/opt.detect_gain);
    if (!opt.single_scale) {
      ScaledInterestPointDetector<LogInterestOperator> detector(interest_operator);
      ip = detect_interest_points(input_image, detector);
    } else {
      InterestPointDetector<LogInterestOperator> detector(interest_operator);
      ip = detect_interest_points(input_image, detector);
    }
  } else if (opt.interest_operator == "obalog" ) {
    OBALoGInterestOperator
      interest_operator(IDEAL_OBALOG_THRESHOLD/opt.detect_gain);
    IntegralInterestPointDetector<OBALoGInterestOperator>
      detector( interest_operator );
    ip = detect_interest_points(input_image, detector);
  } else {
    std::cout << "Unknown interest operator: " << opt.interest_operator << ".  Options are : [ Harris, LoG, OBALoG ]\n";
    exit(0);
  }
  vw_out() << "\t Found " << ip.size() << " interest points.\n";

  // Write out images with interest points marked
  if (opt.debug_images)
    write_point_image(opt.output_prefix + "-debug.tif", input_image, ip);

  // Generate descriptors for interest points.
  vw_out(InfoMessage) << "Generating descriptors using " << opt.descriptor_generator << " generator... ";
  if (opt.descriptor_generator == "patch") {
    PatchDescriptorGenerator descriptor;
    descriptor(input_image, ip);
  } else if ( opt.descriptor_generator == "pca" ) {
    PCASIFTDescriptorGenerator descriptor("pca_basis.exr", "pca_avg.exr");
    descriptor(input_image, ip);
  } else if ( opt.descriptor_generator == "sgrad" ) {
    SGradDescriptorGenerator descriptor;
    descriptor(input_image, ip);
  } else {
    std::cout << "Unknown interest descriptor: " << opt.descriptor_generator << ".  Options are : [ patch, pca, sgrad ]\n";
    exit(0);
  }
  vw_out(InfoMessage) << "done.\n";
  if ( opt.save_intermediate )
    write_binary_ip_file(fs::path(image_name).replace_extension("vwip").string(),ip);
}

void align_images( Options & opt ) {
  // Detect IP from the first image as it is the reference
  InterestPointList ref_ip;
  find_interest_points( opt.input_filenames[0], ref_ip, opt );
  std::vector<InterestPoint> ref_ip_v( ref_ip.size() );
  std::copy( ref_ip.begin(), ref_ip.end(), ref_ip_v.begin() );

  // Determine the size of the reference image
  Vector2i ref_size;
  {
    DiskImageView<PixelRGB<uint8> > ref_image( opt.input_filenames[0] );
    ref_size[0] = ref_image.cols();
    ref_size[1] = ref_image.rows();
  }
  std::string ref_name = opt.input_filenames.front();
  opt.input_filenames.erase( opt.input_filenames.begin() );

  // Detect IP and align all the other images
  BOOST_FOREACH( std::string const& input_name,
                 opt.input_filenames ) {
    InterestPointList input_ip;
    find_interest_points( input_name, input_ip, opt );
    std::vector<InterestPoint> input_ip_v( input_ip.size() );
    std::copy( input_ip.begin(), input_ip.end(), input_ip_v.begin() );

    vw_out(InfoMessage) << "\nInterest Point Matching: "
                        << input_name << "\n";

    DefaultMatcher matcher(opt.matcher_threshold);
    std::vector<InterestPoint> matched_ip1, matched_ip2;
    matcher(ref_ip_v, input_ip_v, matched_ip1, matched_ip2, false,
            TerminalProgressCallback( "tools.ipalign", "Matching:"));
    vw_out(InfoMessage) << "\tFound " << matched_ip1.size() << " putative matches.\n";

    // RANSAC is used to fit a similarity transform between the
    // matched sets of points
    std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matched_ip1);
    std::vector<Vector3> ransac_ip2 = iplist_to_vectorlist(matched_ip2);
    Matrix<double> align_matrix;

    std::vector<size_t> indices;
    if ( opt.homography ) {
      math::RandomSampleConsensus<math::HomographyFittingFunctor, math::InterestPointErrorMetric> ransac(math::HomographyFittingFunctor(), math::InterestPointErrorMetric(), opt.ransac_iterations, opt.inlier_threshold, ransac_ip1.size()/2, true);
      align_matrix = ransac(ransac_ip2,ransac_ip1);
      indices = ransac.inlier_indices(align_matrix,ransac_ip2,ransac_ip1);
    } else {
      math::RandomSampleConsensus<math::AffineFittingFunctor, math::InterestPointErrorMetric> ransac(math::AffineFittingFunctor(), math::InterestPointErrorMetric(), opt.ransac_iterations, opt.inlier_threshold, ransac_ip1.size()/2, true);
      align_matrix = ransac(ransac_ip2,ransac_ip1);
      indices = ransac.inlier_indices(align_matrix,ransac_ip2,ransac_ip1);
    }

    if ( opt.save_intermediate ) {
      std::vector<InterestPoint> final_ip1, final_ip2;
      BOOST_FOREACH( size_t& index, indices ) {
        final_ip1.push_back(matched_ip1[index]);
        final_ip2.push_back(matched_ip2[index]);
      }
      std::string output_filename =
        fs::path(ref_name).replace_extension().string() + "__" +
        fs::path(input_name).stem().string() + ".match";
      write_binary_match_file(output_filename, final_ip1, final_ip2);
    }

    DiskImageView<PixelRGB<uint8> > input_image( input_name );
    ImageViewRef<PixelRGB<uint8> > aligned_image =
      transform(input_image, HomographyTransform(align_matrix),
                ref_size[0], ref_size[1]);
    std::ostringstream ostr;
    ostr << opt.output_prefix << "_" << input_name;
    write_image(ostr.str(), aligned_image,
                TerminalProgressCallback( "tools.ipalign", "Writing:") );
  }
}

// ----------------------------------------------------------------------------

void handle_arguments( int argc, char* argv[], Options& opt ) {
  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h", "Display this help message")
    ("debug-images,d", "Produce additional debugging images as well as the aligned image.")
    ("tile-size", po::value(&opt.tile_size)->default_value(1024),
     "Specify the tile size for detecting interest points.")
    ("output-prefix,o", po::value(&opt.output_prefix)->default_value("aligned"),
     "Specify the output prefix")
    ("save-intermediate,s", "Save working VWIP and Match files")

    // Interest point detector options
    ("detector-gain,g", po::value(&opt.detect_gain)->default_value(1.0),
     "Increasing this number will increase the  gain at which interest points are detected.")
    ("interest-operator", po::value(&opt.interest_operator)->default_value("OBALoG"),
     "Choose an interest metric from [LoG, Harris, OBALoG]")
    ("single-scale", "Do not use the scale-space interest point detector.")

    // Descriptor generator options
    ("descriptor-generator", po::value(&opt.descriptor_generator)->default_value("sgrad"),
     "Choose a descriptor generator from [patch,pca,sgrad]")

    // Alignment options
    ("matcher-threshold,t", po::value(&opt.matcher_threshold)->default_value(0.5),
     "Rejects points during matching if best > matcher_threshold * second_best")
    ("inlier-threshold,i", po::value(&opt.inlier_threshold)->default_value(10), "RANSAC inlier threshold.")
    ("ransac-iterations", po::value(&opt.ransac_iterations)->default_value(100), "Number of RANSAC iterations.")
    ("homography", "Align images using a full projective transform (homography).  By default, aligment uses a more restricted Similarity transform.");

  po::options_description positional("");
  positional.add_options()
    ("input-file",  po::value(&opt.input_filenames),  "Input images");

  po::positional_options_description positional_desc;
  positional_desc.add("input-file", -1);

  po::options_description all_options;
  all_options.add(general_options).add(positional);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << general_options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <ptk-url> <drg-file> <cam-file>\n";

  if ( vm.count("help") )
    vw_throw( ArgumentErr() << usage.str() << general_options );
  if ( opt.input_filenames.empty() )
    vw_throw( ArgumentErr() << "Missing input images!\n"
              << usage.str() << general_options );

  boost::to_lower( opt.interest_operator );
  boost::to_lower( opt.descriptor_generator );
  opt.single_scale = vm.count("single-scale");
  opt.homography = vm.count("homography");
  opt.debug_images = vm.count("debug-images");
  opt.save_intermediate = vm.count("save-intermediate");
}

int main(int argc, char* argv[]) {

  Options opt;
  try {
    handle_arguments( argc, argv, opt );
    align_images( opt );
  } catch ( const ArgumentErr& e ) {
    vw_out() << e.what() << "\n";
    return 1;
  } catch ( const Exception& e ) {
    std::cerr << "\n\nVW Error: " << e.what() << "\n";
    return 1;
  } catch ( const std::bad_alloc& e ) {
    std::cerr << "\n\nError: Ran out of Memory!\n";
    return 1;
  } catch ( const std::exception& e ) {
    std::cerr << "\n\nError: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
