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


/// \file ipalign.cc
///
/// Example program demonstrating how to align two images using the
/// Interest Point module.

#include <vw/Core/Log.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/Geometry.h>
#include <vw/Math/RANSAC.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/Transform.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/InterestPoint/Descriptor.h>
#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/Matcher.h>
#include <vw/InterestPoint/IntegralDetector.h>
#include <vw/InterestPoint/IntegralInterestOperator.h>
#include <vw/Camera/CameraUtilities.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/LensDistortion.h>

using namespace vw;
using namespace vw::ip;
using namespace vw::camera;

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem/path.hpp>
namespace fs = boost::filesystem;
#include <boost/foreach.hpp>

#define MAX_POINTS_TO_DRAW 1000

  // Image Alignment
  //
  // Images are aligned by computing interest points, matching
  // them using a standard 2-Norm nearest-neighor metric, and then
  // rejecting outliers by fitting a similarity between the
  // putative matches using RANSAC.

/// Program parameters
struct Options {
  // Input
  std::string left_image_file,  right_image_file;
  std::string left_camera_file, right_camera_file;

  // Settings
  std::string output_prefix, interest_operator, descriptor_generator;
  std::string align_method;
  float matcher_threshold, detect_gain, tile_size;
  float inlier_threshold;
  int   ransac_iterations;
  bool  single_scale, debug_images;
  bool  save_intermediate;
};

// ------------------------------------------------------------------

/// Draw the interest points and write as an image.
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



/// Produces Interest Points from files
void find_interest_points( std::string const& image_name, InterestPointList& ip, Options const& opt ) {
  const float IDEAL_LOG_THRESHOLD    = .03;
  const float IDEAL_OBALOG_THRESHOLD = .07;
  const float IDEAL_HARRIS_THRESHOLD = 1.2e-5;

  DiskImageView<PixelRGB<uint8> > input_image( image_name );
  ip.clear();

  // Run an interest point detector using the selected parameters.
  // - Templated types require a lot of if statements to process
  vw_out(InfoMessage) << "\nInterest Point Detection: "  << image_name << "\n";
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

  // Generate descriptors for interest points using selected method.
  vw_out(InfoMessage) << "Generating descriptors using " << opt.descriptor_generator << " generator... " << std::flush;
  if (opt.descriptor_generator == "patch") {
    PatchDescriptorGenerator descriptor;
    describe_interest_points( input_image, descriptor, ip );
  } else if ( opt.descriptor_generator == "pca" ) {
    PCASIFTDescriptorGenerator descriptor("pca_basis.exr", "pca_avg.exr");
    describe_interest_points( input_image, descriptor, ip );
  } else if ( opt.descriptor_generator == "sgrad" ) {
    SGradDescriptorGenerator descriptor;
    describe_interest_points( input_image, descriptor, ip );
  } else {
    std::cout << "Unknown interest descriptor: " << opt.descriptor_generator << ".  Options are : [ patch, pca, sgrad ]\n";
    exit(0);
  }
  vw_out(InfoMessage) << "done." << std::endl;
  if ( opt.save_intermediate ) // Optionally write out a binary interest point file
    write_binary_ip_file(fs::path(image_name).replace_extension("vwip").string(), ip);
} // End find_interest_points




/// Find the transform between image pairs and write transformed images to disk
void align_images( Options & opt ) {

  // Determine the size of the reference image
  Vector2i ref_size; 
  DiskImageView<PixelRGB<uint8> > ref_image  ( opt.left_image_file  );
  DiskImageView<PixelRGB<uint8> > input_image( opt.right_image_file );
  ref_size[0] = ref_image.cols();
  ref_size[1] = ref_image.rows();

  // TODO: Refactor some code into functions!

  Matrix<double> align_matrix;  
  if ( opt.align_method != "epipolar" ) {
    // Find IPs
    InterestPointList ref_ip, input_ip;
    find_interest_points( opt.left_image_file, ref_ip, opt );
    find_interest_points( opt.right_image_file, input_ip, opt );
    
    // Convert from InterestPointList to a std::vector of InterestPoints
    std::vector<InterestPoint> ref_ip_v( ref_ip.size() ),
                               input_ip_v( input_ip.size() );
    std::copy( ref_ip.begin(),   ref_ip.end(),   ref_ip_v.begin() );
    std::copy( input_ip.begin(), input_ip.end(), input_ip_v.begin() );

    // Use an interest point matcher to find matched pairs of interest points
    DefaultMatcher matcher(opt.matcher_threshold);
    std::vector<InterestPoint> matched_ip1, matched_ip2;
    matcher(ref_ip_v, input_ip_v, matched_ip1, matched_ip2,
            TerminalProgressCallback( "tools.ipalign", "Matching:"));
    vw_out(InfoMessage) << "\tFound " << matched_ip1.size() << " putative matches.\n";

    // RANSAC is used to fit a similarity transform between the matched sets of points
    // - We generate the transformation matrix and the lists of inlier indices
    std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matched_ip1),
                         ransac_ip2 = iplist_to_vectorlist(matched_ip2);
    std::vector<size_t>  indices;

    if ( opt.align_method == "homography" ) { // Full projective transform
      math::RandomSampleConsensus<math::HomographyFittingFunctor, math::InterestPointErrorMetric> 
                  ransac(math::HomographyFittingFunctor(), math::InterestPointErrorMetric(), opt.ransac_iterations, 
                         opt.inlier_threshold, ransac_ip1.size()/2, true);
      align_matrix = ransac(ransac_ip2, ransac_ip1);
      indices      = ransac.inlier_indices(align_matrix, ransac_ip2, ransac_ip1);
    }
    if ( opt.align_method == "similarity" ) { // Similarity transform
      math::RandomSampleConsensus<math::AffineFittingFunctor, math::InterestPointErrorMetric> 
                  ransac(math::AffineFittingFunctor(), math::InterestPointErrorMetric(), opt.ransac_iterations, 
                         opt.inlier_threshold, ransac_ip1.size()/2, true);
      align_matrix = ransac(ransac_ip2, ransac_ip1);
      indices      = ransac.inlier_indices(align_matrix, ransac_ip2, ransac_ip1);
    }
    
    if ( opt.save_intermediate ) { // Save the list of matched pixels to a binary file
      std::string output_ip_filename = opt.output_prefix + "-left_right_ip.match";
      std::vector<InterestPoint> final_ip1, final_ip2;
      BOOST_FOREACH( size_t& index, indices ) {
        final_ip1.push_back(matched_ip1[index]);
        final_ip2.push_back(matched_ip2[index]);
      }
      write_binary_match_file(output_ip_filename, final_ip1, final_ip2);
    }
  } // End finding align matrix case


  // Write the transformed input image to a file
  ImageViewRef< PixelRGB<uint8> > left_aligned, right_aligned;
  if ( opt.align_method == "epipolar" ) {

    ImageView< PixelRGB<uint8> > left_aligned2, right_aligned2;

    //// Load input pinhole cameras
    PinholeModel left_pin (opt.left_camera_file );
    PinholeModel right_pin(opt.right_camera_file);

    // Generate epipolar-aligned camera models
    boost::shared_ptr<PinholeModel> epipolar_left_pin (new PinholeModel);
    boost::shared_ptr<PinholeModel> epipolar_right_pin(new PinholeModel);
    epipolar(left_pin, right_pin, *(epipolar_left_pin.get()), *(epipolar_right_pin.get()));

    // Expand epipolar cameras to contain the entire source images.
    Vector2i epi_size1, epi_size2;
    BBox2i left_input_roi (Vector2i(0,0), ref_image.get_size()  ),
           right_input_roi(Vector2i(0,0), input_image.get_size());
    resize_epipolar_cameras_to_fit(left_pin, right_pin,
                                   *(epipolar_left_pin.get()), *(epipolar_right_pin.get()),
                                   left_input_roi, right_input_roi,
                                   epi_size1, epi_size2);


    // Zero edge extension. Not ideal. One better create a mask or something.
    PixelRGB<uint8> Zero = PixelRGB<uint8>(0);
    get_epipolar_transformed_pinhole_images(opt.left_camera_file,   opt.right_camera_file,
                                            epipolar_left_pin, epipolar_right_pin,
                                            ref_image, input_image,
                                            left_input_roi, right_input_roi,
                                            epi_size1, epi_size2,
                                            left_aligned2, right_aligned2,
                                            ValueEdgeExtension< PixelRGB<uint8> >(Zero),
                                            BilinearInterpolation());

    std::string output_path = opt.output_prefix + "_right_image.tif";
    write_image(output_path, crop(edge_extend(right_aligned2,ZeroEdgeExtension()),bounding_box(left_aligned2)),
                TerminalProgressCallback( "tools.ipalign", "Writing right:") );
    output_path = opt.output_prefix + "_left_image.tif";    
    write_image(output_path, left_aligned2,  TerminalProgressCallback( "tools.ipalign", "Writing left :") );

  }
  else { // Non-epipolar case
    std::string output_path = opt.output_prefix + "_transformed_right_image.tif";    
    right_aligned = transform(input_image, HomographyTransform(align_matrix),
                              ref_size[0], ref_size[1]);
    write_image(output_path, right_aligned, TerminalProgressCallback( "tools.ipalign", "Writing:") );
  }
                
  
}

// ----------------------------------------------------------------------------

void handle_arguments( int argc, char* argv[], Options& opt ) {
  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h", "Display this help message")
    ("debug-images,d", po::bool_switch(&opt.debug_images),
                       "Produce additional debugging images as well as the aligned image.")
    ("tile-size", po::value(&opt.tile_size)->default_value(1024),
                  "Specify the tile size for detecting interest points.")
    ("output-prefix,o", po::value(&opt.output_prefix)->default_value("aligned"),
                        "Specify the output prefix")
    ("save-intermediate,s", po::bool_switch(&opt.save_intermediate),
                            "Save working VWIP and Match files")

    // Interest point detector options
    ("detector-gain,g", po::value(&opt.detect_gain)->default_value(1.0),
                        "Increasing this number will increase the  gain at which interest points are detected.")
    ("interest-operator", po::value(&opt.interest_operator)->default_value("OBALoG"),
                          "Choose an interest metric from [LoG, Harris, OBALoG]")
    ("single-scale", po::bool_switch(&opt.single_scale),
                     "Do not use the scale-space interest point detector.")

    // Descriptor generator options
    ("descriptor-generator", po::value(&opt.descriptor_generator)->default_value("sgrad"),
                             "Choose a descriptor generator from [patch, pca, sgrad]")

    // Alignment options
    ("matcher-threshold,t", po::value(&opt.matcher_threshold)->default_value(0.5),
                            "Rejects points during matching if best > matcher_threshold * second_best")
    ("inlier-threshold,i", po::value(&opt.inlier_threshold)->default_value(10), 
                           "RANSAC inlier threshold.")
    ("ransac-iterations", po::value(&opt.ransac_iterations)->default_value(100), 
                          "Number of RANSAC iterations.")
    ("align-method", po::value(&opt.align_method)->default_value("similarity"),
                   "Choose similarity, homography, or epipolar image alignment.");

  po::options_description positional("");
  positional.add_options()
    ("left-image",   po::value(&opt.left_image_file),  "Left image"  )
    ("right-image",  po::value(&opt.right_image_file), "Right image" )
    ("left-camera",  po::value(&opt.left_camera_file), "Left camera" )
    ("right-camera", po::value(&opt.right_camera_file),"Right camera");

  po::positional_options_description positional_desc;
  positional_desc.add("left-image",   1);
  positional_desc.add("right-image",  1);
  positional_desc.add("left-camera",  1);
  positional_desc.add("right-camera", 1);

  po::options_description all_options;
  all_options.add(general_options).add(positional);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t" << e.what() << general_options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <left image> <right image> <left camera model> <right camera model>\n";

  if ( vm.count("help") )
    vw_throw( ArgumentErr() << usage.str() << general_options );

  boost::to_lower(opt.interest_operator);
  boost::to_lower(opt.descriptor_generator);
  boost::to_lower(opt.align_method);

  if (opt.left_image_file.empty() || opt.right_image_file.empty())
      vw_throw( ArgumentErr() << "Missing required image files.\n\n" << usage );

  if ((opt.align_method == "epipolar") && (opt.left_camera_file.empty() || opt.right_camera_file.empty()))
      vw_throw( ArgumentErr() << "Missing required camera model files.\n\n" << usage );
}


int main(int argc, char* argv[]) {

  Options opt;
  handle_arguments( argc, argv, opt ); // Load user arguments
  align_images( opt );                 // Do all of the work!
  
  return 0;
}
