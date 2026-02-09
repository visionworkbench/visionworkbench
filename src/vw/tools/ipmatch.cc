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


/// \file ipmatch.cc
///
/// Loads the previously saved interest points for images from files
///  and then computes a set of matches between each input pair.
///
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Log.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Geometry.h>
#include <vw/Math/RANSAC.h>
#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/FileUtils.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Camera/CameraGeometry.h>
#include <vw/InterestPoint/InterestPoint.h>
#include <vw/InterestPoint/Matcher.h>
#include <vw/InterestPoint/MatcherIO.h>
#include <vw/FileIO/GdalWriteOptions.h>
#include <vw/FileIO/GdalWriteOptionsDesc.h>

#include <vector>
#include <string>
#include <sstream>
#include <iostream>

using namespace vw;
using namespace vw::ip;

#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

// Draw the two images side by side with matching interest points shown with lines.
static void write_match_image(std::string const& out_file_name,
                              std::string const& file1,
                              std::string const& file2,
                              std::vector<InterestPoint> const& matched_ip1,
                              std::vector<InterestPoint> const& matched_ip2) {
  // Skip image pairs with no matches.
  if (matched_ip1.empty())
    return;

  boost::scoped_ptr<vw::DiskImageResource> irsrc1( DiskImageResource::open(file1) );
  boost::scoped_ptr<vw::DiskImageResource> irsrc2( DiskImageResource::open(file2) );

  // Work out the scaling to produce the subsampled images. These
  // values are choosen just allow a reasonable rendering time.
  float sub_scale  = sqrt(1500.0 * 1500.0 / float(irsrc1->format().cols * irsrc1->format().rows));
        sub_scale += sqrt(1500.0 * 1500.0 / float(irsrc2->format().cols * irsrc2->format().rows));
  sub_scale /= 2;
  if ( sub_scale > 1 ) 
    sub_scale = 1;

  // Paste the two scaled input images into a composite image
  mosaic::ImageComposite<PixelRGB<uint8> > composite;
  if ( irsrc1->has_nodata_read() ) {
    composite.insert( pixel_cast_rescale<PixelRGB<uint8> >
                      (resample(apply_mask(normalize(create_mask(DiskImageView<PixelGray<float> >
                                                                 (*irsrc1),
                                                                 irsrc1->nodata_read()))), 
                                sub_scale)),
                      0, 0 );
  } else {
    composite.insert( pixel_cast_rescale<PixelRGB<uint8> >
                      (resample(normalize(DiskImageView<PixelGray<float> >(*irsrc1)), 
                                sub_scale)),
                      0, 0 );
  }
  
  if ( irsrc2->has_nodata_read() ) {
    composite.insert(pixel_cast_rescale<PixelRGB<uint8> >
                     (resample(apply_mask(normalize(create_mask(DiskImageView<PixelGray<float> >
                                                                (*irsrc2),
                                                                irsrc2->nodata_read()))), 
                               sub_scale)),
                     int32(irsrc1->format().cols * sub_scale), 0 );
  } else {
    composite.insert(pixel_cast_rescale<PixelRGB<uint8> >
                     (resample(normalize(DiskImageView<PixelGray<float> >(*irsrc2)), 
                               sub_scale)),
                     int32(irsrc1->format().cols * sub_scale), 0 );
  }
  composite.set_draft_mode( true );
  composite.prepare();

  // Rasterize the composite so that we can draw on it.
  ImageView<PixelRGB<uint8> > comp = composite;

  srand(time(NULL));

  // Draw a line between matching interest points in the two images (in the composite)
  for (size_t k = 0; k < matched_ip1.size(); ++k) {

    Vector2f start(matched_ip1[k].x, matched_ip1[k].y);
    Vector2f end(matched_ip2[k].x+irsrc1->format().cols, matched_ip2[k].y);
    start *= sub_scale;
    end   *= sub_scale;
    float inc_amt = 1/norm_2(end-start);

    // Use a random color for each line
    PixelRGB<uint8> line_color((unsigned char)(rand() % 255),
                               (unsigned char)(rand() % 255),
                               (unsigned char)(rand() % 255));

    for (float r = 0; r < 1.0; r += inc_amt){
      int i = (int)(0.5 + start.x() + r*(end.x()-start.x()));
      int j = (int)(0.5 + start.y() + r*(end.y()-start.y()));
      if (i >=0 && j >=0 && i < comp.cols() && j < comp.rows())
        comp(i,j) = line_color;
    }
  }

  boost::scoped_ptr<vw::DiskImageResource> rsrc( DiskImageResource::create(out_file_name, comp.format()) );
  block_write_image( *rsrc, comp, TerminalProgressCallback( "tools.ipmatch", "Writing Debug:" ) );
}

// See --merge-match-files.
void merge_match_files_fun(std::vector<std::string> const& files) {
  if (files.size() < 2) {
    vw_out() << "At least one input and one output match file must exist.\n";
    return;
  }

  // Index by left ip x and y to avoid repetition
  std::map<std::pair<double, double>, std::pair<vw::ip::InterestPoint, vw::ip::InterestPoint>>
    ip_map;

  for (size_t file_it = 0; file_it + 1 < files.size(); file_it++) {
    std::cout << "Reading: " << files[file_it] << std::endl;

    std::vector<vw::ip::InterestPoint> left_ip, right_ip;
    vw::ip::read_binary_match_file(files[file_it], left_ip, right_ip);
    for (size_t ip_it = 0; ip_it < left_ip.size(); ip_it++) {

      ip_map[std::make_pair(left_ip[ip_it].x, left_ip[ip_it].y)] =
        std::make_pair(left_ip[ip_it], right_ip[ip_it]);
    }
  }

  std::vector<vw::ip::InterestPoint> left_ip(ip_map.size()), right_ip(ip_map.size());
  int count = 0;
  for (auto it = ip_map.begin(); it != ip_map.end(); it++) {
    left_ip[count] = (it->second).first;
    right_ip[count] = (it->second).second;
    count++;
  }
  
  vw_out() << "Writing match file: " << files.back() << std::endl;
  vw::create_out_dir(files.back());
  write_binary_match_file(files.back(), left_ip, right_ip);
}

// TODO(oalexan1): Make all options below use the Options structure
struct Options: public vw::GdalWriteOptions {};

int main(int argc, char** argv) {

  Options opt;
  std::vector<std::string> input_file_names;
  double      matcher_threshold;
  std::string ransac_constraint, distance_metric_in, output_prefix, flann_method;
  float       inlier_threshold;
  int         ransac_iterations;
  bool        merge_match_files;

  po::options_description general_options("Options");
  general_options.add_options()
    ("output-prefix,o",     po::value(&output_prefix)->default_value(""), 
     "Write output files using this prefix.")
    ("matcher-threshold,t", po::value(&matcher_threshold)->default_value(0.8), 
     "Threshold for the separation between closest and next closest interest points.")
    ("flann-method",  po::value(&flann_method)->default_value("kmeans"),
     "Choose the FLANN method for matching interest points. The default 'kmeans' is "
     "slower but deterministic, while 'kdtree' is faster but not deterministic.")
    ("non-flann",
     "Use an implementation of the interest matcher that is not reliant on FLANN.")
    ("distance-metric,m",   po::value(&distance_metric_in)->default_value("L2"), 
     "Distance metric to use.  Choose one of: [L2 (default), Hamming (only for binary types like ORB)].")
    ("ransac-constraint,r", po::value(&ransac_constraint)->default_value("similarity"), 
     "RANSAC constraint type.  Choose one of: [similarity, homography, fundamental, or none].")
    ("inlier-threshold,i",  po::value(&inlier_threshold)->default_value(10), 
     "RANSAC inlier threshold.")
    ("ransac-iterations",   po::value(&ransac_iterations)->default_value(100), 
     "Number of RANSAC iterations.")
    ("debug-image,d",       "Write out debug images.")
    ("merge-match-files", po::value(&merge_match_files)->default_value(false)->implicit_value(true),
     "Given several match files for the same image pair, merge them. The input match files and output match file must be specified in this order. This is an undocumented debug option.")
    ;

  general_options.add(vw::GdalWriteOptionsDescription(opt));
  
  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-files", po::value<std::vector<std::string> >(&input_file_names));

  po::options_description options("Allowed options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-files", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <image files> <vwip files>\n\n";
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    std::cout << "An error occurred while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    vw_out() << usage.str();
    return 1;
  }

  // This is a hidden utility for merging match files for the same image pair
  if (merge_match_files) {
    merge_match_files_fun(input_file_names);
    return 0;
  }
  
  opt.setVwSettingsFromOpt();
  
  std::string distance_metric = distance_metric_in;
  boost::algorithm::to_lower(distance_metric);
  if ((distance_metric != "l2") && (distance_metric != "hamming")) {
    vw_out() << "Error: Did not recognize the distance metric." << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }    

  // Separate the images from the .vwip files
  std::vector<std::string> image_paths, vwip_paths;
  for (size_t it = 0; it < input_file_names.size(); it++) {
    if (fs::path(input_file_names[it]).extension().string() != ".vwip")
      image_paths.push_back(input_file_names[it]);
    else
      vwip_paths.push_back(input_file_names[it]);
  }

  if (image_paths.size() < 2) {
    vw_out() << "Error: Must specify at least two image files.\n\n" << usage.str();
    return 1;
  }
  
  // Auto-fill the .vwip files if not provided
  if (vwip_paths.empty()) {
    for (size_t it = 0; it < image_paths.size(); it++) {
      vwip_paths.push_back(fs::path(image_paths[it]).replace_extension(".vwip").string());
    }
  }
  
  if (image_paths.size() != vwip_paths.size()) {
    vw_out() << "Error: Must have as many .vwip files as images.\n\n" << usage.str();
    return 1;
  }

  size_t num_input_images = image_paths.size();

  // Iterate over combinations of the input files and find interest points in each.
  for (size_t i = 0; i < num_input_images; ++i) {
    for (size_t j = i+1; j < num_input_images; ++j) {

      // Read each file off disk
      std::vector<InterestPoint> ip1, ip2;
      ip1 = read_binary_ip_file(vwip_paths[i]);
      ip2 = read_binary_ip_file(vwip_paths[j]);

      vw_out() << "Matching between " << image_paths[i] << " (" << ip1.size() 
               << " points) and "     << image_paths[j] << " (" << ip2.size() << " points).\n";

      std::vector<InterestPoint> matched_ip1, matched_ip2;

      //std::cout << "IP1 --> \n";
      for (size_t k=0; k<ip1.size(); ++k) {
      //  std::cout << ip1[i].to_string() << "\n";  
        ip1[j].polarity = false; // HACK
        ip2[j].polarity = false;
      }

      vw_out() << "Using distance metric: " << distance_metric_in << std::endl;
      if (!vm.count("non-flann")) {
        vw_out() << "Using FLANN method: " << flann_method << std::endl;
        // Run interest point matcher that uses KDTree algorithm.
        if (distance_metric == "l2") {
          InterestPointMatcher<L2NormMetric, NullConstraint> 
            matcher(flann_method, matcher_threshold);
          matcher(ip1, ip2, matched_ip1, matched_ip2, TerminalProgressCallback( "tools.ipmatch","Matching:"));
        } 
        if (distance_metric == "hamming") {
          InterestPointMatcher<HammingMetric, NullConstraint> 
            matcher(flann_method, matcher_threshold);
          matcher(ip1, ip2, matched_ip1, matched_ip2, TerminalProgressCallback( "tools.ipmatch","Matching:"));
        }
      } else {
        vw_out() << "Not using FLANN.\n";
        // Run interest point matcher that does not use FLANN.
        if (distance_metric == "l2") {
          InterestPointMatcherSimple<L2NormMetric, NullConstraint> matcher(matcher_threshold);
          matcher(ip1, ip2, matched_ip1, matched_ip2, TerminalProgressCallback( "tools.ipmatch","Matching:"));
        }
        if (distance_metric == "hamming") {
          InterestPointMatcherSimple<HammingMetric, NullConstraint> matcher(matcher_threshold);
          matcher(ip1, ip2, matched_ip1, matched_ip2, TerminalProgressCallback( "tools.ipmatch","Matching:"));
        }
      } // End non-flann case

      vw_out() << "Found " << matched_ip1.size() << " matches before duplicate removal.\n";

      remove_duplicates(matched_ip1, matched_ip2);
      vw_out() << "Found " << matched_ip1.size() << " matches.\n";

      std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matched_ip1),
                           ransac_ip2 = iplist_to_vectorlist(matched_ip2);
      std::vector<size_t> indices;
      try {
        // RANSAC is used to fit a transform between the matched sets
        // of points.  Points that don't meet this geometric
        // constraint are rejected as outliers.
        if (ransac_constraint == "similarity") {
          math::RandomSampleConsensus<math::SimilarityFittingFunctor, 
                                      math::InterestPointErrorMetric> 
              ransac( math::SimilarityFittingFunctor(),
                      math::InterestPointErrorMetric(),
                      ransac_iterations,
                      inlier_threshold,
                      ransac_ip1.size()/2, true);
          Matrix<double> H(ransac(ransac_ip1,ransac_ip2));
          std::cout << "\t--> Similarity: " << H << "\n";
          indices = ransac.inlier_indices(H,ransac_ip1,ransac_ip2);
        } else if (ransac_constraint == "homography") {
          math::RandomSampleConsensus<math::HomographyFittingFunctor, 
                                      math::InterestPointErrorMetric> 
              ransac( math::HomographyFittingFunctor(),
                      math::InterestPointErrorMetric(),
                      ransac_iterations,
                      inlier_threshold,
                      ransac_ip1.size()/2, true);
          Matrix<double> H(ransac(ransac_ip1,ransac_ip2));
          std::cout << "\t--> Homography: " << H << "\n";
          indices = ransac.inlier_indices(H,ransac_ip1,ransac_ip2);
        } else if (ransac_constraint == "fundamental") {
          math::RandomSampleConsensus<camera::FundamentalMatrix8PFittingFunctor, 
                                      camera::FundamentalMatrixDistanceErrorMetric> 
              ransac( camera::FundamentalMatrix8PFittingFunctor(),
                      camera::FundamentalMatrixDistanceErrorMetric(), 
                      ransac_iterations, 
                      inlier_threshold, 
                      ransac_ip1.size()/2, true );
          Matrix<double> F(ransac(ransac_ip1,ransac_ip2));
          std::cout << "\t--> Fundamental: " << F << "\n";
          indices = ransac.inlier_indices(F,ransac_ip1,ransac_ip2);
        } else if (ransac_constraint == "none") {
          indices.reserve( matched_ip1.size() );
          for ( size_t i = 0; i < matched_ip1.size(); ++i )
            indices.push_back(i);
        } else {
          vw_out() << "Unknown RANSAC constraint type: " << ransac_constraint
                   << ".  Choose one of: [similarity, homography, fundamental, or none]\n";
          return 1;
        }
      } catch (const vw::math::RANSACErr& e ) {
        vw_out() << "RANSAC Failed: " << e.what() << "\n";
        
        vw_out() << "Consider re-running ipfind with a larger value of --ip-per-tile. "
                 << "If ipfind was called with a binary detector, such as 'orb', ensure this "
                 << "tool is called with '--distance-metric hamming'.\n";
        continue;
      }
      vw_out() << "Found " << indices.size() << " final matches.\n";

      std::vector<InterestPoint> final_ip1, final_ip2;
      BOOST_FOREACH(size_t& index, indices) {
        final_ip1.push_back(matched_ip1[index]);
        final_ip2.push_back(matched_ip2[index]);
      }

      std::string match_file = vw::ip::match_filename(output_prefix, image_paths[i],
                                                      image_paths[j]);
        
      vw::create_out_dir(match_file);

      vw_out() << "Writing match file: " << match_file << std::endl;
      write_binary_match_file(match_file, final_ip1, final_ip2);

      if (vm.count("debug-image")) {
        std::string debug_image = fs::path(match_file).replace_extension(".tif").string();

        vw_out() << "Writing debug image: " << debug_image << std::endl;
        write_match_image(debug_image, image_paths[i], image_paths[j],
                          final_ip1, final_ip2);
      }

    }
  }

  return 0;
}
