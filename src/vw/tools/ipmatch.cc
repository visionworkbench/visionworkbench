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
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Camera/CameraGeometry.h>
#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/Matcher.h>

#include <vector>
#include <string>
#include <sstream>
#include <iostream>

using namespace vw;
using namespace vw::ip;

#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/path.hpp>
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
    composite.insert( pixel_cast_rescale<PixelRGB<uint8> >(resample(apply_mask(normalize(create_mask(DiskImageView<PixelGray<float> >(*irsrc1),
                                                                                                     irsrc1->nodata_read()))), 
                                                                    sub_scale)),
                      0, 0 );
  } else {
    composite.insert( pixel_cast_rescale<PixelRGB<uint8> >(resample(normalize(DiskImageView<PixelGray<float> >(*irsrc1)), 
                                                                    sub_scale)),
                      0, 0 );
  }

  if ( irsrc2->has_nodata_read() ) {
    composite.insert(pixel_cast_rescale<PixelRGB<uint8> >(resample(apply_mask(normalize(create_mask(DiskImageView<PixelGray<float> >(*irsrc2),
                                                                                                    irsrc2->nodata_read()))), 
                                                                   sub_scale)),
                     int32(irsrc1->format().cols * sub_scale), 0 );
  } else {
    composite.insert(pixel_cast_rescale<PixelRGB<uint8> >(resample(normalize(DiskImageView<PixelGray<float> >(*irsrc2)), 
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

    for (float r=0; r<1.0; r+=inc_amt ){

      int i = (int)(0.5 + start.x() + r*(end.x()-start.x()));
      int j = (int)(0.5 + start.y() + r*(end.y()-start.y()));
      if (i >=0 && j >=0 && i < comp.cols() && j < comp.rows())
        comp(i,j) = line_color;
    }
  }

  boost::scoped_ptr<vw::DiskImageResource> rsrc( DiskImageResource::create(out_file_name, comp.format()) );
  block_write_image( *rsrc, comp, TerminalProgressCallback( "tools.ipmatch", "Writing Debug:" ) );
}

int main(int argc, char** argv) {
  std::vector<std::string> input_file_names;
  double      matcher_threshold;
  std::string ransac_constraint, distance_metric_in, output_prefix;
  float       inlier_threshold;
  int         ransac_iterations;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h",              "Display this help message")
    ("output-prefix,o",     po::value(&output_prefix)->default_value(""), 
                            "Write output files using this prefix.")
    ("matcher-threshold,t", po::value(&matcher_threshold)->default_value(0.6), 
                            "Threshold for the separation between closest and next closest interest points.")
    ("non-kdtree",          "Use an implementation of the interest matcher that is not reliant on a KDTree algorithm")
    ("distance-metric,m",   po::value(&distance_metric_in)->default_value("L2"), 
                            "Distance metric to use.  Choose one of: [L2 (default), Hamming (only for binary types like ORB)].")
    ("ransac-constraint,r", po::value(&ransac_constraint)->default_value("similarity"), 
                            "RANSAC constraint type.  Choose one of: [similarity, homography, fundamental, or none].")
    ("inlier-threshold,i",  po::value(&inlier_threshold)->default_value(10), 
                            "RANSAC inlier threshold.")
    ("ransac-iterations",   po::value(&ransac_iterations)->default_value(100), 
                            "Number of RANSAC iterations.")
    ("debug-image,d",       "Write out debug images.");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-files", po::value<std::vector<std::string> >(&input_file_names));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-files", -1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <image file> <IP file> ..." << std::endl << std::endl;
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
    vw_out() << usage.str();
    return 1;
  }

  std::string distance_metric = distance_metric_in;
  boost::algorithm::to_lower(distance_metric);
  if ((distance_metric != "l2") && (distance_metric != "hamming")) {
    vw_out() << "Error: Did not recognize the distance metric." << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }    

  if ((input_file_names.size() < 4) || (input_file_names.size() % 2 != 0)){
    vw_out() << "Error: Must specify at least pairs of input files (image and .vwip for each)." 
             << std::endl << std::endl;
    vw_out() << usage.str();
    return 1;
  }
  // Split up the image and IP paths into two vectors
  const size_t num_input_images = input_file_names.size() / 2;
  std::vector<std::string> image_paths(num_input_images),
                           vwip_paths (num_input_images);
  for (size_t i=0; i<num_input_images; ++i) {
    image_paths[i] = input_file_names[2*i  ];
    vwip_paths [i] = input_file_names[2*i+1];
  }

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
      if ( !vm.count("non-kdtree") ) {
        // Run interest point matcher that uses KDTree algorithm.
        if (distance_metric == "l2") {
          InterestPointMatcher< L2NormMetric, NullConstraint> matcher(matcher_threshold);
          matcher(ip1, ip2, matched_ip1, matched_ip2, TerminalProgressCallback( "tools.ipmatch","Matching:"));
        } 
        if (distance_metric == "hamming") {
          InterestPointMatcher< HammingMetric, NullConstraint> matcher(matcher_threshold);
          matcher(ip1, ip2, matched_ip1, matched_ip2, TerminalProgressCallback( "tools.ipmatch","Matching:"));
        }
      } else {
        // Run interest point matcher that does not use KDTree algorithm.
        if (distance_metric == "l2") {
          InterestPointMatcherSimple<L2NormMetric, NullConstraint> matcher(matcher_threshold);
          matcher(ip1, ip2, matched_ip1, matched_ip2, TerminalProgressCallback( "tools.ipmatch","Matching:"));
        }
        if (distance_metric == "hamming") {
          InterestPointMatcherSimple<HammingMetric, NullConstraint> matcher(matcher_threshold);
          matcher(ip1, ip2, matched_ip1, matched_ip2, TerminalProgressCallback( "tools.ipmatch","Matching:"));
        }
      } // End non-KDTree case

      vw_out() << "Found " << matched_ip1.size() << " putative matches before duplicate removal.\n";

      remove_duplicates(matched_ip1, matched_ip2);
      vw_out() << "Found " << matched_ip1.size() << " putative matches.\n";

      std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matched_ip1),
                           ransac_ip2 = iplist_to_vectorlist(matched_ip2);
      std::vector<size_t> indices;
      try {
        // RANSAC is used to fit a transform between the matched sets
        // of points.  Points that don't meet this geometric
        // contstraint are rejected as outliers.
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
        continue;
      }
      vw_out() << "Found " << indices.size() << " final matches.\n";

      std::vector<InterestPoint> final_ip1, final_ip2;
      BOOST_FOREACH( size_t& index, indices ) {
        final_ip1.push_back(matched_ip1[index]);
        final_ip2.push_back(matched_ip2[index]);
      }

      if (output_prefix == "") {
        output_prefix = fs::path(image_paths[i]).replace_extension().string() + "__" +
                        fs::path(image_paths[j]).stem().string();
      }
      vw_out() << "Writing match file: " << output_prefix+".match" << std::endl;
      write_binary_match_file(output_prefix+".match", final_ip1, final_ip2);

      if (vm.count("debug-image")) {
        vw_out() << "Writing debug image: " << output_prefix+".tif" << std::endl;
        write_match_image(output_prefix+".tif",
                          image_paths[i], image_paths[j],
                          final_ip1, final_ip2);
      }
    }
  }

  return 0;
}
