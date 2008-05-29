// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file ipfind.cc
///
/// Finds the interest points in an image and outputs them an Binary
/// (default) or ASCII format.  The ASCII format is compatible with
/// the popular Lowe-SIFT toolchain.
///
#include <vw/InterestPoint.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Math.h>
#include <vw/Mosaic.h>

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

// Duplicate matches for any given interest point probably indicate a
// poor match, so we cull those out here.
void remove_duplicates(std::vector<InterestPoint> &ip1, std::vector<InterestPoint> &ip2) {
  std::vector<InterestPoint> new_ip1, new_ip2;
  
  for (unsigned i = 0; i < ip1.size(); ++i) {
    bool bad_entry = false;
    for (unsigned j = 0; j < ip1.size(); ++j) {
      if (i != j && 
          ((ip1[i].x == ip1[j].x && ip1[i].y == ip1[j].y) || 
           (ip2[i].x == ip2[j].x && ip2[i].y == ip2[j].y)) ) {
        bad_entry = true;
      }
    }
    if (!bad_entry) {
      new_ip1.push_back(ip1[i]);
      new_ip2.push_back(ip2[i]);
    }
  }
  
  ip1 = new_ip1;
  ip2 = new_ip2;
}

// Draw the two images side by side with matching interest points
// shown with lines.
static void write_match_image(std::string out_file_name, 
                              std::string const& file1,
                              std::string const& file2,
                              std::vector<InterestPoint> matched_ip1, 
                              std::vector<InterestPoint> matched_ip2) {
  // Skip image pairs with no matches.
  if (matched_ip1.size() == 0)
    return;

  DiskImageView<PixelRGB<uint8> > src1(file1);
  DiskImageView<PixelRGB<uint8> > src2(file2);

  mosaic::ImageComposite<PixelRGB<uint8> > composite;
  composite.insert(pixel_cast<PixelRGB<uint8> >(channel_cast_rescale<uint8>(src1.impl())),0,0);
  composite.insert(pixel_cast<PixelRGB<uint8> >(channel_cast_rescale<uint8>(src2.impl())),src1.impl().cols(),0);
  composite.set_draft_mode( true );
  composite.prepare();

  // Rasterize the composite so that we can draw on it.
  ImageView<PixelRGB<uint8> > comp = composite;
  
  // Draw a red line between matching interest points
  for (unsigned int i = 0; i < matched_ip1.size(); ++i) {
    Vector2 start(matched_ip1[i].x, matched_ip1[i].y);
    Vector2 end(matched_ip2[i].x+src1.impl().cols(), matched_ip2[i].y);
    for (float r=0; r<1.0; r+=1/norm_2(end-start)){
      int i = (int)(0.5 + start.x() + r*(end.x()-start.x()));
      int j = (int)(0.5 + start.y() + r*(end.y()-start.y()));
      if (i >=0 && j >=0 && i < comp.cols() && j < comp.rows())
        comp(i,j) = PixelRGB<uint8>(255, 0, 0);
    }
  }
  
  write_image(out_file_name, comp);
}

int main(int argc, char** argv) {
  std::vector<std::string> input_file_names;
  double matcher_threshold;
  std::string ransac_constraint;
  int inlier_threshold;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help", "Display this help message")
    ("matcher-threshold,t", po::value<double>(&matcher_threshold)->default_value(0.8), "Threshold for the interest point matcher.")
    ("ransac-constraint,r", po::value<std::string>(&ransac_constraint)->default_value("similarity"), "RANSAC constraint type.  Choose one of: [similarity, homograhy].")
    ("inlier-threshold,i", po::value<int>(&inlier_threshold)->default_value(10), "RANSAC inlier threshold.")
    ("debug-image,d", "Write out debug images.");

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

  if( input_file_names.size() < 2 ) {
    vw_out(0) << "Error: Must specify at least two input files!" << std::endl << std::endl;
    vw_out(0) << usage.str();
    return 1;
  }

  // Iterate over combinations of the input files and find interest points in each.
  for (unsigned i = 0; i < input_file_names.size(); ++i) {
    for (unsigned j = i+1; j < input_file_names.size(); ++j) {
      
      // Read each file off disk
      std::vector<InterestPoint> ip1, ip2;
      ip1 = read_binary_ip_file(prefix_from_filename(input_file_names[i])+".vwip");
      ip2 = read_binary_ip_file(prefix_from_filename(input_file_names[j])+".vwip");
      vw_out(0) << "Matching between " << input_file_names[i] << " (" << ip1.size() << " points) and " << input_file_names[j] << " (" << ip2.size() << " points).\n"; 

      // Run the basic interest point matcher.
      InterestPointMatcher<L2NormMetric,NullConstraint> matcher(matcher_threshold);
      //      DefaultMatcher matcher(matcher_threshold);
      std::vector<InterestPoint> matched_ip1, matched_ip2;
      matcher(ip1, ip2, matched_ip1, matched_ip2, false, TerminalProgressCallback());
      remove_duplicates(matched_ip1, matched_ip2);
      vw_out(InfoMessage) << "Found " << matched_ip1.size() << " putative matches.\n";

      std::vector<Vector3> ransac_ip1 = iplist_to_vectorlist(matched_ip1);
      std::vector<Vector3> ransac_ip2 = iplist_to_vectorlist(matched_ip2);
      std::vector<int> indices;
      try {
        // RANSAC is used to fit a transform between the matched sets
        // of points.  Points that don't meet this geometric
        // contstraint are rejected as outliers.
        if (ransac_constraint == "similarity") {
          RandomSampleConsensus<math::SimilarityFittingFunctor, InterestPointErrorMetric> ransac( vw::math::SimilarityFittingFunctor(),
                                                                                                  InterestPointErrorMetric(), 
                                                                                                  inlier_threshold ); // inlier_threshold
          Matrix<double> H = ransac(ransac_ip1,ransac_ip2);
          indices = ransac.inlier_indices(H,ransac_ip1,ransac_ip2);
        } else if (ransac_constraint == "homography") {
          RandomSampleConsensus<math::HomographyFittingFunctor, InterestPointErrorMetric> ransac( vw::math::HomographyFittingFunctor(),
                                                                                                  InterestPointErrorMetric(), 
                                                                                                  inlier_threshold ); // inlier_threshold
          Matrix<double> H = ransac(ransac_ip1,ransac_ip2);
          indices = ransac.inlier_indices(H,ransac_ip1,ransac_ip2);
        } else {
          std::cout << "Unknown RANSAC constraint type: " << ransac_constraint << ".  Choose one of: [similarity, homography]\n";
          exit(0);
        }
      } catch (vw::ip::RANSACErr &e) {
        std::cout << "RANSAC Failed: " << e.what() << "\n";
        exit(0);
      }
      vw_out(InfoMessage) << "Found " << indices.size() << " final matches.\n";
      
      std::vector<InterestPoint> final_ip1, final_ip2;
      for (unsigned idx=0; idx < indices.size(); ++idx) {
        final_ip1.push_back(matched_ip1[indices[idx]]);
        final_ip2.push_back(matched_ip2[indices[idx]]);
      }

      std::string output_filename = 
        prefix_from_filename(input_file_names[i]) + "__" +
        prefix_from_filename(input_file_names[j]) + ".match";
      write_binary_match_file(output_filename, final_ip1, final_ip2);

      if (vm.count("debug-image")) {
        std::string matchimage_filename = 
          prefix_from_filename(input_file_names[i]) + "__" +
          prefix_from_filename(input_file_names[j]) + ".jpg";
        write_match_image(matchimage_filename, 
                          input_file_names[i], input_file_names[j],
                          final_ip1, final_ip2);
      }
    }
  }
}   

