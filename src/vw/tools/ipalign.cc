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

/// \file ipmatch.cc
///
/// Example program demonstrating how to align two images using the
/// Interest Point module.
///
/// Usage:
/// ./ipmatch [image 1] [image 2] [output prefix]

#include <vw/InterestPoint.h>
#include <vw/Math.h>
#include <vw/Image.h>
#include <vw/Mosaic.h>
#include <vw/FileIO.h>

using namespace vw;
using namespace vw::ip;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#define MAX_POINTS_TO_DRAW 1000

// ------------------------------------------------------------------

static std::string prefix_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1) 
    result.erase(index, result.size());
  return result;
}

/// Erases a file suffix if one exists and returns the base string
static std::string suffix_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1) 
    result.erase(0, index);
  return result;
}

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

// Draw the two images side by side with matching interest points
// shown with lines.
template <class ViewT>
void write_match_image(std::string out_file_name, 
                       ImageViewBase<ViewT> const& src1,
                       ImageViewBase<ViewT> const& src2,
                       std::vector<InterestPoint> matched_ip1, 
                       std::vector<InterestPoint> matched_ip2) {
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
      // red, green, blue
      comp(i,j) = PixelRGB<uint8>(255, 0, 0);
    }
  }
  
  write_image(out_file_name, comp);
}

// --------------------------------------------------------------------------------------

int main(int argc, char** argv) {

  std::vector<std::string> input_file_names;
  std::string output_file_name, interest_operator;
  float matcher_threshold, log_threshold, harris_threshold;
  int tile_size;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help", "Display this help message")
    ("verbose", "Verbose output")
    ("debug-images,d", "Produce additional debugging images as well as the aligned image.")
    ("tile-size,t", po::value<int>(&tile_size)->default_value(2048), "Specify the tile size for detecting interest points.")
    ("output-file,o", po::value<std::string>(&output_file_name)->default_value("aligned-image.jpg"), "Specify the output filename")


    // Interest point detector options
    ("matcher-threshold,t", po::value<float>(&matcher_threshold)->default_value(0.5), "Rejects points during matching if best > matcher_threshold * second_best")
    ("log-threshold", po::value<float>(&log_threshold)->default_value(0.03), "Sets the threshold for the Laplacian of Gaussian Interest Operator")
    ("harris-threshold", po::value<float>(&harris_threshold)->default_value(1e-5), "Sets the threshold for the Harris Interest Operator")
    ("interest-operator", po::value<std::string>(&interest_operator)->default_value("LoG"), "Choose an interest metric from [LoG, Harris]")
    ("single-scale", "Do not use the scale-space interest point detector.")

    // Alignment options
    ("homography", "Align images using a full projective transform (homography).  By default, aligment uses a more restricted Similarity transform.");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value<std::vector<std::string> >(&input_file_names));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", -1);

  po::variables_map vm;
  po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
  po::notify( vm );

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <filename>..." << std::endl << std::endl;
  usage << general_options << std::endl;

  if( vm.count("help") ) {
    std::cout << usage.str();
    return 1;
  }

  if( input_file_names.size() != 2 ) {
    std::cout << "Error: Must specify exactly two input files!" << std::endl << std::endl;
    std::cout << usage.str();
    return 1;
  }

  // Set Vision Workbench debug output level
  if ( vm.count("verbose") != 0 )
    set_debug_level(DebugMessage);
  else 
    set_debug_level(InfoMessage);

  // Load the two images
  DiskImageView<PixelRGB<uint8> > left_image(input_file_names[0]);
  DiskImageView<PixelRGB<uint8> > right_image(input_file_names[1]);

  // Image Alignment
  //
  // Images are aligned by computing interest points, matching
  // them using a standard 2-Norm nearest-neighor metric, and then
  // rejecting outliers by fitting a similarity between the
  // putative matches using RANSAC.
  
  vw_out(InfoMessage) << "\nInterest Point Detection:\n";
  InterestPointList ip1, ip2;
  if (interest_operator == "Harris") {
    HarrisInterestOperator interest_operator(harris_threshold);
    if (!vm.count("single-scale")) {
      ScaledInterestPointDetector<HarrisInterestOperator> detector(interest_operator);
      ip1 = detector(left_image, tile_size);
      ip2 = detector(right_image, tile_size);
    } else {
      InterestPointDetector<HarrisInterestOperator> detector(interest_operator);
      ip1 = detector(left_image, tile_size);
      ip2 = detector(right_image, tile_size); 
    }
  } else if (interest_operator == "LoG") {
    // Use a scale-space Laplacian of Gaussian feature detector. The
    // associated threshold is abs(interest) > interest_threshold.
    LogInterestOperator interest_operator(log_threshold);
    if (!vm.count("single-scale")) {
      ScaledInterestPointDetector<LogInterestOperator> detector(interest_operator);
      ip1 = detector(left_image, tile_size);
      ip2 = detector(right_image, tile_size); 
    } else {
      InterestPointDetector<LogInterestOperator> detector(interest_operator);
      ip1 = detector(left_image, tile_size);
      ip2 = detector(right_image, tile_size); 
    }
  } else {
    std::cout << "Unknown interest operator: " << interest_operator << ".  Options are : [ Harris, LoG ]\n";
    exit(0);
  }

  // Write out images with interest points marked
  std::string prefix = prefix_from_filename(output_file_name);
  std::string suffix = suffix_from_filename(output_file_name);
  if (vm.count("debug-images")) {
    write_point_image(prefix + "-PL" + suffix, left_image, ip1);
    write_point_image(prefix + "-PR" + suffix, right_image, ip2);
  }

  // Generate descriptors for interest points.
  vw_out(InfoMessage) << "Generating descriptors... ";
  PatchDescriptorGenerator descriptor;
  //PCASIFTDescriptorGenerator descriptor("pca_basis.exr", "pca_avg.exr");
  descriptor(left_image, ip1);
  descriptor(right_image, ip2);
  vw_out(InfoMessage) << "done.\n";
  
  // The basic interest point matcher does not impose any
  // constraints on the matched interest points.
  vw_out(InfoMessage) << "\nInterest Point Matching:\n";
  DefaultMatcher matcher(matcher_threshold);

  // RANSAC needs the matches as a vector, and so does the matcher.
  // this is messy, but for now we simply make a copy.
  std::vector<InterestPoint> ip1_copy(ip1.size()), ip2_copy(ip2.size());
  std::copy(ip1.begin(), ip1.end(), ip1_copy.begin());
  std::copy(ip2.begin(), ip2.end(), ip2_copy.begin());

  std::vector<InterestPoint> matched_ip1, matched_ip2;
  matcher(ip1_copy, ip2_copy, matched_ip1, matched_ip2);
  vw_out(InfoMessage) << "\tFound " << matched_ip1.size() << " putative matches.\n";

  // Write out the putative point correspondence image
  if (vm.count("debug-images"))
    write_match_image(prefix+"-putative-match" + suffix, left_image, right_image, matched_ip1, matched_ip2);
  
  // RANSAC is used to fit a similarity transform between the
  // matched sets of points
  Matrix<double> align_matrix;
  if (vm.count("homography")) 
    align_matrix = ransac(matched_ip2, matched_ip1, 
                          vw::math::HomographyFittingFunctor(),
                          InterestPointErrorMetric());
  else
    align_matrix = ransac(matched_ip2, matched_ip1, 
                          vw::math::AffineFittingFunctor(),
                          InterestPointErrorMetric());

  vw_out(InfoMessage) << "Writing out aligned pair of images\n";

  // Write out the aligned pair of images
  ImageViewRef<PixelRGB<uint8> > aligned_image = transform(right_image, HomographyTransform(align_matrix),
                                                           left_image.cols(), left_image.rows());
  write_image(output_file_name, aligned_image); 

  return 0;
}
