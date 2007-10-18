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

/// \file match_pair.cc
///
/// Example program demonstrating how to align two images using the
/// Interest Point module.
///
/// Usage:
/// ./match_pair [image 1] [image 2] [output prefix]

#include <vw/InterestPoint.h>
#include <vw/Image.h>
#include <vw/Math.h>
#include <vw/Mosaic.h>
#include <vw/FileIO.h>

using namespace vw;
using namespace vw::ip;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#define TOP_POINTS 1000


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

// Draws the interest points as green squares, with orientation and
// scale indicated by a red line.
int draw_points( ImageView<PixelRGB<float> >& image,
                 KeypointList const& points ) {
  // Draw points into color planes
  int n = 0;
  for (KeypointList::const_iterator pt = points.begin();
       pt != points.end() && n < TOP_POINTS; ++pt, ++n) {
    // Draw a red line from the point outward along the orientation
    for (int r=0; r<(int)(8*(*pt).scale); ++r){
      int i = (int)(0.5 + (*pt).x + r*cos((*pt).orientation));
      int j = (int)(0.5 + (*pt).y + r*sin((*pt).orientation));
      // red, green, blue
      image(i,j)[0]=1.0; image(i,j)[1]=0.0; image(i,j)[2]=0.0;
    }
    // Draw a green 3x3 filled square at the point to indicate center
    int i0 = (int)(0.5 + (*pt).x);
    int j0 = (int)(0.5 + (*pt).y);
    for (int j=j0-1; j<=j0+1; ++j){
      for (int i=i0-1; i<=i0+1; ++i){
        // red, green, blue
        image(i,j)[0]=0.0; image(i,j)[1]=1.0; image(i,j)[2]=0.0;
      }
    }
  }
  return 0;
}

// Draw the interest points and write as an image.
template <class ViewT>
void write_point_image(std::string out_file_name, ImageViewBase<ViewT> const& src,
                       KeypointList const& points) {
  ImageView<PixelRGB<float> > viz = copy(src);
  draw_points(viz, points);
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
  mosaic::ImageComposite<typename ViewT::pixel_type> composite;
  composite.insert(src1.impl(),0,0);
  composite.insert(src2.impl(),src1.impl().cols(),0);
  composite.set_draft_mode( true );
  composite.prepare();

  // Rasterize the composite so that we can draw on it.
  ImageView<PixelRGB<float> > comp = composite;
  
  // Draw a red line between matching keypoints
  for (unsigned int i = 0; i < matched_ip1.size(); ++i) {
    Vector2 start(matched_ip1[i][0], matched_ip1[i][1]);
    Vector2 end(matched_ip2[i][0]+src1.impl().cols(), matched_ip2[i][1]);
    for (float r=0; r<1.0; r+=1/norm_2(end-start)){
      int i = (int)(0.5 + start.x() + r*(end.x()-start.x()));
      int j = (int)(0.5 + start.y() + r*(end.y()-start.y()));
      // red, green, blue
      comp(i,j) = PixelRGB<float>(1.0, 0.0, 0.0);
    }
  }
  
  write_image(out_file_name, comp);
}

template <class PixelT> 
void do_match(std::vector<std::string> input_file_names, std::string output_file_name, bool write_debug_images, bool use_homography, double matcher_threshold) {
  // Load the two images
  std::string input_file1 = input_file_names[0];
  std::string input_file2 = input_file_names[1];
  DiskImageView<PixelT> left_disk_image(input_file1);
  DiskImageView<PixelT> right_disk_image(input_file2);

  // Using DiskImageViews directly does not quite work yet.
  ImageView<PixelGray<float> > left = left_disk_image;
  ImageView<PixelGray<float> > right = right_disk_image;
  
  // Image Alignment
  //
  // Images are aligned by computing interest points, matching
  // them using a standard 2-Norm nearest-neighor metric, and then
  // rejecting outliers by fitting a similarity between the
  // putative matches using RANSAC.
  
  // Interest points are matched in image chunk of <= 2048x2048
  // pixels to conserve memory.
  vw_out(InfoMessage) << "\nInterest Point Detection:\n";

  static const int MAX_KEYPOINT_IMAGE_DIMENSION = 2048;
  typedef DefaultThresholdT<LoGInterest>::type Threshold;
  typedef ScaledInterestPointDetector<LoGInterest> Detector;

  // Use a scale-space Laplacian of Gaussian feature detector. The associated
  // default threshold is abs(interest) > 0.03.
  Detector detector;
  // If too few interest points, try something like:
  // Detector detector(Threshold(0.01));

  // KeypointList := std::list<InterestPoint>
  KeypointList ip1 = interest_points(left, detector, MAX_KEYPOINT_IMAGE_DIMENSION);
  KeypointList ip2 = interest_points(right, detector, MAX_KEYPOINT_IMAGE_DIMENSION);

  static const int NUM_POINTS = 1000;
  vw_out(InfoMessage) << "Truncating to " << NUM_POINTS << " points:\n";
  cull_interest_points(ip1, NUM_POINTS);
  cull_interest_points(ip2, NUM_POINTS);

  // Write out images with interest points marked
  std::string prefix = prefix_from_filename(output_file_name);
  std::string suffix = suffix_from_filename(output_file_name);
  if (write_debug_images) {
    write_point_image(prefix + "-PL" + suffix, left, ip1);
    write_point_image(prefix + "-PR" + suffix, right, ip2);
  }

  // Generate descriptors for interest points.
  vw_out(InfoMessage) << "Generating descriptors:\n";
  compute_descriptors(left, ip1, PatchDescriptor() );
  compute_descriptors(right, ip2, PatchDescriptor() );
    
  // The basic interest point matcher does not impose any
  // constraints on the matched interest points.
  vw_out(InfoMessage) << "\nInterest Point Matching:\n";
  DefaultMatcher matcher(matcher_threshold);
  // RANSAC needs the matches as a vector.
  std::vector<InterestPoint> matched_ip1, matched_ip2;
  matcher.match(ip1, ip2, matched_ip1, matched_ip2);
  vw_out(InfoMessage) << "\tFound " << matched_ip1.size() << " putative matches.\n";

  // Write out the putative point correspondence image
  if (write_debug_images)
    write_match_image(prefix+"-putative-match" + suffix, left, right, matched_ip1, matched_ip2);
  
  // RANSAC is used to fit a similarity transform between the
  // matched sets of points
  Matrix<double> align_matrix;
  if (use_homography) 
    align_matrix = ransac(matched_ip2, matched_ip1, 
                          vw::math::HomographyFittingFunctor(),
                          KeypointErrorMetric());
  else
    align_matrix = ransac(matched_ip2, matched_ip1, 
                          vw::math::SimilarityFittingFunctor(),
                          KeypointErrorMetric());

  // Write out the aligned pair of images
  ImageViewRef<PixelT> aligned_image = transform(right_disk_image, HomographyTransform(align_matrix),
                                                   left_disk_image.cols(), left_disk_image.rows());
  write_image(output_file_name, channel_cast_rescale<uint8>(aligned_image)); 
}

int main(int argc, char** argv) {

  std::vector<std::string> input_file_names;
  std::string output_file_name;
  double matcher_threshold;

  po::options_description general_options("Options");
  general_options.add_options()
    ("help", "Display this help message")
    ("output-file,o", po::value<std::string>(&output_file_name)->default_value("output.jpg"), "Specify the output filename")
    ("matcher-threshold,t", po::value<double>(&matcher_threshold)->default_value(0.5), "Rejects points during matching if best > matcher_threshold * second_best")
    ("debug-images,d", "Produce additional debugging images as well as the aligned image.")
    ("grayscale", "Use grayscale image processing.")
    ("homography", "Align images using a full projective transform (homography).  By default, aligment uses a more restricted Similarity transform.")
    ("verbose", "Verbose output");
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

  if ( vm.count("grayscale") != 0 )
    do_match<PixelGray<float> >(input_file_names, output_file_name, vm.count("debug-images"), vm.count("homography"), matcher_threshold);
  else 
    do_match<PixelRGB<float> >(input_file_names, output_file_name, vm.count("debug-images"), vm.count("homography"), matcher_threshold);

  return 0;
}
