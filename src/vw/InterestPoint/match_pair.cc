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
#include <vw/FileIO.h>

using namespace vw;
using namespace vw::ip;

#define TOP_POINTS 1000

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

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "Usage: match_pair [image 1] [image 2] [output prefix]\n";
    return -1;
  }

  // Set Vision Workbench debug output level
  set_debug_level(DebugMessage);

  // Load the two images
  std::string input_file1 = argv[1];
  std::string input_file2 = argv[2];
  std::string out_prefix = argv[3];
  DiskImageView<PixelGray<float> > left_disk_image(input_file1);
  DiskImageView<PixelGray<float> > right_disk_image(input_file2);

  // Using DiskImageViews directly does not quite work yet.
  ImageView<float> left = channels_to_planes(left_disk_image);
  ImageView<float> right = channels_to_planes(right_disk_image);
  
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
  write_point_image(out_prefix + "-PL.jpg", left_disk_image, ip1);
  write_point_image(out_prefix + "-PR.jpg", right_disk_image, ip2);

  // Generate descriptors for interest points.
  vw_out(InfoMessage) << "Generating descriptors:\n";
  PatchDescriptor<float> desc1(left);
  PatchDescriptor<float> desc2(right);
  desc1.compute_descriptors(ip1);
  desc2.compute_descriptors(ip2);
    
  // The basic interest point matcher does not impose any
  // constraints on the matched interest points.
  vw_out(InfoMessage) << "\nInterest Point Matching:\n";
  DefaultMatcher matcher;
  // RANSAC needs the matches as a vector.
  std::vector<InterestPoint> matched_ip1, matched_ip2;
  matcher.match(ip1, ip2, matched_ip1, matched_ip2);
  vw_out(InfoMessage) << "\tFound " << matched_ip1.size() << " putative matches.\n";
  
  // RANSAC is used to fit a similarity transform between the
  // matched sets of points
  Matrix<double> align_matrix = ransac(matched_ip2, matched_ip1, 
                                       vw::math::SimilarityFittingFunctor(),
                                       KeypointErrorMetric());

  // Write out the aligned pair of images
  ImageViewRef<PixelGray<float> > Limg = left_disk_image;
  ImageViewRef<PixelGray<float> > Rimg = transform(right_disk_image, HomographyTransform(align_matrix),
                                                   left_disk_image.cols(), left_disk_image.rows());

  std::string output_file1 = out_prefix + "-L.jpg";
  std::string output_file2 = out_prefix + "-R.jpg";

  write_image(output_file1, channel_cast_rescale<uint8>(Limg));
  write_image(output_file2, channel_cast_rescale<uint8>(Rimg)); 
  
  return 0;
}
