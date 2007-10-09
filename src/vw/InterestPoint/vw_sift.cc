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

/// \file vw_sift.cc
///
/// Finds the SIFT keypoints in an image and outputs them in the
/// ASCII format used by David Lowe's sift binary.
///
/// Example usage:
/// vw_sift book.jpg book.key
///
/// As opposed to usage for Lowe's sift:
/// sift <book.pgm >book.key
///
/// NOTE: Currently uses scale-invariant LoG interest point detector
/// instead of actual SIFT detector

#include <vw/InterestPoint.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <stdio.h>

using namespace vw;
using namespace vw::ip;

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: vw_sift [image] [output file]\n";
    return -1;
  }

  // Set Vision Workbench debug output level
  set_debug_level(VerboseDebugMessage + 1);

  DiskImageView<PixelGray<float> > disk_image(argv[1]);
  ImageView<float> image = channels_to_planes(disk_image);

  vw_out(InfoMessage) << "Finding keypoints...\n";

  static const int MAX_KEYPOINT_IMAGE_DIMENSION = 2048;
  // FIXME: change to SIFT_Descriptor when able
  typedef ScaledInterestPointDetector<LoGInterest> Detector;

  KeypointList ip = interest_points(image, Detector(), MAX_KEYPOINT_IMAGE_DIMENSION);
  unsigned num_pts = ip.size();
  if (num_pts == 0) {
    vw_out(InfoMessage) << "0 keypoints found.\n";
    return 0;
  }
  unsigned size = ip.front().descriptor.size();

  // Write out detected keypoints to file.
  FILE *out = fopen(argv[2], "w");
  fprintf(out, "%u %u\n", num_pts, size);
  for (KeypointList::iterator i = ip.begin(); i != ip.end(); ++i) {
    float orientation = i->orientation;
    while (orientation > M_PI) orientation -= 2 * M_PI;
    while (orientation < -M_PI) orientation += 2 * M_PI;
    fprintf(out, "%.2f %.2f %.2f %.3f", i->y, i->x, i->scale, orientation);
    for (unsigned element = 0; element < size; ++element) {
      if (element % 20 == 0) fprintf(out, "\n");
      fprintf(out, " %u", (unsigned)(i->descriptor[element] * 255.0));
    }
    fprintf(out, "\n");
  }
  fclose(out);

  vw_out(InfoMessage) << ip.size() << " keypoints found.\n";

  return 0;
}
