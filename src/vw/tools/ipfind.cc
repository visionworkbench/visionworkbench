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
/// Finds the interest points in an image and outputs them an 
/// ASCII format.
///
/// Example usage:
/// ipfind book.jpg book.ip
///

#include <vw/InterestPoint.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <stdio.h>

using namespace vw;
using namespace vw::ip;

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: ipfind [image] [output file]\n";
    return -1;
  }

  // Set Vision Workbench debug output level
  set_debug_level(VerboseDebugMessage + 1);

  DiskImageView<PixelGray<float> > disk_image(argv[1]);
  ImageView<float> image = channels_to_planes(disk_image);

  vw_out(InfoMessage) << "Finding interest points...\n";

  static const int MAX_INTERESTPOINT_IMAGE_DIMENSION = 2048;
  // FIXME: change to a better descriptor, like PCA
  ScaledInterestPointDetector<LogInterestOperator> detector;
  InterestPointList ip = detector(image, MAX_INTERESTPOINT_IMAGE_DIMENSION);

  unsigned num_pts = ip.size();
  if (num_pts == 0) {
    vw_out(InfoMessage) << "0 interest points found.\n";
    return 0;
  }
  unsigned size = ip.front().descriptor.size();

  // Write out detected interest points to file.
  FILE *out = fopen(argv[2], "w");
  fprintf(out, "%u %u\n", num_pts, size);
  for (InterestPointList::iterator i = ip.begin(); i != ip.end(); ++i) {
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

  vw_out(InfoMessage) << ip.size() << " interest points found.\n";

  return 0;
}
