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
#include <vw/Camera/EXIF.h>
#include <iostream>

using namespace vw::camera;

int main(int argc, char** argv) {

  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " <filename>.\n\n";
    std::cout << "Prints ouf the basic metadata in a TIFF or JPEG file with EXIF information.\n";
    exit(0);
  }

  ExifView exif;
  if (!exif.load_exif(argv[1])) {
    printf("Failed to load EXIF data from %s\n", argv[1]);
  }

  printf("Camera make: %s\n", exif.get_make().c_str());
  printf("Camera model: %s\n", exif.get_model().c_str());
  printf("F number: f/%f\n", exif.get_f_number());
  printf("Exposure time: %f s\n", exif.get_exposure_time());
  printf("ISO equivalent: %i\n", exif.get_iso());
  printf("ApertureValue: %f\n", exif.get_aperture_value());
  printf("ShutterSpeedValue: %f\n", exif.get_shutter_speed_value());
  printf("BrightnessValue: %f\n", exif.get_brightness_value());

  return 0;
}
