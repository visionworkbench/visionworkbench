// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Camera/Exif.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace vw::camera;

int main(int argc, char** argv) {

  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " <filename>.\n\n";
    std::cout << "Prints ouf the basic metadata in a TIFF or JPEG file with EXIF information.\n";
    exit(0);
  }

  try {
    ExifView exif(argv[1]);

    printf("Camera Settings:\n");
    printf("Camera make: %s\n", exif.get_make().c_str());
    printf("Camera model: %s\n", exif.get_model().c_str());
    printf("F number: f/%f\n", exif.get_f_number());
    printf("Exposure time: %f s\n", exif.get_exposure_time());
    printf("ISO equivalent: %f\n", exif.get_iso());
    double focal_length;
    exif.query_by_tag(EXIF_FocalLength, focal_length);
    printf("Focal Length: %f mm\n", focal_length);

    printf("\nDerived Values:\n");
    printf("Average Luminance: %f cd/m^2\n", exif.get_average_luminance());

    printf("\nValues in the APEX logarithmic system:\n");
    printf("ApertureValue: %f\n", exif.get_aperture_value());
    printf("ShutterSpeedValue: %f\n", exif.get_time_value());
    printf("ExposureValue: %f\n", exif.get_exposure_value());
    printf("FilmSpeedValue: %f\n", exif.get_film_speed_value());
    printf("LuminanceValue: %f\n", exif.get_luminance_value());


  } catch (ExifErr &e) {
    std::cout << "An EXIF error occurred: " << e.what() << "\n";
  }

  return 0;
}
