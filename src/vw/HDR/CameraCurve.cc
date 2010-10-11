// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file CameraCurve.cc
///
/// Functions for deducing the camera response curve by comparing
/// relative brightness values from several images of the same scene.

#include <vw/Camera/Exif.h>
#include <vw/HDR/CameraCurve.h>
#include <vw/Math/LinearAlgebra.h>
#include <boost/algorithm/string.hpp>

namespace vw {
namespace hdr {

  /// Generate a set of brightness values based on a known ratio of
  /// exposures between images.  Note that this brightness value will
  /// only be correct in a relative sense.  Given this, we choose
  /// (arbitrarily) a set of brightness values that correspond to a
  /// "typical" exposure of 1/60th of a second at ISO 100 and an
  /// aperture of f/5.6 (B = 235.2 in this case)
  ///
  ///
  std::vector<double> brightness_values_from_exposure_ratio(double exposure_ratio, int size) {
    const double base_brightness = 235.2; // average luminance for
                                          // 1/60th s exposure at ISO
                                          // 100 and f/5.6
    std::vector<double> brightness_values(size);
    for ( unsigned i=0; i < brightness_values.size(); ++i ) {
      brightness_values[i] = base_brightness*pow(exposure_ratio, (double)i);
    }
    return brightness_values;
  }

  /// Generate a list of brightness values from EXIF information
  /// stored in a list of files.
  std::vector<double> brightness_values_from_exif(std::vector<std::string> const& filenames) {
    int num_images = filenames.size();
    std::vector<double> brightness_values(num_images);

    for (int i = 0; i < num_images; i++) {
      vw::camera::ExifView exif(filenames[i].c_str());
      brightness_values[i] = exif.get_average_luminance();
    }
    return brightness_values;
  }

  // A simple gaussian weighting function for use in camera curve
  // estimation (below)
  inline double gaussian_weighting_func(double x) {
    return exp(-pow((x-0.5),2)/(0.07));
  }

  Vector<double> estimate_camera_curve(vw::Matrix<double> const& pixels,
                                       std::vector<double> const& brightness_values) {

    const int n = 256;                    // Create a solution with 256 points
    const int smoothing_factor = 10;      // Smoothing by a factor of 10

    // Initialize storage
    Matrix<double> A(pixels.rows()*pixels.cols()+n+1, n+pixels.rows());
    Vector<double> b(A.rows());

    // Include data fitting equations
    int k = 0;
    for (unsigned i = 0; i < pixels.rows(); ++i) {
      for (unsigned j = 0; j < pixels.cols(); ++j) {
        int idx = int(pixels(i,j)*(n-1));
        double wij = gaussian_weighting_func(pixels(i,j));
        A(k,idx) = wij;
        A(k,n+i) = -wij;
        b(k) = wij * log(1/brightness_values[j]);
        ++k;
      }
    }

    // Fix the curve by setting its middle value to 0
    A(k,n/2) = 1;
    ++k;

    // Include the smoothness equations
    for (int i = 0; i < n-2; ++i) {
      double weight = gaussian_weighting_func(double(i+1)/(n-1));
      A(k,i) = smoothing_factor*weight;
      A(k,i+1) = -2*smoothing_factor*weight;
      A(k,i+2) = smoothing_factor*weight;
      ++k;
    }

    // Solve the system using linear least squares
    Vector<double> x = least_squares(A, b);
    return subvector(x,0,n);
  }

  void write_curves(std::string const& curves_file,
                    CameraCurveFn const &curves) {

    FILE* output_file = fopen(curves_file.c_str(), "w");
    if ( !output_file ) vw_throw( IOErr() << "write_curves: failed to open file for writing." );
    for (unsigned i = 0; i < curves.num_channels(); ++i) {
      for ( unsigned j = 0; j < curves.lookup_table(0).size(); ++j ) {
        fprintf(output_file, "%f ", curves.lookup_table(i)[j]);
      }
      fprintf(output_file, "\n");
    }
    fclose(output_file);
  }

  CameraCurveFn read_curves(std::string const& curves_file) {
    FILE* input_file = fopen(curves_file.c_str(), "r");
    if ( !input_file ) vw_throw( IOErr() << "read_curves: failed to open file for reading." );

    char c_line[10000];

    std::vector<vw::Vector<double> > lookup_tables;
    while ( !feof(input_file) ) {
      if ( !fgets(c_line, 10000, input_file) )
        break;
      std::string line = c_line;
      boost::trim_left(line);
      boost::trim_right(line);

      std::vector< std::string > split_vec; // #2: Search for individual values
      boost::split( split_vec, line, boost::is_any_of(" ") );
      Vector<double> curve(split_vec.size());
      for ( unsigned i = 0; i < split_vec.size(); ++i ) {
        curve[i] = atof(split_vec[i].c_str());
      }
      lookup_tables.push_back(curve);
    }
    fclose(input_file);

    return CameraCurveFn(lookup_tables);
  }


}} // namespace vw::hdr
