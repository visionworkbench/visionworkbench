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

#include <vw/Camera/Exif.h>

#include <math.h>

// --------------------------------------------------------------
//                   ExifView
// --------------------------------------------------------------

vw::camera::ExifView::ExifView(std::string const& filename) {
  if (!m_data.import_data(filename)) {
    vw_throw(ExifErr() << "Could not parse EXIF data out of \"" << filename << "\".");
  }
}

// Query the data by tag ID (common tags are enumerated at the top if
// Exif.h)
void vw::camera::ExifView::query_by_tag(const uint16 tag, int& value) {
  bool success = m_data.get_tag_value(tag, value);
  if (!success) vw_throw(ExifErr() << "Could not read EXIF tag " << tag << ".");
}

// Query the data by tag ID (common tags are enumerated at the top if
// Exif.h)
void vw::camera::ExifView::query_by_tag(const uint16 tag, double& value) {
  bool success = m_data.get_tag_value(tag, value);
  if (!success) vw_throw(ExifErr() << "Could not read EXIF tag: " << tag << ".");
}

// Query the data by tag ID (common tags are enumerated at the top if
// Exif.h)
void vw::camera::ExifView::query_by_tag(const uint16 tag, std::string& value) {
  bool success = m_data.get_tag_value(tag, value);
  if (!success) vw_throw(ExifErr() << "Could not read EXIF tag: " << tag << ".");
}

// Camera info
std::string vw::camera::ExifView::get_make() {
  std::string make;
  query_by_tag(EXIF_Make, make);
  return make;
}

std::string vw::camera::ExifView::get_model() {
  std::string model;
  query_by_tag(EXIF_Model, model);
  return model;
}

// Camera settings
double vw::camera::ExifView::get_f_number() {
  double value;
  try { 
    query_by_tag(EXIF_FNumber, value); 
    return value;
  } catch (ExifErr &e) {
    query_by_tag(EXIF_ApertureValue, value); 
    return pow(2.0, value * 0.5);
  }
}

double vw::camera::ExifView::get_exposure_time() {
  double value;
  try { 
    query_by_tag(EXIF_ExposureTime, value); 
    return value;
  } catch (ExifErr &e) {
    query_by_tag(EXIF_ShutterSpeedValue, value); 
    return pow(2.0, -value);
  }
}

double vw::camera::ExifView::get_iso() {
  double value;
  try { 
    query_by_tag(EXIF_ISOSpeedRatings, value); 
    return value;
  } catch (ExifErr &e) {
    query_by_tag(EXIF_ExposureIndex, value); 
    return value;
  }
  // otherwise probably have to dig through MakerNote
}

// FIXME: report in some logical way when value doesn't exist
double vw::camera::ExifView::get_aperture_value() {
  double value;
  try { 
    query_by_tag(EXIF_ApertureValue, value); 
    return value;
  } catch (ExifErr &e) {
    query_by_tag(EXIF_FNumber, value); 
    return 2 * log(value)/log(2.);  // log2(value) = log(value)/log(2)
  }
}

double vw::camera::ExifView::get_time_value() {
  double value;
  try { 
    query_by_tag(EXIF_ShutterSpeedValue, value); 
    return value;
  } catch (ExifErr &e) {
    query_by_tag(EXIF_ExposureTime, value); 
    return log(1/value)/log(2.); // log2(value) = log(value)/log(2)
  }
}

double vw::camera::ExifView::get_exposure_value() {
  return get_time_value() + get_aperture_value();
}

// Film speed value is log_2(N * iso)
// 
// N is a constant that establishes the relationship between the ASA
// arithmetic film speed and the ASA speed value. The value of N
// 1/3.125, as defined by the EXIF 2.2 spec.  The constant K is the
// reflected light meter calibration constant.  See
// http://en.wikipedia.org/wiki/Light_meter#Exposure_meter_calibration
//
double vw::camera::ExifView::get_film_speed_value() {
  double iso = (double)get_iso();
  const double N = 1/3.125;
  //  const double K = 12.5; 
  return log(iso * N)/log(2.); // log2(value) = log(value)/log(2)
}

double vw::camera::ExifView::get_luminance_value() {
  double Bv;
  try{
    query_by_tag(EXIF_BrightnessValue, Bv);
    return Bv;
  } catch (ExifErr &e) { 
    try {
      double Av = get_aperture_value();
      double Tv = get_time_value();
      double Sv = get_film_speed_value();
      return (Av + Tv - Sv);
    } catch (ExifErr &e) {
      vw_throw(ExifErr() << "Insufficient EXIF information to compute brightness value.");
      return 0; // never reached
    }
  }
}


// Film speed value is log_2(N * iso)
// 
// N is a constant that establishes the relationship between the ASA
// arithmetic film speed and the ASA speed value. The value of N
// 1/3.125, as defined by the EXIF 2.2 spec.  The constant K is the
// reflected light meter calibration constant.  See
// http://en.wikipedia.org/wiki/Light_meter#Exposure_meter_calibration
//
double vw::camera::ExifView::get_average_luminance() {
  // const double N = 1/3.125;
  const double K = 12.5; 

  try {
    double A = get_f_number();
    double T = get_exposure_time();
    double S = get_iso();
    double B = (A*A*K)/(T*S);
    return B;
  } catch (ExifErr &e) {
    vw_throw(ExifErr() << "Insufficient EXIF information to compute average scene luminance.");
    return 0; // never reached
  }
}

