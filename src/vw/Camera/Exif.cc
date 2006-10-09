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

bool vw::camera::ExifView::load_exif(std::string const& filename) {
  return data.import_data(filename.c_str());
}

// Camera info
std::string vw::camera::ExifView::get_make() {
  char * make;
  if (data.get_tag_value(TAG_Make, make)) return make;
  return std::string();
}

std::string vw::camera::ExifView::get_model() {
  char * model;
  if (data.get_tag_value(TAG_Model, model)) return model;
  return std::string();
}

// Camera settings
double vw::camera::ExifView::get_f_number() {
  double value;
  if (data.get_tag_value(TAG_FNumber, value)) return value;
  if (data.get_tag_value(TAG_ApertureValue, value))
    return pow(2.0, value * 0.5);
  return -1;
}

double vw::camera::ExifView::get_exposure_time() {
  double value;
  if (data.get_tag_value(TAG_ExposureTime, value)) return value;
  if (data.get_tag_value(TAG_ShutterSpeedValue, value))
    return pow(2.0, -value);
  return -1;
}

int vw::camera::ExifView::get_iso() {
  int value;
  if (data.get_tag_value(TAG_ISOSpeedRatings, value)) return value;
  if (data.get_tag_value(TAG_ExposureIndex, value)) return value;
  // otherwise probably have to dig through MakerNote
  return -1;
}

// APEX equivalents
// TODO: report in some logical way when value doesn't exist
double vw::camera::ExifView::get_aperture_value() {
  double value;
  if (data.get_tag_value(TAG_ApertureValue, value)) return value;
  if (data.get_tag_value(TAG_FNumber, value))
    return 2 * log2(value);
  return -1;
}

double vw::camera::ExifView::get_shutter_speed_value() {
  double value;
  if (data.get_tag_value(TAG_ShutterSpeedValue, value)) return value;
  if (data.get_tag_value(TAG_ExposureTime, value))
    return -log2(value);
  return -1;
}

double vw::camera::ExifView::get_film_sensitivity() {
  double iso = (double)get_iso();
  if (iso < 0) return -1;
  return log2(iso / 3.125);
}

double vw::camera::ExifView::get_brightness_value() {
  double Bv;
  if (data.get_tag_value(TAG_BrightnessValue, Bv)) return Bv;
  double Av = get_aperture_value();
  double Tv = get_shutter_speed_value();
  double Sv = get_film_sensitivity();
  if (Sv == -1) Sv = 5.0; // assume ISO 100
 
 return (Av + Tv - Sv);
}

