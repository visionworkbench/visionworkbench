// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifndef __VW_FLOOD_COMMON_H__
#define __VW_FLOOD_COMMON_H__

#include <stdlib.h>
#include <algorithm>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Transform.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReferenceUtils.h>

/**
  Common flood detection tools.
*/

namespace vw {

//=============================================================================
// Constants

const uint8 FLOOD_DETECT_WATER  = 255;
const uint8 FLOOD_DETECT_LAND   = 1;
const uint8 FLOOD_DETECT_NODATA = 0;

const double RAD_TO_DEG = 180.0 / M_PI; // TODO: Where does GeoReference.cc get this from?
const double DEG_TO_RAD = M_PI / 180.0;

//=============================================================================
// Functions

/// Return the first item in a list of strings that contains the given substring.
std::string find_string_in_list(std::vector<std::string> const& input_strings,
                                std::string const& search_string) {
  for (size_t f=0; f<input_strings.size(); ++f) {
    if (input_strings[f].find(search_string) != std::string::npos) {
      return input_strings[f];
    }
  }
  return "";
}

/// Extract numeric the value from one line of the metadata file.
float parse_metadata_line(std::string const& line) {
  // Parse the line containing the info
  size_t eqpos = line.find("=");
  std::string num = line.substr(eqpos+1);
  return atof(num.c_str());
}

std::string num2str(int n) {
  std::stringstream s;
  s << n;
  return s.str();
}

/// Helper function to compute one of the many band derived indices.
float compute_index( float a, float b) {
  float denom = a + b;
  if (denom == 0)
    return 100; // Avoid divide-by-zero
  return (a - b) / denom;
}


// TODO: MOVE THIS
template <typename T>
T clamp(T value, T min, T max) {
  if (value > max) return max;
  if (value < min) return min;
  return value;
}

/// Shortcut to clamp between 0 and 1
template <typename T>
T clamp01(T value) {
  return clamp<T>(value, 0, 1);
}

/// Scale the given input range to the 0 to 1 range.
float rescale_to_01(float value, float min, float max) {
  float rng = max - min;
  return ((value - min) / rng);
}

/// Computes the distance between the Earth and the Sun at a given time in AU.
/// - Copied from "Radiometric Use of WorldView-2 Imagery"
double compute_earth_sun_distance(int year, int month, int day, int hour, int minute, double second) {
  double ut = (double)hour + (double)minute/60.0 + second/3600.0;
  if ((month == 1) or (month == 2)) {
    --year;
    month += 12;
  }
  int a = year/100.0;
  int b = 2 -a + (a/4);
  double julian_day = floor(365.25*(year+4716)) + floor(30.6001*(month+1)) 
                      + day + ut/24.0 + b - 1524.5;                       
  double d = julian_day - 2451545.0;
  double g = 357.529 + 0.98560028*d;
  double earth_sun_distance = 1.00014 - 0.01671*cos(g) - 0.00014*cos(2*g);
  //std::cout << "julian_day = " << julian_day << std::endl;
  //std::cout << "earth_sun_distance = " << earth_sun_distance << std::endl;
  return earth_sun_distance;
}


} // end namespace vw

#endif

