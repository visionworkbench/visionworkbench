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

#include <vw/Cartography/Utm.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>

#include <limits>
#include <cmath>

namespace vw {
namespace cartography {

// Auxiliary function, it does not give the final zone
int getLongitudeZone(double longitude) {
  
  // Check for nan and inf
  if (std::isnan(longitude) || std::isinf(longitude)) 
    vw::vw_throw(vw::ArgumentErr() << "Invalid longitude: " << longitude << ".\n");
    
  while (longitude < -180) longitude += 360;
  while (longitude > 180) longitude -= 360;
  return static_cast<int>((longitude + 180) / 6) + 1;
}
  
char getLatitudeBand(double latitude) {

  if (latitude < -90) latitude = -90;
  if (latitude > 90) latitude = 90;
  
  const char bands[] = "CDEFGHJKLMNPQRSTUVWX";
  
  if (latitude < -80) return 'A';
  if (latitude > 84) return 'Z';
  
  int bandIndex = static_cast<int>((latitude + 80) / 8);
  return bands[bandIndex];
}

// Find the UTM zone. Must be called for the WGS84 datum,
// with lat <= 84 && lat >= -80.
int getUTMZone(double latitude, double longitude) {
  
  int zone = getLongitudeZone(longitude);
  char band = getLatitudeBand(latitude);
  
  // Special case for Norway
  if (band == 'V' && zone == 31 && longitude >= 3)
      zone = 32;

  // Special case for Svalbard
  if (band == 'X') {
    if (longitude >= 0 && longitude < 9) zone = 31;
    else if (longitude >= 9 && longitude < 21) zone = 33;
    else if (longitude >= 21 && longitude < 33) zone = 35;
    else if (longitude >= 33 && longitude < 42) zone = 37;
  }
  
  // Return the zone only, without the band
  return zone;
}

// Compare with https://www.geoplaner.com/
void testUtm() {

  double lat = 0.0, lon = 0.0;
  
  // California
  lat = 35.2; lon = -118.2;
  vw::vw_out() << "California: expect: " << 11 
    << " got: " << getUTMZone(lat, lon) << "\n";
  
  // Grand Mesa1
  lat = 39.026667; lon = -108.081389;
  vw::vw_out() << "Grand Mesa 1: expect: " << 12 
    << " got: " << getUTMZone(lat, lon) << "\n";

  // Grand Mesa2
  lat = 39.026667; lon = -108.0;
  vw::vw_out() << "Grand Mesa 2: expect: " << 13 
    << " got: " << getUTMZone(lat, lon) << "\n";

  // Norway exception
  lat = 60.5; lon = 4.0;
  vw::vw_out() << "Norway special case: expect: " << 32 
    << " got: " << getUTMZone(lat, lon) << "\n";  

  // Svalbard exception 1
  lat = 75.0; lon = 8.0;
  vw::vw_out() << "Svalbard special case 1, expect: " << 31 
  << " got: " << getUTMZone(lat, lon) << "\n";
  
  // Svalbard exception 2
  lat = 75.0; lon = 15.0;
  vw::vw_out() << "Svalbard special case 2, expect: " << 33
    << " got: " << getUTMZone(lat, lon) << "\n";

  // Svalbard exception 3
  lat = 75.0; lon = 25.0;
  vw::vw_out() << "Svalbard special case 3, expect: " << 35
    << " got: " << getUTMZone(lat, lon) << "\n"; 
 
 // Svalbard exception 4
 lat = 75.0; lon = 35.0;
  vw::vw_out() << "Svalbard special case 4, expect: " << 37
    << " got: " << getUTMZone(lat, lon) << "\n";
} 

}} // end namespace vw::cartography
