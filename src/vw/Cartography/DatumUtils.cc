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

#include <vw/Cartography/DatumUtils.h>
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>

namespace vw {

// This can catch user mistakes. If warn_only is false, and the error is huge,
// this will throw an exception. Otherwise it will print a warning.
void checkDatumConsistency(vw::cartography::Datum const& datum1,                        
                           vw::cartography::Datum const& datum2,
                           bool warn_only) {

  double err1 = std::abs(datum1.semi_major_axis() - datum2.semi_major_axis());
  double err2 = std::abs(datum1.semi_minor_axis() - datum2.semi_minor_axis());
  double err = std::max(err1, err2);
  
  // Small bodies can be tricky, for those allow various radii
  double rad = 1.0e+6; // 1000 km
  bool small_body = (datum1.semi_major_axis() < rad || datum1.semi_minor_axis() < rad ||
                     datum2.semi_major_axis() < rad || datum2.semi_minor_axis() < rad);
                 
  // Do not warn for small errors, like between WGS84 and 3D CRS. This is meant
  // to catch gross problems, such as planet mixup.                   
  if (err >= 10.0) {
    std::ostringstream oss;
    oss.precision(8);
    oss << "Found two distinct datums. The difference in semi-axes is: " 
        << err << " meters.\n"
        << "Datum 1: " << datum1 << "\n"
        << "Datum 2: " << datum2 << "\n";
    // Because in difference in semi-axes, datums can easily differ by 50 km or so,
    // so be be generous with the threshold.    
    if (err < 200000.0 || warn_only || small_body) // this is mild
       vw::vw_out(vw::WarningMessage) 
       << oss.str() 
       << "This is likely harmless, but check your inputs.\n";
    else // this is severe
      vw::vw_throw(vw::ArgumentErr() << oss.str());
  }
}

} // end namespace vw
