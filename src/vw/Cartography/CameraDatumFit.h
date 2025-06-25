// __BEGIN_LICENSE__
//  Copyright (c) 2006-2025, United States Government as represented by the
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

#ifndef __VW_CARTOGRAPHY_CAMERA_DATUM_FIT__
#define __VW_CARTOGRAPHY_CAMERA_DATUM_FIT__

/// \file CameraDatumFit.h

#include <vw/Cartography/Datum.h>

/// Fit a camera to given pixels and ground points. This 
// makes use of the datum, so cannot be in the Camera module.

namespace vw { 

namespace camera {              
  class PinholeModel;
  class OpticalBarModel;
}

namespace cartography {

// Find the best-fitting pinhole model given xyz points and pixel values.
void fitPinhole(std::vector<Vector3> const& xyz_vec,
                double cam_height, double cam_weight, double cam_ctr_weight,
                vw::cartography::Datum const& datum,
                std::vector<double> const& pixel_values,
                vw::camera::PinholeModel & out_cam);

// Find the best-fitting optical bar model given xyz points and pixel values.
void fitOpticalBar(std::vector<Vector3> const& xyz_vec,
                   double cam_height, double cam_weight, double cam_ctr_weight,
                   vw::cartography::Datum const& datum,
                   std::vector<double> const& pixel_values,
                   vw::camera::OpticalBarModel & out_cam);

} // namespace cartography
} // namespace vw

#endif // __VW_CARTOGRAPHY_CAMERA_DATUM_FIT__
