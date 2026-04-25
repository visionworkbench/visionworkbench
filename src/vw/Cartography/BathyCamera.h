// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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

// Camera-aware bathymetry helpers. Each function takes a camera model and
// a vw::BathyPlane and answers a single ray/pixel question with Snell's
// law refraction applied where appropriate. Pure plane and refraction
// math lives in SnellLaw and BathyData; this file is the layer that joins
// those primitives to a vw::camera::CameraModel. The BathyStereoModel
// triangulator is a separate consumer and lives in its own file.

#ifndef __VW_CARTOGRAPHY_BATHYCAMERA_H__
#define __VW_CARTOGRAPHY_BATHYCAMERA_H__

#include <vw/Math/Vector.h>
#include <vw/Cartography/BathyData.h>

namespace vw {

namespace camera {
  class CameraModel;
}

// Intersect a ray from camera center along camera direction with the datum at
// given semi-axes, with optional bathymetry correction. If the ray passes
// through the bathy plane (water surface) before it meets the datum, apply
// Snell's law refraction and continue with the new bent ray until reaching the
// datum. Returns the intersection point, or zero vector on failure.
vw::Vector3 datumBathyIntersection(vw::Vector3 const& cam_ctr,
                                   vw::Vector3 const& cam_dir,
                                   double major_axis, double minor_axis,
                                   BathyPlane const& bathy_plane,
                                   double refraction_index);

// Project an ECEF point to pixel, accounting for bathymetry if the point
// is below the bathy plane (water surface).
vw::Vector2 point_to_pixel(vw::camera::CameraModel const* cam,
                           vw::BathyPlane const& bathy_plane,
                           double refraction_index,
                           vw::Vector3 const& ecef_point);

} // namespace vw

#endif // __VW_CARTOGRAPHY_BATHYCAMERA_H__
