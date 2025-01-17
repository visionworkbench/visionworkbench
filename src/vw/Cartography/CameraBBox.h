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


#ifndef __VW_CARTOGRAPHY_CAMERABBOX_H__
#define __VW_CARTOGRAPHY_CAMERABBOX_H__

/// \file CameraBBox.h Contains logic for bounding box, pixel intersection, and
/// misc utilities.

#include <vw/config.h>
#if defined(VW_HAVE_PKG_CAMERA)

#include <vw/Image/ImageViewRef.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Cartography/GeoReference.h>

#include <boost/shared_ptr.hpp>

namespace vw { namespace cartography {

  // Intersect the ray back-projected from the camera with the datum.
  Vector3 datum_intersection(Datum const& datum,
                             camera::CameraModel const* model,
                             Vector2 const& pix);

  // Return the intersection between the ray emanating from the
  // current camera pixel with the datum ellipsoid. The return value
  // is a map-projected point location (the intermediate between
  // lon-lat-altitude and pixel).
  Vector2 geospatial_intersect(GeoReference const& georef,
                               Vector3 const& camera_ctr, Vector3 const& camera_vec,
                               bool& has_intersection);

  // Find a handful of valid DEM values and average them. It helps later when
  // intersecting with the DEM, especially for Mars, where the DEM heights ca be
  // very far from the datum. 
  double demHeightGuess(vw::ImageViewRef<vw::PixelMask<float>> const& dem);

  // Intersect the ray going from the given camera pixel with a DEM.
  // The return value is a Cartesian point. If the ray goes through a
  // hole in the DEM where there is no data, we return no-intersection
  // or intersection with the datum, depending on whether the variable
  // treat_nodata_as_zero is false or true.
  Vector3 camera_pixel_to_dem_xyz(Vector3 const& camera_ctr, Vector3 const& camera_vec,
                                  vw::ImageViewRef<vw::PixelMask<float>> const& dem_image,
                                  GeoReference const& georef,
                                  bool treat_nodata_as_zero,
                                  bool & has_intersection,
                                  double height_error_tol = 1e-1,  // error in DEM height
                                  double max_abs_tol      = 1e-14, // abs cost fun change b/w iters
                                  double max_rel_tol      = 1e-14,
                                  int num_max_iter        = 100,
                                  Vector3 xyz_guess       = Vector3(),
                                  double height_guess     = 
                                  std::numeric_limits<double>::quiet_NaN());

  /// Compute the bounding box in points (georeference space) that is
  /// defined by georef. Scale is MPP as georeference space is in meters.
  /// - If coords is provided the intersection coordinates will be stored there.
  BBox2 camera_bbox(GeoReference const& georef,
                    boost::shared_ptr<vw::camera::CameraModel> camera_model,
                    int32 cols, int32 rows, float &scale,
                    std::vector<Vector2> *coords=0);

  /// Overload with no scale return
  BBox2 camera_bbox(GeoReference const& dem_georef,
                    boost::shared_ptr<vw::camera::CameraModel> camera_model,
                    int32 cols, int32 rows);

  /// Intersections that take into account DEM topography
  /// - Returns a bounding box in Georeference coordinate system (projected if available)
  ///    containing everything visible in the camera image.
  /// - Computes mean_gsd which is the estimated mean ground resolution of the camera.
  ///   Note that ground resolution in row and col directions can be different for LRO NAC.
  ///   This will just return a mean of the two. 
  ///   The mean_gsd is in GeoReference measurement units (not necessarily meters!)
  /// - If the quick option is enabled, only rays along the image borders will be used
  ///   to perform the computation.
  /// - If coords is provided the intersection coordinates will be stored there.
  // Camera footprint on the ground. See the .h file for details.
  BBox2 camera_bbox(vw::ImageViewRef<vw::PixelMask<float>> const& dem,
                    GeoReference const& dem_georef,
                    GeoReference const& target_georef, // return box in this projection
                    boost::shared_ptr<vw::camera::CameraModel> camera_model,
                    int32 cols, int32 rows, float &mean_gsd,
                    bool quick=false,
                    std::vector<Vector3> *coords = NULL,
                    int num_samples = 1000);
  
  /// Overload of camera_bbox when we don't care about getting the mean_gsd back.
  BBox2 camera_bbox(vw::ImageViewRef<vw::PixelMask<float>> const& dem,
                    GeoReference const& dem_georef,
                    GeoReference const& target_georef,
                    boost::shared_ptr<vw::camera::CameraModel> camera_model,
                    int32 cols, int32 rows);

} // namespace cartography
} // namespace vw

#endif // VW_HAVE_PKG_CAMERA

#endif // __VW_CARTOGRAPHY_CAMERABBOX_H__
