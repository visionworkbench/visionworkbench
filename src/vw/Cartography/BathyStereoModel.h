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

// This file is here because it depends on the Cartography module.

#ifndef __VW_CARTOGRAPHY_BATHYSTEREOMODEL_H__
#define __VW_CARTOGRAPHY_BATHYSTEREOMODEL_H__

#include <vw/Math/Vector.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelMask.h>

namespace vw {

// Forward declaration
namespace camera {
  class CameraModel;
}

// Check if the given left and right pixels are in the masked region (invalid in
// the mask). That will mean bathymetry correction should be applied.
bool areMasked(ImageViewRef<PixelMask<float>> const& left_mask,
               ImageViewRef<PixelMask<float>> const& right_mask,
               Vector2 const& lpix, Vector2 const& rpix);

struct BathyPlane {
  std::vector<double> bathy_plane;
  bool use_curved_water_surface;
  vw::cartography::GeoReference water_surface_projection;

  BathyPlane(): use_curved_water_surface(false) {}
};

// A struct to hold the bathymetry settings and data
struct BathyData {
  std::vector<vw::ImageViewRef<vw::PixelMask<float>>> bathy_masks;
  std::vector<BathyPlane> bathy_planes;
  float refraction_index;
  BathyData(): refraction_index(1.0) {}
};

// Read a bathy mask. Return the nodata value.
vw::ImageViewRef<vw::PixelMask<float>> read_bathy_mask(std::string const& filename,
                                                       float & nodata_val);

// Read a set of bathy masks.
void read_bathy_masks(std::vector<std::string> const& mask_filenames,
                      std::vector<vw::ImageViewRef<vw::PixelMask<float>>> & bathy_masks);

// Read the bathy planes and associated data. More often than not they will be
// identical. If there is more than one bathy plane file, they are all kept in
// the same string, separated by space.
void readBathyPlanes(std::string const& bathy_plane_files,
                     int num_images,
                     std::vector<BathyPlane> & bathy_plane_vec);

// Given a ray going down towards Earth, starting at point in_xyz and
// with unit direction in_dir, a plane 'p' to the water surface with four
// coefficients such that the plane equation is p[0] * x + p[1] * y
// + p[2] * z + p[3] = 0, the normal (p[0], p[1], p[2]) pointing
// upwards away from Earth, and water refraction index, find where
// this ray meets the water plane named out_xyz, and the ray direction out_dir
// after it bends according to Snell's law. Return true on success.
bool snellLaw(vw::Vector3 const& in_xyz, vw::Vector3 const& in_dir,
              std::vector<double> const& plane, double refraction_index,
              vw::Vector3 & out_xyz, vw::Vector3 & out_dir);

class BathyStereoModel: public vw::stereo::StereoModel {
public:

  //------------------------------------------------------------------
  // Constructors / Destructors
  //------------------------------------------------------------------
  BathyStereoModel(std::vector<const vw::camera::CameraModel *> const& cameras,
                    double angle_tol = 0.0):
    vw::stereo::StereoModel(cameras, angle_tol),
    m_bathy_correct(false) {}

  BathyStereoModel(vw::camera::CameraModel const* camera_model1,
                    vw::camera::CameraModel const* camera_model2,
                    double angle_tol = 0.0):
    vw::stereo::StereoModel(camera_model1, camera_model2, angle_tol),
    m_bathy_correct(false), m_single_bathy_plane(true) {}

  virtual ~BathyStereoModel() {}

  //------------------------------------------------------------------
  // Public Methods
  //------------------------------------------------------------------

  /// Apply a stereo model to multiple or just two image coordinates.
  /// Returns an xyz point. The error is set to 0 if triangulation
  /// did not succeed, otherwise it is the vector between the closest points on the rays.
  virtual vw::Vector3 operator()(std::vector<vw::Vector2> const& pixVec,
                                  vw::Vector3& errorVec, bool do_bathy, bool & did_bathy) const;
  virtual vw::Vector3 operator()(std::vector<vw::Vector2> const& pixVec,
                                  double & error) const;
  virtual vw::Vector3 operator()(vw::Vector2 const& pix1, vw::Vector2 const& pix2,
                                  vw::Vector3& errorVec) const;
  virtual vw::Vector3 operator()(vw::Vector2 const& pix1, vw::Vector2 const& pix2,
                                  double & error) const;

  // Settings used for bathymetry correction. The left and right images
  // get individual bathy plane settings, but they may be identical.
  void set_bathy(double refraction_index,
                 std::vector<BathyPlane> const& bathy_plane_vec);

private:
  bool m_bathy_correct;                        // If to do bathy correction
  bool m_single_bathy_plane;                   // if the left and right images use same plane 
  double m_refraction_index;                   // Water refraction index
  std::vector<BathyPlane> m_bathy_plane_vec;   // Bathy plane settings
};

} // namespace vw

#endif  // __VW_CARTOGRAPHY_BATHYSTEREOMODEL_H__
