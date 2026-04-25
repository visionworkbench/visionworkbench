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

// The BathyStereoModel class. Camera-aware bathy ray helpers
// (datumBathyIntersection, point_to_pixel) are declared in
// vw/Cartography/BathyCamera.h, which this header transitively
// pulls in for backward compatibility with consumers that used
// to get them from BathyStereoModel.h.

#ifndef __VW_CARTOGRAPHY_BATHYSTEREOMODEL_H__
#define __VW_CARTOGRAPHY_BATHYSTEREOMODEL_H__

#include <vw/Math/Vector.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Cartography/BathyData.h>
#include <vw/Cartography/BathyCamera.h>

namespace vw {

namespace camera {
  class CameraModel;
}

typedef boost::shared_ptr<camera::CameraModel> CamPtr;

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
