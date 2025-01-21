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

#include <vw/Cartography/CameraDatumFit.h>
#include <vw/Math/LevenbergMarquardt.h>

#include <vw/Camera/CameraParamsPack.h>

namespace vw { namespace cartography {

using namespace vw::camera;
using namespace vw;

// Solve for best fitting camera that projects given xyz locations at
// given pixels. If cam_weight > 0, try to constrain the camera height
// above datum at the value of cam_height.
// If camera center weight is given, use that to constrain the
// camera center, not just its height.
template <class CAM>
class CameraSolveLmaHt: public vw::math::LeastSquaresModelBase<CameraSolveLmaHt<CAM>> {
  std::vector<vw::Vector3> const& m_xyz;
  CAM m_camera_model;
  double m_cam_height, m_cam_weight, m_cam_ctr_weight;
  vw::cartography::Datum m_datum;
  Vector3 m_input_cam_ctr;

public:

  typedef vw::Vector<double>    result_type;   // pixel residuals
  typedef vw::Vector<double, 6> domain_type;   // camera parameters (camera center and axis angle)
  typedef vw::Matrix<double> jacobian_type;

  /// Instantiate the solver with a set of xyz to pixel pairs and a pinhole model
  CameraSolveLmaHt(std::vector<vw::Vector3> const& xyz,
                    CAM const& camera_model,
                    double cam_height, double cam_weight, double cam_ctr_weight,
                    vw::cartography::Datum const& datum):
    m_xyz(xyz),
    m_camera_model(camera_model),
    m_cam_height(cam_height), m_cam_weight(cam_weight), m_cam_ctr_weight(cam_ctr_weight),
    m_datum(datum), m_input_cam_ctr(m_camera_model.camera_center(vw::Vector2())) {}

  /// Given the camera, project xyz into it
  inline result_type operator()(domain_type const& C) const {

    // Create the camera model
    CAM camera_model = m_camera_model;  // make a copy local to this function
    vector_to_camera(camera_model, C);  // update its parameters

    int xyz_len = m_xyz.size();
    size_t result_size = xyz_len * 2;
    if (m_cam_weight > 0)
      result_size += 1;
    else if (m_cam_ctr_weight > 0)
      result_size += 3; // penalize deviation from original camera center

    // See where the xyz coordinates project into the camera.
    result_type result;
    result.set_size(result_size);
    for (size_t i = 0; i < xyz_len; i++) {
      Vector2 pixel = camera_model.point_to_pixel(m_xyz[i]);
      result[2*i  ] = pixel[0];
      result[2*i+1] = pixel[1];
    }

    if (m_cam_weight > 0) {
      // Try to make the camera stay at given height
      Vector3 cam_ctr = subvector(C, 0, 3);
      Vector3 llh = m_datum.cartesian_to_geodetic(cam_ctr);
      result[2*xyz_len] = m_cam_weight*(llh[2] - m_cam_height);
    } else if (m_cam_ctr_weight > 0) {
      // Try to make the camera stay close to given center
      Vector3 cam_ctr = subvector(C, 0, 3);
      for (int it = 0; it < 3; it++)
        result[2*xyz_len + it] = m_cam_ctr_weight*(m_input_cam_ctr[it] - cam_ctr[it]);
    }

    return result;
  }
}; // End class CameraSolveLmaHt

// Find the best-fitting camera given xyz points and pixel values.
template<class CAM>
void fitCam(std::vector<Vector3> const& xyz_vec,
            double cam_height, double cam_weight, double cam_ctr_weight,
            vw::cartography::Datum const& datum,
            std::vector<double> const& pixel_values,
            CAM & out_cam) {

  Vector<double> out_vec; // must copy to this structure
  int residual_len = pixel_values.size();

  if (cam_weight > 0.0)
    residual_len += 1; // for camera height residual
  else if (cam_ctr_weight > 0)
    residual_len += 3; // for camera center residual

  // Copy the image pixels
  out_vec.set_size(residual_len);
  for (size_t corner_it = 0; corner_it < pixel_values.size(); corner_it++)
    out_vec[corner_it] = pixel_values[corner_it];

  // Use 0 for the remaining fields corresponding to camera height or
  // camera center constraint
  for (int it = pixel_values.size(); it < residual_len; it++)
    out_vec[it] = 0.0;

  double abs_tolerance  = 1e-24;
  double rel_tolerance  = 1e-24;
  int    max_iterations = 2000;
  int status = 0;
  Vector<double> final_params;
  Vector<double> seed;

  CameraSolveLmaHt<CAM>
    lma_model(xyz_vec, out_cam, cam_height, cam_weight, cam_ctr_weight, datum);
  camera_to_vector(out_cam, seed);

  final_params = math::levenberg_marquardt(lma_model, seed, out_vec,
              status, abs_tolerance, rel_tolerance, max_iterations);

  vector_to_camera(out_cam, final_params);

  if (status < 1)
    vw_out() << "The Levenberg-Marquardt solver failed. Results may be inaccurate.\n";

  return;
}

// Find the best-fitting optical bar model given xyz points and pixel values.
// This is wrapper around the above function. It is put here to avoid frequent
// recompilation of templates.
void fitOpticalBar(std::vector<Vector3> const& xyz_vec,
                   double cam_height, double cam_weight, double cam_ctr_weight,
                   vw::cartography::Datum const& datum,
                   std::vector<double> const& pixel_values,
                   vw::camera::OpticalBarModel & out_cam) {
  fitCam<vw::camera::OpticalBarModel>(xyz_vec, cam_height, cam_weight, cam_ctr_weight, 
                                      datum, pixel_values, out_cam);
}

// Find the best-fitting pinhole model given xyz points and pixel values.
// Also a wrapper, as above.
void fitPinhole(std::vector<Vector3> const& xyz_vec,
                double cam_height, double cam_weight, double cam_ctr_weight,
                vw::cartography::Datum const& datum,
                std::vector<double> const& pixel_values,
                vw::camera::PinholeModel & out_cam) {

  fitCam<vw::camera::PinholeModel>(xyz_vec, cam_height, cam_weight, cam_ctr_weight, 
                                   datum, pixel_values, out_cam);
}
}} // end namespace vw::cartography
