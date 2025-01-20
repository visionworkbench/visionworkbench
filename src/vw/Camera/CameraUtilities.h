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


/// \file CameraUtilities.h
///
/// This file contains miscellaneous functions for working with camera models.
///
#ifndef __VW_CAMERAMODEL_CAMERAUTILITIES_H__
#define __VW_CAMERAMODEL_CAMERAUTILITIES_H__

#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/CAHVModel.h>
#include <vw/Camera/CAHVORModel.h>
#include <vw/Camera/CAHVOREModel.h>
#include <vw/Camera/CameraTransform.h>
#include <vw/Camera/LensDistortion.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Camera/OpticalBarModel.h>
#include <vw/Image/Transform.h>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace vw {
namespace camera {

// These functions are defined in the .cc file:

/// Unpack a vector into a rotation + translation + scale
void vector_to_transform(Vector<double> const & C, 
			 Matrix3x3            & rotation,
			 Vector3              & translation,
			 double               & scale);

// Pack a rotation + translation + scale into a vector
void transform_to_vector(Vector<double>  & C, 
			 Matrix3x3 const & rotation,
			 Vector3   const & translation,
			 double    const & scale);

/// Adjust a given camera so that the xyz points project project to
/// the pixel values.
void fit_camera_to_xyz(std::string const& camera_type,
		       bool refine_camera, 
		       std::vector<Vector3> const& xyz_vec,
		       std::vector<double> const& pixel_values,
		       bool verbose, boost::shared_ptr<CameraModel> & out_cam);

// Given some GCP so that at least two images have at at least three GCP each,
// but each GCP is allowed to show in one image only, use the GCP
// to transform cameras to ground coordinates.
void align_cameras_to_ground(std::vector<std::vector<Vector3>> const& xyz,
                             std::vector<std::vector<Vector2>> const& pix,
                             std::vector<PinholeModel> & sfm_cams,
                             vw::Matrix3x3 & rotation, 
                             vw::Vector3 & translation,
                             double & scale);
  
/// Load a pinhole camera model of any supported type
boost::shared_ptr<CameraModel> load_pinhole_camera_model(std::string const& path);

/// Load a pinhole, CAHV, CAHVOR, or CAHVORE model and convert to CAHV.
boost::shared_ptr<CAHVModel>
load_cahv_pinhole_camera_model(std::string const& image_path,
                               std::string const& camera_path);

/// Compute a good sample spacing for the given input image so that 
/// the undistortion functions below won't take too long.
int auto_compute_sample_spacing(Vector2i const image_size);
 
/// Adjust a pair of epipolar-aligned cameras so that the input images are fully
/// contained in the transformed images.
// TODO(oalexan1): Move to EpipolarTransform.h.
void resize_epipolar_cameras_to_fit(PinholeModel const& cam1,
                                    PinholeModel const& cam2,
                                    PinholeModel      & epi_cam1,
                                    PinholeModel      & epi_cam2,
                                    BBox2i       const& roi1,
                                    BBox2i       const& roi2,
                                    Vector2i          & epi_size1,
                                    Vector2i          & epi_size2);

// Convert an optical model to a pinhole model without distortion
// (The distortion will be taken care of later.)
PinholeModel opticalbar2pinhole(OpticalBarModel const& opb_model, int sample_spacing,
                                double camera_to_ground_dist);

// Apply a given rotation + translation + scale to a pinhole model. We
// assume the model already has set its optical center, focal length,
// etc.
template <class CAM>
void apply_rot_trans_scale(CAM & P, Vector<double> const& C){

  if (C.size() != 7) 
    vw_throw( LogicErr() << "Expecting a vector of size 7.\n" );

  // Parse the transform
  Matrix3x3 rotation;
  Vector3   translation;
  double    scale;
  vector_to_transform(C, rotation, translation, scale);

  // Apply the transform
  P.apply_transform(rotation, translation, scale);
}

// These template functions are defined inline:

// Pack a pinhole or optical bar model's rotation and camera center into a vector
template <class CAM>
void camera_to_vector(CAM const& P, Vector<double> & C){

  Vector3 ctr = P.camera_center();
  Vector3 axis_angle = P.camera_pose().axis_angle();
  C.set_size(6);
  subvector(C, 0, 3) = ctr;
  subvector(C, 3, 3) = axis_angle;
}

// Pack a vector into a pinhole or optical bar model. We assume the model
// already has set its optical center, focal length, etc. 
template <class CAM>
void vector_to_camera(CAM & P, Vector<double> const& C){

  if (C.size() != 6) 
    vw_throw( LogicErr() << "Expecting a vector of size 6.\n" );

  Vector3 ctr      = subvector(C, 0, 3);
  Quat    rotation = axis_angle_to_quaternion(subvector(C, 3, 3));
  P.set_camera_center(ctr);
  P.set_camera_pose(rotation);
}

/// Class to use with the LevenbergMarquardt solver to optimize the parameters of a desired lens
///  to match the passed in pairs of undistorted (input) and distorted (output) points.
template <class DistModelT>
struct DistortionOptimizeFunctor:
    public math::LeastSquaresModelBase< DistortionOptimizeFunctor<DistModelT> > {
  typedef Vector<double> result_type;
  typedef Vector<double> domain_type;
  typedef Matrix<double> jacobian_type;

  const camera::PinholeModel& m_cam;
  const std::vector<Vector2>& m_undist_coords;
  
  /// Init the object with the pinhole model and the list of undistorted image coordinates.
  /// - The list of distorted image coordinates in a Vector<double> (packed alternating col, row) 
  ///    must be passed to the solver function.
  DistortionOptimizeFunctor(const camera::PinholeModel& cam,
                            const std::vector<Vector2>& undist_coords):
    m_cam(cam), m_undist_coords(undist_coords) {}
  
  /// Return a Vector<double> of all the distorted pixel coordinates.
  inline result_type operator()( domain_type const& x ) const {
    DistModelT lens(x); // Construct lens distortion model with given parameters
    result_type out_vec;
    
    out_vec.set_size(m_undist_coords.size()*2);
    for (size_t i=0; i<m_undist_coords.size(); ++i) {
      Vector2 loc = lens.distorted_coordinates(m_cam, m_undist_coords[i]);
      out_vec[2*i  ] = loc[0]; // The col and row values are packed in successive indices.
      out_vec[2*i+1] = loc[1];
    }

    return out_vec;
  }
}; // End class LensOptimizeFunctor

// Fit a pinhole model with one that has a fast distortion function needed for
// quick computation of the point_to_pixel function. Does not replace the camera
// if the approximation error is too high, unless forced. Does not force the
// replacement of distortion for fast models, unless forced.
PinholeModel fitPinholeModel(CameraModel const* in_model, 
                             vw::Vector2 const& image_size,
                             std::string const& out_distortion_type,
                             bool force_conversion,
                             int sample_spacing = 0,
                             int rpc_degree = 0,
                             double camera_to_ground_dist = 0);

// Given two epipolar rectified CAHV camera models and input images,
// get the aligned version of the images suitable for writing to disk and processing.
// TODO(oalexan1): De-templatize. Move to EpipolarTransform.h.
template <class ImageInT, class ImageOutT>
void get_epipolar_transformed_images(std::string const& left_camera_file,
                                     std::string const& right_camera_file,
                                     boost::shared_ptr<CameraModel> const left_cahv_camera,
                                     boost::shared_ptr<CameraModel> const right_cahv_camera,
                                     ImageInT  const& left_image_in,
                                     ImageInT  const& right_image_in,
                                     ImageOutT      & left_image_out,
                                     ImageOutT      & right_image_out,
                                     ValueEdgeExtension<typename ImageOutT::pixel_type> edge_ext) {

  // In the epipolar alignment case, the "camera_models" function returns the CAHVModel type!
  CAHVModel* left_epipolar_cahv  = dynamic_cast<CAHVModel*>(unadjusted_model(&(*left_cahv_camera )));
  CAHVModel* right_epipolar_cahv = dynamic_cast<CAHVModel*>(unadjusted_model(&(*right_cahv_camera)));
  if (!left_epipolar_cahv || !right_epipolar_cahv) {
    vw_throw(ArgumentErr() << "load_cahv_pinhole_camera_model: CAHVModel cast failed!\n");
  }

  // Lens distortion is corrected in the transform from the input camera files
  //  to the input distortionless CAHV camera models.

  std::string lcase_file = boost::to_lower_copy(left_camera_file);

  // Remove lens distortion and create epipolar rectified images.
  if (boost::ends_with(lcase_file, ".cahvore")) {
    CAHVOREModel left_cahvore (left_camera_file );
    CAHVOREModel right_cahvore(right_camera_file);
    left_image_out  = camera_transform(left_image_in,  left_cahvore,  *left_epipolar_cahv,
                                       edge_ext, BilinearInterpolation()
                                       );
    right_image_out = camera_transform(right_image_in, right_cahvore, *right_epipolar_cahv,
                                       edge_ext, BilinearInterpolation()
                                       );
    
  } else if (boost::ends_with(lcase_file, ".cahvor") ||
             boost::ends_with(lcase_file, ".cmod") ) {
    CAHVORModel left_cahvor (left_camera_file );
    CAHVORModel right_cahvor(right_camera_file);
    left_image_out  = camera_transform(left_image_in,  left_cahvor,  *left_epipolar_cahv,
                                       edge_ext, BilinearInterpolation()
                                       );
    right_image_out = camera_transform(right_image_in, right_cahvor, *right_epipolar_cahv,
                                       edge_ext, BilinearInterpolation()
                                       );

  } else if ( boost::ends_with(lcase_file, ".cahv") ||
              boost::ends_with(lcase_file, ".pin" )) {
    CAHVModel left_cahv (left_camera_file );
    CAHVModel right_cahv(right_camera_file);
    left_image_out  = camera_transform(left_image_in,  left_cahv,  *left_epipolar_cahv,
                                       edge_ext, BilinearInterpolation()
                                       );
    right_image_out = camera_transform(right_image_in, right_cahv, *right_epipolar_cahv,
                                       edge_ext, BilinearInterpolation()
                                       );

  } else if ( boost::ends_with(lcase_file, ".pinhole") ||
              boost::ends_with(lcase_file, ".tsai") ) {
    PinholeModel left_pin (left_camera_file );
    PinholeModel right_pin(right_camera_file);

    // If necessary, replace the lens distortion models with a approximated models
    //  that will be much faster in the camera_transform calls below.
    Vector2i left_image_size (left_image_in.cols(),  left_image_in.rows() );
    Vector2i right_image_size(right_image_in.cols(), right_image_in.rows());
    bool force_conversion = false;
    left_pin = fitPinholeModel(&left_pin, left_image_size, "TsaiLensDistortion", 
                               force_conversion);
    right_pin = fitPinholeModel(&right_pin, right_image_size, "TsaiLensDistortion", 
                                force_conversion);
    
    left_image_out  = camera_transform(left_image_in,  left_pin,  *left_epipolar_cahv,
                                       edge_ext, BilinearInterpolation());
    right_image_out = camera_transform(right_image_in, right_pin, *right_epipolar_cahv,
                                       edge_ext, BilinearInterpolation());

  } else {
    vw_throw(ArgumentErr() << "get_epipolar_transformed_images: unsupported camera file type.\n");
  }

}

// TODO(oalexan1): De-templatize. Move to EpipolarTransform.h.
template <class ImageInT, class ImageOutT,  class EdgeT, class InterpT>
void get_epipolar_transformed_pinhole_images(std::string const& left_camera_file,
                                             std::string const& right_camera_file,
                                             boost::shared_ptr<CameraModel> const left_epi_cam,
                                             boost::shared_ptr<CameraModel> const right_epi_cam,
                                             ImageInT  const& left_image_in,
                                             ImageInT  const& right_image_in,
                                             BBox2i    const& left_image_in_roi,
                                             BBox2i    const& right_image_in_roi,
                                             Vector2i       & left_image_out_size,
                                             Vector2i       & right_image_out_size,
                                             ImageOutT      & left_image_out,
                                             ImageOutT      & right_image_out,
                                             EdgeT const& edge_func, InterpT const& interp_func) {

  // Load the original camera models with distortion
  boost::shared_ptr<PinholeModel> left_pin (new PinholeModel(left_camera_file ));
  boost::shared_ptr<PinholeModel> right_pin(new PinholeModel(right_camera_file));

  // If necessary, replace the lens distortion models with a approximated models
  //  that will be much faster in the camera_transform calls below.
  Vector2i left_image_in_size (left_image_in.cols(),  left_image_in.rows() );
  Vector2i right_image_in_size(right_image_in.cols(), right_image_in.rows());
  bool force_conversion = false;
  *left_pin.get() = fitPinholeModel(left_pin.get(), left_image_in_size, 
                                    "TsaiLensDistortion", force_conversion);
  *right_pin.get() = fitPinholeModel(right_pin.get(), right_image_in_size, 
                                     "TsaiLensDistortion", force_conversion);
 
  // Make sure there are no adjustments on the aligned camera models
  PinholeModel* left_epi_pin  = dynamic_cast<PinholeModel*>(unadjusted_model(&(*left_epi_cam)));
  PinholeModel* right_epi_pin = dynamic_cast<PinholeModel*>(unadjusted_model(&(*right_epi_cam)));
  
  // If so, create an adjusted version of the input files from disk and then transform, otherwise just transform.
  Vector3 ZERO_TRANSLATION(0,0,0);
  Quat    ZERO_ROTATION(1,0,0,0);
  double  NO_SCALE = 1.0;
  if (left_image_in_roi.min() != Vector2i(0,0)) {
    // Get the input raw camera model with the crop applied.
    AdjustedCameraModel left_adj_cam(left_pin, ZERO_TRANSLATION,  ZERO_ROTATION, left_image_in_roi.min(), NO_SCALE);
    
    // Convert from the cropped input to the epipolar-aligned camera model.
    left_image_out  = camera_transform(left_image_in,  left_adj_cam,  *left_epi_pin,
                                     left_image_out_size, edge_func, interp_func);
  }
  else { // Don't need to handle crops.
    left_image_out  = camera_transform(left_image_in,  *left_pin,  *left_epi_pin,
                                       left_image_out_size, edge_func, interp_func);
  }
  // Repeat for the right camera.
  if (right_image_in_roi.min() != Vector2i(0,0)) {
    AdjustedCameraModel right_adj_cam(right_pin, ZERO_TRANSLATION,  ZERO_ROTATION, right_image_in_roi.min(), NO_SCALE);

    right_image_out = camera_transform(right_image_in, right_adj_cam, *right_epi_pin,
                                       right_image_out_size, edge_func, interp_func);
  }
  else { // No crops
        right_image_out = camera_transform(right_image_in, *right_pin, *right_epi_pin,
                                       right_image_out_size, edge_func, interp_func);
  }
}

}}      // namespace vw::camera

#endif  //__CAMERAMODEL_CAMERAUTILITIES_H__
