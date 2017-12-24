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
#include <vw/Image/Transform.h>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace vw {
namespace camera {

/// Load a pinhole camera model of any supported type
boost::shared_ptr<vw::camera::CameraModel> load_pinhole_camera_model(std::string const& path);

/// Load a pinhole, CAHV, CAHVOR, or CAHVORE model and convert to CAHV.
boost::shared_ptr<vw::camera::CAHVModel>
load_cahv_pinhole_camera_model(std::string const& image_path,
                               std::string const& camera_path);
 
/// Class to use with the LevenbergMarquardt solver to optimize the parameters of a desired lens
///  to match the passed in pairs of undistorted (input) and distorted (output) points.
template <class DistModelT, int NumModelParamsT>
struct DistortionOptimizeFunctor :  public math::LeastSquaresModelBase< DistortionOptimizeFunctor<DistModelT, NumModelParamsT> > {
  typedef Vector<double> result_type;
  typedef Vector<double, NumModelParamsT> domain_type;
  typedef Matrix<double> jacobian_type;

  const camera::PinholeModel& m_cam;
  const std::vector<Vector2>& m_undist_coords;
  
  /// Init the object with the pinhole model and the list of undistorted image coordinates.
  /// - The list of distorted image coordinates in a Vector<double> (packed alternating col, row) 
  ///    must be passed to the solver function.
  DistortionOptimizeFunctor(const camera::PinholeModel& cam, const std::vector<Vector2>& undist_coords)
    : m_cam(cam), m_undist_coords(undist_coords) {}
  
  /// Return a Vector<double> of all the distorted pixel coordinates.
  inline result_type operator()( domain_type const& x ) const {
    DistModelT lens(x); // Construct lens distortion model with given parameters
    result_type        out_vec;
    
    out_vec.set_size(m_undist_coords.size()*2);
    for (size_t i=0; i<m_undist_coords.size(); ++i) {
      Vector2 loc = lens.distorted_coordinates(m_cam, m_undist_coords[i]);
      out_vec[2*i  ] = loc[0]; // The col and row values are packed in successive indices.
      out_vec[2*i+1] = loc[1];
    }
    return out_vec;
  }
}; // End class LensOptimizeFunctor

/// Class to solve for undistortion coefficients. Applicable only for RPCLensDistortion.
template <class DistModelT, int NumModelParamsT>
struct UndistortionOptimizeFunctor :  public math::LeastSquaresModelBase< UndistortionOptimizeFunctor<DistModelT, NumModelParamsT> > {
  typedef Vector<double> result_type;
  typedef Vector<double, NumModelParamsT> domain_type;
  typedef Matrix<double> jacobian_type;

  const camera::PinholeModel& m_cam;
  const std::vector<Vector2>& m_dist_coords;
  
  /// Init the object with the pinhole model and the list of undistorted image coordinates.
  /// - The list of distorted image coordinates in a Vector<double> (packed alternating col, row) 
  ///    must be passed to the solver function.
  UndistortionOptimizeFunctor(const camera::PinholeModel& cam, const std::vector<Vector2>& dist_coords)
    : m_cam(cam), m_dist_coords(dist_coords) {}
  
  /// Return a Vector<double> of all the distorted pixel coordinates.
  inline result_type operator()( domain_type const& x ) const {
    DistModelT lens;
    lens.set_undistortion_parameters(x); // Construct lens distortion model with given parameters
    result_type        out_vec;
    
    out_vec.set_size(m_dist_coords.size()*2);
    for (size_t i=0; i< m_dist_coords.size(); ++i) {
      Vector2 loc = lens.undistorted_coordinates(m_cam, m_dist_coords[i]);
      out_vec[2*i  ] = loc[0]; // The col and row values are packed in successive indices.
      out_vec[2*i+1] = loc[1];
    }
    return out_vec;
  }
}; // End class LensOptimizeFunctor

// A very analogous function to
// update_pinhole_for_fast_point2pixel. This one computes the
// undistortion coefficients.  Only applicable for RPC.
template<class DistModelT, int NumModelParamsT>
double compute_undistortion(PinholeModel& pin_model, Vector2i image_size,
                            int sample_spacing=50) {

  // Get info on existing distortion model
  const vw::camera::LensDistortion* input_distortion = pin_model.lens_distortion();
  std::string lens_name = input_distortion->name();
  
  // Check for all of the models that currently support a fast distortion function.
  // - The other models use a solver for this function, greatly increasing the run time.
  if (lens_name != RPCLensDistortion::class_name() && lens_name != RPCLensDistortion5::class_name() && lens_name != RPCLensDistortion6::class_name())  
    vw_throw(ArgumentErr() << "Undistortion can only be computed for RPC models.\n");
  
  // Get input camera information
  const double pixel_pitch = pin_model.pixel_pitch();

  // Generate a set of point pairs
  std::vector<Vector2> distorted_coords;
  Vector<double> undistorted_coords;
  int num_cols   = image_size[0]/sample_spacing; // Grab point pairs for the solver at every
  int num_rows   = image_size[1]/sample_spacing; // interval of "sample_spacing"
  int num_coords = num_rows*num_cols;
  distorted_coords.resize(num_coords);
  undistorted_coords.set_size(num_coords*2);

  // TODO: Need to get boxes of where the RPC model is valid!!!
  
  int index = 0;
  for (int r=0; r<num_rows; ++r) {
    int row = r*sample_spacing;
    for (int c=0; c<num_cols; ++c) {
      int col = c*sample_spacing;

      // TODO: Need more care below to fill up the domain of validity!!!
      // Generate an distorted/undistorted point pair using the input model.
      // - Note that the pixel pairs need to be corrected for pitch here.
      Vector2 undist_pixel = Vector2(col, row) * pixel_pitch;
      Vector2 dist_pixel   = input_distortion->distorted_coordinates(pin_model, undist_pixel);
      //std::cout << "Undist: " << undist_pixel << " ---- Distorted: " << dist_pixel << std::endl;

      distorted_coords[index] = dist_pixel;
      // Pack these points into 1D vector alternating col,row      
      undistorted_coords[2*index  ] = undist_pixel[0]; 
      undistorted_coords[2*index+1] = undist_pixel[1];
      ++index;
    }
  } // End loop through sampled pixels

  // Now solve for a complementary lens distortion scheme

  // Init solver object with the distorted coordinates
  UndistortionOptimizeFunctor<DistModelT, NumModelParamsT> solver_model(pin_model,
                                                                       distorted_coords);
  int status;
  Vector<double, NumModelParamsT> seed; // Start with all zeros (no distortion)

  // Solve for the best new model params that give us the undistorted coordinates from the distorted coordinates.
  Vector<double, NumModelParamsT> undist_params
    = math::levenberg_marquardt(solver_model, seed, undistorted_coords, status);

  // Make a copy of the model
  boost::shared_ptr<LensDistortion> dist_copy = input_distortion->copy();
  DistModelT new_model = *((DistModelT*)(dist_copy.get()));

  // Set the new undistortion params
  new_model.set_undistortion_parameters(undist_params);
  new_model.set_image_size(image_size);
  
  // Check the error
  double diff = 0;
  for (size_t i=0; i < distorted_coords.size(); ++i) {
    Vector2 undistorted        = new_model.undistorted_coordinates(pin_model, distorted_coords[i]);
    Vector2 actual_undistorted = Vector2(undistorted_coords[2*i], undistorted_coords[2*i+1]);
    //std::cout << "New: " << undistorted << " ---- Actual: " << actual_undistorted
    //          << " error: " << norm_2(undistorted - actual_undistorted) << std::endl;
    diff += norm_2(undistorted - actual_undistorted);
  }
  diff /= static_cast<double>(distorted_coords.size());
  vw_out() << DistModelT::class_name() << " undistortion approximation mean error: "
           << diff << ".\n";
  
  pin_model.set_lens_distortion(new_model);
  return diff;
} 
  
// For RPC, must always ensure the undistortion coefficients are up to date.
// This function will be used before saving or displaying a pinhole model.
// This function is not in the .cc file as it is related to the above.
inline void update_rpc_undistortion(PinholeModel const& model){

  const vw::camera::LensDistortion* distortion = model.lens_distortion();
  std::string lens_name = distortion->name();
  if (lens_name != RPCLensDistortion::class_name() && lens_name != RPCLensDistortion5::class_name() && lens_name != RPCLensDistortion6::class_name())
    return;
  
  // Have to cast away the const-ness. Not nice. Only happens for RPC
  // distortion.
  PinholeModel * pin_ptr = const_cast<PinholeModel*>(&model);

  if (lens_name == RPCLensDistortion::class_name()) {
    RPCLensDistortion * rpc_dist = dynamic_cast<RPCLensDistortion*>
      (const_cast<LensDistortion*>(distortion));
    if (rpc_dist == NULL) 
      vw_throw( ArgumentErr() << "PinholeModel::expecting an " + RPCLensDistortion::class_name() +
                " model." );
    
    // Only update this if we have to
    if (!rpc_dist->can_undistort()) 
      compute_undistortion<RPCLensDistortion, RPCLensDistortion::num_distortion_params>
        (*pin_ptr, rpc_dist->image_size());
  }

  if (lens_name == RPCLensDistortion5::class_name()) {
    RPCLensDistortion5 * rpc_dist = dynamic_cast<RPCLensDistortion5*>
      (const_cast<LensDistortion*>(distortion));
    if (rpc_dist == NULL) 
      vw_throw( ArgumentErr() << "PinholeModel::expecting an " + RPCLensDistortion5::class_name() +
                " model." );
    
    // Only update this if we have to
    if (!rpc_dist->can_undistort()) 
      compute_undistortion<RPCLensDistortion5, RPCLensDistortion5::num_distortion_params>
        (*pin_ptr, rpc_dist->image_size());
  }
  
  if (lens_name == RPCLensDistortion6::class_name()) {
    RPCLensDistortion6 * rpc_dist = dynamic_cast<RPCLensDistortion6*>
      (const_cast<LensDistortion*>(distortion));
    if (rpc_dist == NULL) 
      vw_throw( ArgumentErr() << "PinholeModel::expecting an " + RPCLensDistortion6::class_name() +
                " model." );
    
    // Only update this if we have to
    if (!rpc_dist->can_undistort()) 
      compute_undistortion<RPCLensDistortion6, RPCLensDistortion6::num_distortion_params>
        (*pin_ptr, rpc_dist->image_size());
  }

}

/// If necessary, replace the lens distortion model in the pinhole camera model
///  with an approximated model that has a fast distortion function needed for
///  quick computation of the point_to_pixel function.
/// - Does not replace the camera if the approximation error is too high.
/// - Returns the approximation error.
template<class DistModelT, int NumModelParamsT>
double update_pinhole_for_fast_point2pixel(PinholeModel& pin_model, Vector2i image_size,
                                           int sample_spacing=50) {

  // Get info on existing distortion model
  const vw::camera::LensDistortion* input_distortion = pin_model.lens_distortion();
  std::string lens_name = input_distortion->name();
  
  // Check for all of the models that currently support a fast distortion function.
  // - The other models use a solver for this function, greatly increasing the run time.
  if ( (lens_name == "NULL")                           ||
       (lens_name == "TSAI")                           ||
       (lens_name == "AdjustableTSAI")                 ||
       (lens_name == RPCLensDistortion::class_name())  ||
       (lens_name == RPCLensDistortion5::class_name()) ||
       (lens_name == RPCLensDistortion6::class_name())
       ) {
    //vw_out() << "Input distortion is: " << lens_name << ". Refusing to run.\n";
    return 0;
  }
  
  // Otherwise we will try to approximate the input lens with a lens model of type DistModelT.
  
  // Get input camera information
  const double pixel_pitch = pin_model.pixel_pitch();

  // Generate a set of point pairs
  std::vector<Vector2> undistorted_coords;
  Vector<double> distorted_coords;
  int num_cols   = image_size[0]/sample_spacing; // Grab point pairs for the solver at every
  int num_rows   = image_size[1]/sample_spacing; // interval of "sample_spacing"
  int num_coords = num_rows*num_cols;
  undistorted_coords.resize(num_coords);
  distorted_coords.set_size(num_coords*2);
  
  int index = 0;
  for (int r=0; r<num_rows; ++r) {
    int row = r*sample_spacing;
    for (int c=0; c<num_cols; ++c) {
      int col = c*sample_spacing;

      // Generate an undistorted/distorted point pair using the input model.
      // - Note that the pixel pairs need to be corrected for pitch here.
      Vector2 undist_pixel = Vector2(col, row) * pixel_pitch;
      Vector2 dist_pixel   = input_distortion->distorted_coordinates(pin_model, undist_pixel);
      //std::cout << "Undist: " << undist_pixel << " ---- Distorted: " << dist_pixel << std::endl;

      undistorted_coords[index    ] = undist_pixel;
      // Pack these points into 1D vector alternating col,row
      distorted_coords  [2*index  ] = dist_pixel[0];
      distorted_coords  [2*index+1] = dist_pixel[1];
      ++index;
    }
  } // End loop through sampled pixels

  // Now solve for a complementary lens distortion scheme

  // Init solver object with the undistorted coordinates
  DistortionOptimizeFunctor<DistModelT, NumModelParamsT> solver_model(pin_model,
                                                                       undistorted_coords);
  int status;
  Vector<double, NumModelParamsT> seed; // Start with all zeros (no distortion)

  // Solve for the best new model params that give us the distorted coordinates from the undistorted coordinates.
  Vector<double, NumModelParamsT> model_params
    = math::levenberg_marquardt(solver_model, seed, distorted_coords, status);
  
  // Check the error
  DistModelT new_model(model_params);
  double diff = 0;
  for (size_t i=0; i<undistorted_coords.size(); ++i) {
    Vector2 distorted        = new_model.distorted_coordinates(pin_model, undistorted_coords[i]);
    Vector2 actual_distorted = Vector2(distorted_coords[2*i], distorted_coords[2*i+1]);
    //std::cout << "New: " << distorted << " ---- Actual: " << actual_distorted << std::endl;
    diff += norm_2(distorted - actual_distorted);
  }
  diff /= static_cast<double>(undistorted_coords.size());
  vw_out() << "Approximated " << lens_name << " distortion model using a model of type " <<
    new_model.name() <<  " with a " << diff << " mean error.\n";

  // If the approximation is not very good, keep the original model and warn
  //  the user that things might take a long time.
  const double MAX_ERROR = 0.3;
  if (diff > MAX_ERROR)
    vw_out() << "Warning: Failed to approximate reverse pinhole lens distortion, using the original (slow) model.\n";
  else {

    pin_model.set_lens_distortion(new_model);

    // This must be invoked here, when we know the image size
    if (new_model.name() == RPCLensDistortion::class_name()) 
      compute_undistortion<RPCLensDistortion, RPCLensDistortion::num_distortion_params>
	(pin_model, image_size, sample_spacing);
    if (new_model.name() == RPCLensDistortion5::class_name()) 
      compute_undistortion<RPCLensDistortion5, RPCLensDistortion5::num_distortion_params>
	(pin_model, image_size, sample_spacing);
    if (new_model.name() == RPCLensDistortion6::class_name()) 
      compute_undistortion<RPCLensDistortion6, RPCLensDistortion6::num_distortion_params>
	(pin_model, image_size, sample_spacing);
  }
  
  return diff;
} // End function update_pinhole_for_fast_point2pixel

// Given two epipolar rectified CAHV camera models and input images,
// get the aligned version of the images suitable for writing to disk and processing.
template <class ImageInT, class ImageOutT>
void get_epipolar_transformed_images(std::string const& left_camera_file,
                                     std::string const& right_camera_file,
                                     boost::shared_ptr<camera::CameraModel> const left_cahv_camera,
                                     boost::shared_ptr<camera::CameraModel> const right_cahv_camera,
                                     ImageInT  const& left_image_in,
                                     ImageInT  const& right_image_in,
                                     ImageOutT      & left_image_out,
                                     ImageOutT      & right_image_out,
                                     ValueEdgeExtension<typename ImageOutT::pixel_type> edge_ext) {

  // In the epipolar alignment case, the "camera_models" function returns the CAHVModel type!
  CAHVModel* left_epipolar_cahv  = dynamic_cast<CAHVModel*>(vw::camera::unadjusted_model(&(*left_cahv_camera )));
  CAHVModel* right_epipolar_cahv = dynamic_cast<CAHVModel*>(vw::camera::unadjusted_model(&(*right_cahv_camera)));
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
    update_pinhole_for_fast_point2pixel<TsaiLensDistortion, TsaiLensDistortion::num_distortion_params>(left_pin,  left_image_size );
    update_pinhole_for_fast_point2pixel<TsaiLensDistortion, TsaiLensDistortion::num_distortion_params>(right_pin, right_image_size);
    
    left_image_out  = camera_transform(left_image_in,  left_pin,  *left_epipolar_cahv,
                                       edge_ext, BilinearInterpolation());
    right_image_out = camera_transform(right_image_in, right_pin, *right_epipolar_cahv,
                                       edge_ext, BilinearInterpolation());

  } else {
    vw_throw(ArgumentErr() << "get_epipolar_transformed_images: unsupported camera file type.\n");
  }

}

/// Adjust a pair of epipolar-aligned cameras so that the input images are fully
/// contained in the transformed images.
inline void 
resize_epipolar_cameras_to_fit(PinholeModel const& cam1,      PinholeModel const& cam2,
                               PinholeModel      & epi_cam1,  PinholeModel      & epi_cam2,
                               Vector2i     const& size1,     Vector2i     const& size2,
                               Vector2i          & epi_size1, Vector2i          & epi_size2) {
  // Get transforms from input images to epipolar images
  CameraTransform<PinholeModel, PinholeModel> in_to_epi1(cam1, epi_cam1);
  CameraTransform<PinholeModel, PinholeModel> in_to_epi2(cam2, epi_cam2);
  
  // Figure out the bbox needed to contain the transformed image
  BBox2 epi_bbox1 = compute_transformed_bbox_fast(size1, in_to_epi1);
  BBox2 epi_bbox2 = compute_transformed_bbox_fast(size2, in_to_epi2);
 
  //std::cout << "epi_bbox1: " << epi_bbox1 << std::endl;
  //std::cout << "epi_bbox2: " << epi_bbox2 << std::endl;

  // Figure out leftmost and uppermost pixel coordinates in resampled images 
  double min_col = std::min(epi_bbox1.min().x(), epi_bbox2.min().x());
  double min_row = std::min(epi_bbox1.min().y(), epi_bbox2.min().y());
  
  // Compute an adjustment of the camera center point (CCD point below focal point)
  //  such that leftmost and uppermost pixels fall at col 0 and row 0 respectively.
  // - Shift the center by the number of pixels converted to physical CCD units.
  // - We can freely adjust the intrinsic portions of the epipolar cameras as long as
  //   we do the same modification to both cameras.
  Vector2 point_offset  = epi_cam1.point_offset();
  Vector2 center_adjust = Vector2(min_col, min_row)*epi_cam1.pixel_pitch();
  
  //std::cout << "Point offset = " << point_offset << std::endl;
  //std::cout << "center_adjust = " << center_adjust << std::endl;

  // Apply the adjustments to the input epipolar cameras and recompute the bounding boxes  
  epi_cam1.set_point_offset(point_offset - center_adjust);
  epi_cam2.set_point_offset(point_offset - center_adjust);
  CameraTransform<PinholeModel, PinholeModel> in_to_epi1_new(cam1, epi_cam1);
  CameraTransform<PinholeModel, PinholeModel> in_to_epi2_new(cam2, epi_cam2);
  epi_bbox1 = compute_transformed_bbox_fast(size1, in_to_epi1_new);
  epi_bbox2 = compute_transformed_bbox_fast(size2, in_to_epi2_new);
  
  //std::cout << "epi_bbox1 NEW: " << epi_bbox1 << std::endl;
  //std::cout << "epi_bbox2 NEW: " << epi_bbox2 << std::endl;
 
  // Just return the size so that the output images are written 
  // from (0,0) all the way to the maximum position which is what the rest of the code assumes.
  epi_size1 = epi_bbox1.max();
  epi_size2 = epi_bbox2.max();
} // End resize_epipolar_cameras_to_fit

// Get aligned, epipolar rectified pinhole images for stereo processing.
  template <class ImageInT, class ImageOutT,  class EdgeT, class InterpT>
void get_epipolar_transformed_pinhole_images(std::string const& left_camera_file,
                                             std::string const& right_camera_file,
                                             boost::shared_ptr<camera::CameraModel> const left_camera,
                                             boost::shared_ptr<camera::CameraModel> const right_camera,
                                             ImageInT  const& left_image_in,
                                             ImageInT  const& right_image_in,
                                             Vector2i  const& left_image_out_size,
                                             Vector2i  const& right_image_out_size,
                                             ImageOutT      & left_image_out,
                                             ImageOutT      & right_image_out,
                                             EdgeT const& edge_func, InterpT const& interp_func){

  // Cast input camera models to Pinhole
  PinholeModel* left_epipolar_pin  = dynamic_cast<PinholeModel*>(vw::camera::unadjusted_model(&(*left_camera )));
  PinholeModel* right_epipolar_pin = dynamic_cast<PinholeModel*>(vw::camera::unadjusted_model(&(*right_camera)));

  // Load the original camera models with distortion
  PinholeModel left_pin (left_camera_file );
  PinholeModel right_pin(right_camera_file);

  // If necessary, replace the lens distortion models with a approximated models
  //  that will be much faster in the camera_transform calls below.
  Vector2i left_image_in_size (left_image_in.cols(),  left_image_in.rows() );
  Vector2i right_image_in_size(right_image_in.cols(), right_image_in.rows());
  update_pinhole_for_fast_point2pixel<TsaiLensDistortion, TsaiLensDistortion::num_distortion_params>(left_pin,  left_image_in_size );
  update_pinhole_for_fast_point2pixel<TsaiLensDistortion, TsaiLensDistortion::num_distortion_params>(right_pin, right_image_in_size);
  
  // Transform the images
  left_image_out  = camera_transform(left_image_in,  left_pin,  *left_epipolar_pin,
                                     left_image_out_size,
                                     edge_func, interp_func);
  right_image_out = camera_transform(right_image_in, right_pin, *right_epipolar_pin,
                                     right_image_out_size,
                                     edge_func, interp_func);
}


}}      // namespace vw::camera

#endif  //__CAMERAMODEL_CAMERAUTILITIES_H__
