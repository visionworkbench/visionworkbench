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

/// Load a pinhole camera model of any supported type
boost::shared_ptr<vw::camera::CameraModel> load_pinhole_camera_model(std::string const& path);

/// Load a pinhole, CAHV, CAHVOR, or CAHVORE model and convert to CAHV.
boost::shared_ptr<vw::camera::CAHVModel>
load_cahv_pinhole_camera_model(std::string const& image_path,
                               std::string const& camera_path);

/// Compute a good sample spacing for the given input image so that 
/// the undistortion functions below won't take too long.
int auto_compute_sample_spacing(Vector2i const image_size);
 
// For RPC, must always ensure the undistortion coefficients are up to date.
// This function will be used before saving or displaying a pinhole model.
// This function is not in the .cc file as it is related to the above.
void update_rpc_undistortion(PinholeModel const& model);


/// Adjust a pair of epipolar-aligned cameras so that the input images are fully
/// contained in the transformed images.
void resize_epipolar_cameras_to_fit(PinholeModel const& cam1,      PinholeModel const& cam2,
                                    PinholeModel      & epi_cam1,  PinholeModel      & epi_cam2,
                                    BBox2i       const& roi1,      BBox2i       const& roi2,
                                    Vector2i          & epi_size1, Vector2i          & epi_size2);


// These template functions are defined inline:


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

/// Class to solve for undistortion coefficients. Applicable only for RPCLensDistortion.
template <class DistModelT>
struct UndistortionOptimizeFunctor:
    public math::LeastSquaresModelBase< UndistortionOptimizeFunctor<DistModelT> > {
  typedef Vector<double> result_type;
  typedef Vector<double> domain_type;
  typedef Matrix<double> jacobian_type;

  const camera::PinholeModel& m_cam;
  const std::vector<Vector2>& m_dist_coords;
  
  /// Init the object with the pinhole model and the list of undistorted image coordinates.
  /// - The list of distorted image coordinates in a Vector<double> (packed alternating col, row) 
  ///    must be passed to the solver function.
  UndistortionOptimizeFunctor(const camera::PinholeModel& cam,
                              const std::vector<Vector2>& dist_coords):
    m_cam(cam), m_dist_coords(dist_coords) {}
  
  /// Return a Vector<double> of all the distorted pixel coordinates.
  inline result_type operator()( domain_type const& x ) const {
    DistModelT lens;
    lens.set_distortion_parameters(x); // Won't be used, but is required by the API
    lens.set_undistortion_parameters(x); // Construct lens distortion model with given parameters
    result_type out_vec;
    
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
// TODO: Move this to the class RPCLensDistortion.
template<class DistModelT>
double compute_undistortion(PinholeModel& pin_model, Vector2i image_size, int sample_spacing) {

  // Get info on existing distortion model
  const vw::camera::LensDistortion* input_distortion = pin_model.lens_distortion();
  std::string lens_name = input_distortion->name();

  // Check for all of the models that currently support a fast distortion function.
  // - The other models use a solver for this function, greatly increasing the run time.
  if (lens_name != RPCLensDistortion::class_name())  
    vw_throw(ArgumentErr() << "Undistortion can only be computed for RPC models.\n");
  
  // Get input camera information
  const double pixel_pitch = pin_model.pixel_pitch();
  // Generate a set of point pairs between undistorted and distorted coordinates.
  std::vector<Vector2> distorted_coords;
  Vector<double> undistorted_coords;
#if 0
  int num_cols   = image_size[0]/sample_spacing; // Grab point pairs for the solver at every
  int num_rows   = image_size[1]/sample_spacing; // interval of "sample_spacing"
  num_cols = std::max(num_cols, 2);
  num_rows = std::max(num_rows, 2);
#else
  int num_cols   = image_size[0]/sample_spacing; // Grab point pairs for the solver at every
  int num_rows   = image_size[1]/sample_spacing; // interval of "sample_spacing"
#endif
  int num_coords = num_rows*num_cols;
  
  distorted_coords.resize(num_coords);
  undistorted_coords.set_size(num_coords*2);

  // TODO: Need to get boxes of where the RPC model is valid.
  // TODO: The logic below will break down for huge distortion.
  int index = 0;
#if 0
  for (int r = 0; r <= num_rows - 1; r++) {
    double row = (image_size[0] - 1) * double(r)/(num_rows - 1); // sample carefully

    for (int c = 0; c <= num_cols - 1; c++) {
      double col = (image_size[1] - 1) * double(c)/(num_cols - 1); // sample carefully
#else
  for (int r = 0; r < num_rows; r++) {

    // TODO: This does not sample well close to num_rows
    int row = r*sample_spacing;
    
    for (int c = 0; c < num_cols; c++) {

      // TODO: This does not sample well close to num_cols
      int col = c*sample_spacing;
#endif      
      // TODO: Need more care below to fill up the domain of validity.
      // Generate an distorted/undistorted point pair using the input model.
      // Note that the pixel pairs need to be corrected for pitch here.
      Vector2 undist_pixel = Vector2(col, row) * pixel_pitch;
      Vector2 dist_pixel   = input_distortion->distorted_coordinates(pin_model, undist_pixel);

      distorted_coords[index] = dist_pixel;
      // Pack these points into 1D vector alternating col, row      
      undistorted_coords[2*index  ] = undist_pixel[0]; 
      undistorted_coords[2*index+1] = undist_pixel[1];
      ++index;
    }
  } // End loop through sampled pixels

  
  // Now solve for a complementary lens distortion scheme

  // Init solver object with the distorted coordinates
  UndistortionOptimizeFunctor<DistModelT> solver_model(pin_model, distorted_coords);
  int status;
  Vector<double> seed;
  seed.set_size(input_distortion->num_dist_params());
  if (DistModelT::class_name() != RPCLensDistortion::class_name()) 
    seed.set_all(0); // Start with all zeros (no distortion)
  else
    RPCLensDistortion::init_as_identity(seed); 

  // Solve for the best new model params that give us the undistorted
  // coordinates from the distorted coordinates.
  Vector<double> undist_params = math::levenberg_marquardt(solver_model,
                                                           seed, undistorted_coords, status);

  // Make a copy of the model
  boost::shared_ptr<LensDistortion> dist_copy = input_distortion->copy();
  DistModelT* new_model = dynamic_cast<DistModelT*>(dist_copy.get()); // Child class access to dist_copy
  if (new_model == NULL) 
    vw_throw( ArgumentErr() << "PinholeModel::expecting an " +
              DistModelT::class_name() + " model." );

  // Set the new undistortion params
  new_model->set_undistortion_parameters(undist_params);
  new_model->set_image_size(image_size);
  
  // Check the error
  double diff = 0;
  for (size_t i=0; i < distorted_coords.size(); ++i) {
    Vector2 undistorted        = new_model->undistorted_coordinates(pin_model, distorted_coords[i]);
    Vector2 actual_undistorted = Vector2(undistorted_coords[2*i], undistorted_coords[2*i+1]);
    diff += norm_2(undistorted - actual_undistorted);
  }
  diff /= static_cast<double>(distorted_coords.size());
  diff /= pixel_pitch; // convert to pixels
  vw_out() << DistModelT::class_name() << " undistortion approximation mean pixel error: "
           << diff << ".\n";
  
  pin_model.set_lens_distortion(new_model); // Set updated distortion model
  return diff;
} 

///  Given a camera model (pinhole or optical bar), create an approximate pinhole model
/// of the desired type.
template<class DistModelT>
double create_approx_pinhole_model(CameraModel * const input_model,
                                   PinholeModel& out_model, Vector2i image_size,
                                   int sample_spacing, bool force_conversion,
                                   int rpc_degree) {

  double pixel_pitch;
  std::string lens_name;
  OpticalBarModel * opb_model = dynamic_cast<OpticalBarModel*>(input_model);
  PinholeModel    * pin_model = dynamic_cast<PinholeModel*>(input_model);

  // Handle the case when the input model is Pinhole
  if (pin_model != NULL) {
    
    out_model   = *pin_model;
    pixel_pitch = pin_model->pixel_pitch();
    
    // Get info on existing distortion model
    lens_name = pin_model->lens_distortion()->name();
    
    // Check for all of the models that currently support a fast distortion function.
    // - The other models use a solver for this function, greatly increasing the run time.
    if ( (!force_conversion) && 
         (lens_name == "NULL"                           ||
          lens_name == "TSAI"                           ||
          lens_name == "AdjustableTSAI"                 ||
          lens_name == RPCLensDistortion::class_name())) {
      //vw_out() << "Input distortion is: " << lens_name << ". Refusing to run.\n";
      return 0;
    }
  }

  // Handle the case when the input model is OpticalBarModel
  if (opb_model != NULL) {
    
    pixel_pitch = opb_model->get_pixel_size();
    double  f = opb_model->get_focal_length(); 
    Vector2 c = opb_model->get_optical_center();
    c *= pixel_pitch; // because OpticalBarModel optical center is in pixels

    // Cook up a pinhole model copying the optical center, focal length,
    // pitch, and orientation. This has no distortion, on purpose.
    Vector2i image_size = opb_model->get_image_size();
    out_model = PinholeModel(opb_model->camera_center(image_size/2.0),
                             opb_model->camera_pose(image_size/2.0).rotation_matrix(),
                             f, f, c[0], c[1]);
    out_model.set_pixel_pitch(pixel_pitch);
    lens_name = "OpticalBar";
  }
  
  // We will try to approximate the input lens with a lens model of type DistModelT.
  
  // Generate a set of point pairs
  std::vector<Vector2> undistorted_coords;
  Vector<double> distorted_coords;
  int num_cols   = image_size[0]/sample_spacing; // Grab point pairs for the solver at every
  int num_rows   = image_size[1]/sample_spacing; // interval of "sample_spacing"
  int num_coords = num_rows*num_cols;
  undistorted_coords.resize(num_coords);
  distorted_coords.set_size(num_coords*2);
  
  int index = 0;
  for (int r = 0; r < num_rows; r++) {

    // TODO: This does not sample well close to num_rows
    int row = r*sample_spacing;
    
    for (int c = 0; c < num_cols; c++) {

      // TODO: This does not sample well close to num_cols
      int col = c*sample_spacing;
      
      // Generate an undistorted/distorted point pair using the input model.
      // - Note that the pixel pairs need to be corrected for pitch here.
      Vector2 pix(col, row);
      Vector2 undist_pixel, dist_pixel;

      // Compute the distorted pixel for the pinhole camera
      if (pin_model != NULL) {

        undist_pixel = pix * pixel_pitch;
        dist_pixel   = pin_model->lens_distortion()->distorted_coordinates(*pin_model, undist_pixel);
        
      }else if (opb_model != NULL) {

        // We will sample the image, which has distorted pixels. Then
        // we will determine the undistorted pixels.
        // The same logic may need to be used for pinhole too.
        
        dist_pixel = pix * pixel_pitch; // Distortion is in units of pixel_pitch

        // Find a point on the ray emanating from the current pixel
        double dist_to_ground = 100000; // 100 km 
        Vector3 xyz = opb_model->camera_center(pix) + dist_to_ground*opb_model->pixel_to_vector(pix);

        // Project into the output pinhole model. For now, it has no distortion.
        undist_pixel = out_model.point_to_pixel(xyz);
        undist_pixel *= pixel_pitch; // Undistortion is in units of pixel pitch
      }
      
      undistorted_coords[index] = undist_pixel;
      // Pack these points into 1D vector alternating col,row
      distorted_coords  [2*index]   = dist_pixel[0];
      distorted_coords  [2*index+1] = dist_pixel[1];
      ++index;
    }
  } // End loop through sampled pixels

  // Now solve for a complementary lens distortion scheme

  // Init solver object with the undistorted coordinates
  DistortionOptimizeFunctor<DistModelT> solver_model(out_model, undistorted_coords);
  int status;
  Vector<double> seed; // Start with all zeros (no distortion)

  // For AdjustableTSAI and RPC, the number of distortion parameters is not fixed.
  // We don't handle AdjustableTSAI at all here.
  int num_distortion_params = 0;
  DistModelT init_model;
  if (DistModelT::class_name() == "AdjustableTSAI")
    vw_throw(ArgumentErr() << "Cannot create an AdjustableTSAI pinhole model.\n");
  else if (DistModelT::class_name() != RPCLensDistortion::class_name()) {
    num_distortion_params = init_model.num_dist_params();
    seed.set_size(num_distortion_params);
    seed.set_all(0);
  } else{
    num_distortion_params = RPCLensDistortion::num_dist_params(rpc_degree); 
    seed.set_size(num_distortion_params);
    RPCLensDistortion::init_as_identity(seed); 
  }
  
  // Solve for the best new model params that give us the distorted
  // coordinates from the undistorted coordinates.
  Vector<double> model_params = math::levenberg_marquardt(solver_model, seed,
                                                          distorted_coords, status);
  
  // Check the error
  DistModelT new_model(model_params);
  double diff = 0;
  for (size_t i = 0; i < undistorted_coords.size(); ++i) {
    Vector2 distorted        = new_model.distorted_coordinates(out_model, undistorted_coords[i]);
    Vector2 actual_distorted = Vector2(distorted_coords[2*i], distorted_coords[2*i+1]);
    diff += norm_2(distorted - actual_distorted);
  }
  diff /= static_cast<double>(undistorted_coords.size());
  diff /= pixel_pitch; // convert the errors to pixels
  
  vw_out() << "Approximated an " << lens_name << " distortion model using a model of type " <<
    new_model.name() <<  " with a " << diff << " mean pixel error.\n";

  // If the approximation is not very good, keep the original model and warn
  //  the user that things might take a long time.
  const double MAX_ERROR = 0.3;
  if ( (!force_conversion) && diff > MAX_ERROR)
    vw_out() << "Warning: Failed to approximate reverse pinhole lens distortion, "
             << "using the original (slow) model.\n";
  else {

    out_model.set_lens_distortion(&new_model);

    // This must be invoked here, when we know the image size
    if (new_model.name() == RPCLensDistortion::class_name()) 
      compute_undistortion<RPCLensDistortion>(out_model, image_size, sample_spacing);
  }
  
  return diff;
} 

/// If necessary, replace the lens distortion model in the pinhole camera model
///  with an approximated model that has a fast distortion function needed for
///  quick computation of the point_to_pixel function.
/// - Does not replace the camera if the approximation error is too high.
/// - Returns the approximation error.
template<class DistModelT>
double update_pinhole_for_fast_point2pixel(PinholeModel& pin_model, Vector2i image_size,
                                           int sample_spacing = 0,
                                           bool force_conversion = false) {
  
  if (sample_spacing <= 0) 
    sample_spacing = auto_compute_sample_spacing(image_size);
  
  PinholeModel in_model = pin_model;
  int rpc_degree = 0;
  return create_approx_pinhole_model<DistModelT>
    (&in_model, pin_model, image_size, sample_spacing, force_conversion, rpc_degree);
}

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
    update_pinhole_for_fast_point2pixel<TsaiLensDistortion>(left_pin,  left_image_size );
    update_pinhole_for_fast_point2pixel<TsaiLensDistortion>(right_pin, right_image_size);
    
    left_image_out  = camera_transform(left_image_in,  left_pin,  *left_epipolar_cahv,
                                       edge_ext, BilinearInterpolation());
    right_image_out = camera_transform(right_image_in, right_pin, *right_epipolar_cahv,
                                       edge_ext, BilinearInterpolation());

  } else {
    vw_throw(ArgumentErr() << "get_epipolar_transformed_images: unsupported camera file type.\n");
  }

}



template <class ImageInT, class ImageOutT,  class EdgeT, class InterpT>
void get_epipolar_transformed_pinhole_images(std::string const& left_camera_file,
                                             std::string const& right_camera_file,
                                             boost::shared_ptr<camera::CameraModel> const left_epi_cam,
                                             boost::shared_ptr<camera::CameraModel> const right_epi_cam,
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
  update_pinhole_for_fast_point2pixel<TsaiLensDistortion>(*left_pin,  left_image_in_size );
  update_pinhole_for_fast_point2pixel<TsaiLensDistortion>(*right_pin, right_image_in_size);
 
  // Make sure there are no adjustments on the aligned camera models
  PinholeModel* left_epi_pin  = dynamic_cast<PinholeModel*>(vw::camera::unadjusted_model(&(*left_epi_cam)));
  PinholeModel* right_epi_pin = dynamic_cast<PinholeModel*>(vw::camera::unadjusted_model(&(*right_epi_cam)));
  
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
