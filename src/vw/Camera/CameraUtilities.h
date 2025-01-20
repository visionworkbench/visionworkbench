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

///  Given a camera model (pinhole or optical bar), create an approximate pinhole model
/// of the desired type.
template<class DistModelT>
double create_approx_pinhole_model(CameraModel const* input_model,
                                   PinholeModel& out_model, Vector2i image_size,
                                   int sample_spacing, bool force_conversion,
                                   int rpc_degree, double camera_to_ground_dist) {

  std::cout << "--sample spacing before: " << sample_spacing << std::endl;
  if (sample_spacing <= 0) 
    sample_spacing = auto_compute_sample_spacing(image_size);
  std::cout << "--sample spacing after: " << sample_spacing << std::endl;
  
  double pixel_pitch;
  std::string lens_name;
  OpticalBarModel const* opb_model = dynamic_cast<OpticalBarModel const*>(input_model);
  PinholeModel    const* pin_model = dynamic_cast<PinholeModel const*>(input_model);

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
      std::cout << "---tmp reusing to run!\n";
      return 0;
    }
  }

  // Handle the case when the input model is OpticalBarModel
  if (opb_model != NULL) {
    lens_name = "OpticalBar";
    out_model = opticalbar2pinhole(*opb_model, sample_spacing, camera_to_ground_dist);
    pixel_pitch = out_model.pixel_pitch();
  }
  
  // We will try to approximate the input lens with a lens model of type DistModelT.
  
  // Generate a set of point pairs
  std::vector<Vector2> undistorted_coords;
  Vector<double> distorted_coords;
  // Grab point pairs for the solver at every
  // interval of "sample_spacing"
  int num_col_samples = image_size[0]/sample_spacing;
  int num_row_samples = image_size[1]/sample_spacing;
  num_col_samples = std::max(num_col_samples, 2);
  num_row_samples = std::max(num_row_samples, 2);
  int num_coords = num_row_samples*num_col_samples;
  undistorted_coords.resize(num_coords);
  distorted_coords.set_size(num_coords*2);

  std::vector<Vector3> input_xyz;
  std::vector<double> input_pixels;
  
  int index = 0;
  for (int r = 0; r < num_row_samples; r++) {

    // sample the full interval [0, image_size[1]-1] uniformly
    double row = (image_size[1] - 1.0) * double(r)/(num_row_samples - 1.0); 
    
    for (int c = 0; c < num_col_samples; c++) {

      // sample the full interval [0, image_size[0]-1] uniformly
      double col = (image_size[0] - 1.0) * double(c)/(num_col_samples - 1.0);

      // Generate an undistorted/distorted point pair using the input model.
      // - Note that the pixel pairs need to be corrected for pitch here.
      Vector2 pix(col, row);

      input_pixels.push_back(pix[0]);
      input_pixels.push_back(pix[1]);
      
      Vector2 undist_pixel, dist_pixel;

      // Compute the distorted pixel for the pinhole camera
      if (pin_model != NULL) {

        undist_pixel = pix * pixel_pitch;
        dist_pixel 
          = pin_model->lens_distortion()->distorted_coordinates(*pin_model, undist_pixel);
        
      } else if (opb_model != NULL) {

        // We will sample the image, which has distorted pixels. Then
        // we will determine the undistorted pixels.
        // The same logic may need to be used for pinhole too.
        
        dist_pixel = pix * pixel_pitch; // Distortion is in units of pixel_pitch

        Vector3 xyz = opb_model->camera_center(pix)
          + camera_to_ground_dist*opb_model->pixel_to_vector(pix);

        input_xyz.push_back(xyz);
        
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

  int status = -1;
  Vector<double> seed; // Start with all zeros (no distortion)

  // For AdjustableTSAI and RPC, the number of distortion parameters is not fixed.
  // We don't handle AdjustableTSAI at all here.
  int num_distortion_params = 0;
  DistModelT init_model;
  bool do_rpc = (DistModelT::class_name() == RPCLensDistortion::class_name());
  int initial_rpc_degree = 1;
  if (DistModelT::class_name() == "AdjustableTSAI")
    vw_throw(ArgumentErr() << "Cannot create an AdjustableTSAI pinhole model.\n");
  else if (!do_rpc) {
    num_distortion_params = init_model.num_dist_params();
    seed.set_size(num_distortion_params);
    seed.set_all(0);
  } else{
    // rpc model
    num_distortion_params = RPCLensDistortion::num_dist_params(initial_rpc_degree); 
    seed.set_size(num_distortion_params);
    RPCLensDistortion::init_as_identity(seed); 
  }

  // When desired to find RPC distortion of given degree, first find the distortion
  // of degree 1, use that as a seed to find the distortion of degree 2, until
  // arriving at the desired degree. This  produces better results than aiming
  // right away for the desired degree, as in that case one gets stuck in a bad
  // local minimum.
  int num_passes = 1;
  if (do_rpc) 
    num_passes = std::max(rpc_degree, 1);
    
  // Solve for the best new model params that give us the distorted
  // coordinates from the undistorted coordinates.
  double mean_error = 0.0, max_error = 0.0, norm = 0.0;
  DistModelT new_model;
  Vector<double> model_params;

  vw_out() << "Approximating a distortion model of type: " << lens_name << "\n";
  if (do_rpc) 
    vw_out() << "Compute the RPC distortion model starting at degree 1 "
             << "and then refine it until reaching degree " << rpc_degree << ".\n";
  
  for (int pass = 1; pass <= num_passes; pass++) {

    if (do_rpc && pass >= 2) {
      // Use the previously solved model as a seed. Increment its degree by adding
      // a new power with a zero coefficient in front.
      seed = model_params;
      RPCLensDistortion::increment_degree(seed);
    }
    
    // Init solver object with the undistorted coordinates
    DistortionOptimizeFunctor<DistModelT> solver_model(out_model, undistorted_coords);

    // Find model_params by doing a best fit
    model_params = math::levenberg_marquardt(solver_model, seed,
                                             distorted_coords, status);
    // Check the error
    new_model = DistModelT(model_params);
    mean_error = 0.0;
    max_error = 0.0;
    for (size_t i = 0; i < undistorted_coords.size(); ++i) {
      Vector2 distorted = new_model.distorted_coordinates(out_model, undistorted_coords[i]);
      Vector2 actual_distorted = Vector2(distorted_coords[2*i], distorted_coords[2*i+1]);
      double error = norm_2(distorted - actual_distorted);
      norm += error * error; 
      mean_error += error;
      if (error > max_error)
        max_error = error;
    }
    mean_error /= static_cast<double>(undistorted_coords.size());
    max_error /= pixel_pitch; // convert the errors to pixels
    mean_error /= pixel_pitch; // convert the errors to pixels
    norm = sqrt(norm) / pixel_pitch; // take the square root and convert to pixels
    
    
    vw_out() << "Using a distortion model of type " << new_model.name();
    
    if (do_rpc)
      vw_out() << " of degree " << pass;
    
    vw_out() << " with mean and max pixel error of "
             << mean_error << ", " << max_error << ".\n";
  }
  
  // If the approximation is not very good, keep the original model and warn
  //  the user that things might take a long time.
  const double MAX_ERROR = 0.3;
  if ( (!force_conversion) && mean_error > MAX_ERROR)
    vw_out() << "Warning: Failed to approximate reverse pinhole lens distortion, "
             << "using the original (slow) model.\n";
  else
    out_model.set_lens_distortion(&new_model);
  
  return mean_error;
} 

// Approximate a given model with a pinhole model with given distortion
PinholeModel fitPinholeModel(CameraModel const* in_model, 
                             vw::Vector2 const& image_size,
                             std::string const& out_distortion_type,
                             bool force_conversion,
                             int sample_spacing = 0,
                             int rpc_degree = 0,
                             double camera_to_ground_dist = 0);

/// If necessary, replace the lens distortion model in the pinhole camera model
///  with an approximated model that has a fast distortion function needed for
///  quick computation of the point_to_pixel function.
/// - Does not replace the camera if the approximation error is too high.
/// - Returns the approximation error.
template<class DistModelT>
double update_pinhole_for_fast_point2pixel(PinholeModel& pin_model, Vector2i image_size,
                                           int sample_spacing = 0,
                                           bool force_conversion = false) {
  
  std::cout << "--now in update_pinhole_for_fast_point2pixel\n";
  
  PinholeModel in_model = pin_model;
  int rpc_degree = 0;
  double camera_to_ground_dist = 0;
  return create_approx_pinhole_model<DistModelT>
    (&in_model, pin_model, image_size, sample_spacing, force_conversion,
     rpc_degree, camera_to_ground_dist);
}

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
  update_pinhole_for_fast_point2pixel<TsaiLensDistortion>(*left_pin,  left_image_in_size );
  update_pinhole_for_fast_point2pixel<TsaiLensDistortion>(*right_pin, right_image_in_size);
 
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
