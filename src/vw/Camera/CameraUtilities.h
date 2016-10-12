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
#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace vw {
namespace camera {

// TODO: Move code into a .cc file!

/// Load a pinhole camera model of any supported type
inline boost::shared_ptr<vw::camera::CameraModel>
load_pinhole_camera_model(std::string const& path){

  std::string lcase_file = boost::to_lower_copy(path);
  if (boost::ends_with(lcase_file,".cahvore") ) {
    return boost::shared_ptr<vw::camera::CameraModel>( new vw::camera::CAHVOREModel(path) );
  } else if (boost::ends_with(lcase_file,".cahvor") ||
             boost::ends_with(lcase_file,".cmod"  )   ) {
    return boost::shared_ptr<vw::camera::CameraModel>( new vw::camera::CAHVORModel(path) );
  } else if ( boost::ends_with(lcase_file,".cahv") ||
              boost::ends_with(lcase_file,".pin" )   ) {
    return boost::shared_ptr<vw::camera::CameraModel>( new vw::camera::CAHVModel(path) );
  } else if ( boost::ends_with(lcase_file,".pinhole") ||
              boost::ends_with(lcase_file,".tsai"   )   ) {
    return boost::shared_ptr<vw::camera::CameraModel>( new vw::camera::PinholeModel(path) );
  } else {
    vw::vw_throw(vw::ArgumentErr() << "PinholeStereoSession: unsupported camera file type.\n");
  }
}


/// Load a pinhole, CAHV, CAHVOR, or CAHVORE model and convert to CAHV.
inline boost::shared_ptr<vw::camera::CAHVModel>
load_cahv_pinhole_camera_model(std::string const& image_path,
                               std::string const& camera_path){
  // Get the image size
  vw::DiskImageView<float> disk_image(image_path);
  vw::Vector2i image_size(disk_image.cols(), disk_image.rows());

  // Load the appropriate camera model object and if necessary
  // convert it to the CAHVModel type.
  std::string lcase_file = boost::to_lower_copy(camera_path);
  boost::shared_ptr<vw::camera::CAHVModel> cahv(new vw::camera::CAHVModel);
  if (boost::ends_with(lcase_file, ".cahvore") ) {
    vw::camera::CAHVOREModel cahvore(camera_path);
    *(cahv.get()) = vw::camera::linearize_camera(cahvore, image_size, image_size);
  } else if (boost::ends_with(lcase_file, ".cahvor")  ||
             boost::ends_with(lcase_file, ".cmod"  )   ) {
    vw::camera::CAHVORModel cahvor(camera_path);
    *(cahv.get()) = vw::camera::linearize_camera(cahvor, image_size, image_size);

  } else if ( boost::ends_with(lcase_file, ".cahv") ||
              boost::ends_with(lcase_file, ".pin" )) {
    *(cahv.get()) = vw::camera::CAHVModel(camera_path);

  } else if ( boost::ends_with(lcase_file, ".pinhole") ||
              boost::ends_with(lcase_file, ".tsai"   )   ) {
    vw::camera::PinholeModel left_pin(camera_path);
    *(cahv.get()) = vw::camera::strip_lens_distortion(left_pin);

  } else {
    vw_throw(vw::ArgumentErr() << "load_cahv_pinhole_camera_model - unsupported camera file type.\n");
  }

  return cahv;
}


/// Class to use with the LevenbergMarquardt solver to optimize the parameters of a TSAI lens
///  to match the passed in pairs of undistorted (input) and distorted (output) points.
struct TsaiLensOptimizeFunctor :  public math::LeastSquaresModelBase<TsaiLensOptimizeFunctor> {
  typedef Vector<double> result_type;
  typedef Vector4 domain_type;
  typedef Matrix<double> jacobian_type;

  const camera::PinholeModel& m_cam;
  const std::vector<Vector2>& m_raw_coords;
  
  /// Init the object with the pinhole model and the list of undistorted image coordinates.
  /// - The list of distorted image coordinates in a Vector<double> (packed alternating col, row) 
  ///    must be passed to the solver function.
  TsaiLensOptimizeFunctor(const camera::PinholeModel& cam, const std::vector<Vector2>& raw_coords)
    : m_cam(cam), m_raw_coords(raw_coords) {}
  
  /// Return a Vector<double> of all the distorted pixel coordinates.
  inline result_type operator()( domain_type const& x ) const {
    TsaiLensDistortion lens(x); // Construct lens distortion model with given parameters
    result_type        out_vec;
    
    out_vec.set_size(m_raw_coords.size()*2);
    for (size_t i=0; i<m_raw_coords.size(); ++i) {
      Vector2 loc = lens.distorted_coordinates(m_cam, m_raw_coords[i]);
      out_vec[2*i  ] = loc[0]; // The col and row values are packed in successive indices.
      out_vec[2*i+1] = loc[1];
    }
    return out_vec;
  }
}; // End class LensOptimizeFunctor


/// If necessary, replace the lens distortion model in the pinhole camera model
///  with an approximated model that has a fast distortion function needed for
///  quick computation of the point_to_pixel function.
/// - Does not replace the camera if the approximation error is too high.
/// - Returns the approximation error.
inline
double update_pinhole_for_fast_point2pixel(PinholeModel& pin_model, Vector2i image_size,
                                         int sample_spacing=50) {

  // Get info on existing distortion model
  const vw::camera::LensDistortion* input_distortion = pin_model.lens_distortion();
  std::string lens_name = input_distortion->name();
  
  // Check for all of the models that currently support a fast distortion function.
  // - The other models use a solver for this function, greatly increasing the run time.
  if ((lens_name == "NULL") || (lens_name == "TSAI") || (lens_name == "AdjustableTSAI"))
    return 0;

  // Otherwise we will try to approximate the input lens with a TSAI lens model.
  
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
      Vector2 raw_pixel       = Vector2(col, row) * pixel_pitch;
      Vector2 distorted_pixel = input_distortion->distorted_coordinates(pin_model, raw_pixel);
      //std::cout << "Raw: " << raw_pixel << " ---- Distorted: " << distorted_pixel << std::endl;

      undistorted_coords[index    ] = raw_pixel;
      distorted_coords  [2*index  ] = distorted_pixel[0]; // Pack these points into 1D vector alternating col,row
      distorted_coords  [2*index+1] = distorted_pixel[1];
      ++index;
    }
  } // End loop through sampled pixels

  // Now solve for a complementary lens distortion scheme

  // Init solver object with the undistorted coordinates
  TsaiLensOptimizeFunctor solver_model(pin_model, undistorted_coords);
  int status;
  Vector4 seed; // Start with all zeros (no distortion)
  // Solve for the best tsai params that give us the distorted coordinates from the undistorted coordinates.
  Vector4 tsai_params = math::levenberg_marquardt( solver_model, seed, distorted_coords, status);
  
  // Check the error
  TsaiLensDistortion new_tsai(tsai_params);
  double diff = 0;
  for (size_t i=0; i<undistorted_coords.size(); ++i) {
    Vector2 distorted        = new_tsai.distorted_coordinates(pin_model, undistorted_coords[i]);
    Vector2 actual_distorted = Vector2(distorted_coords[2*i], distorted_coords[2*i+1]);
    //std::cout << "New: " << distorted << " ---- Actual: " << actual_distorted << std::endl;
    diff += norm_2(distorted - actual_distorted);
  }
  diff /= static_cast<double>(undistorted_coords.size());
  VW_OUT(InfoMessage, "camera") << "Approximated " << lens_name 
      << " distortion model using a Tsai model with " << diff << " mean error.\n";

  // If the approximation is not very good, keep the original model and warn
  //  the user that things might take a long time.
  const double MAX_ERROR = 0.3;
  if (diff > MAX_ERROR)
    vw_out() << "Warning: Failed to approximate reverse pinhole lens distortion, using the original (slow) model.\n";
  else
    pin_model.set_lens_distortion(new_tsai);
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
                                     ImageOutT      & right_image_out) {

    // In the epipolar alignment case, the "camera_models" function returns the CAHVModel type!
    CAHVModel* left_epipolar_cahv  = dynamic_cast<CAHVModel*>(vw::camera::unadjusted_model(&(*left_cahv_camera )));
    CAHVModel* right_epipolar_cahv = dynamic_cast<CAHVModel*>(vw::camera::unadjusted_model(&(*right_cahv_camera)));
    if (!left_epipolar_cahv || !right_epipolar_cahv) {
      vw_throw(ArgumentErr() << "load_cahv_pinhole_camera_model: CAHVModel cast failed!\n");
    }

    std::string lcase_file = boost::to_lower_copy(left_camera_file);

    // Remove lens distortion and create epipolar rectified images.
    if (boost::ends_with(lcase_file, ".cahvore")) {
      CAHVOREModel left_cahvore (left_camera_file );
      CAHVOREModel right_cahvore(right_camera_file);
      left_image_out  = camera_transform(left_image_in,  left_cahvore,  *left_epipolar_cahv );
      right_image_out = camera_transform(right_image_in, right_cahvore, *right_epipolar_cahv);
      
    } else if (boost::ends_with(lcase_file, ".cahvor") ||
               boost::ends_with(lcase_file, ".cmod") ) {
      CAHVORModel left_cahvor (left_camera_file );
      CAHVORModel right_cahvor(right_camera_file);
      left_image_out  = camera_transform(left_image_in,  left_cahvor,  *left_epipolar_cahv );
      right_image_out = camera_transform(right_image_in, right_cahvor, *right_epipolar_cahv);

    } else if ( boost::ends_with(lcase_file, ".cahv") ||
                boost::ends_with(lcase_file, ".pin" )) {
      CAHVModel left_cahv (left_camera_file );
      CAHVModel right_cahv(right_camera_file);
      left_image_out  = camera_transform(left_image_in,  left_cahv,  *left_epipolar_cahv );
      right_image_out = camera_transform(right_image_in, right_cahv, *right_epipolar_cahv);

    } else if ( boost::ends_with(lcase_file, ".pinhole") ||
                boost::ends_with(lcase_file, ".tsai") ) {
      PinholeModel left_pin (left_camera_file );
      PinholeModel right_pin(right_camera_file);

      // If necessary, replace the lens distortion models with a approximated models
      //  that will be much faster in the camera_transform calls below.
      Vector2i left_image_size (left_image_in.cols(),  left_image_in.rows() );
      Vector2i right_image_size(right_image_in.cols(), right_image_in.rows());
      update_pinhole_for_fast_point2pixel(left_pin,  left_image_size );
      update_pinhole_for_fast_point2pixel(right_pin, right_image_size);
      
      left_image_out  = camera_transform(left_image_in,  left_pin,  *left_epipolar_cahv );
      right_image_out = camera_transform(right_image_in, right_pin, *right_epipolar_cahv);

    } else {
      vw_throw(ArgumentErr() << "load_cahv_pinhole_camera_model: unsupported camera file type.\n");
    }

}



}}      // namespace vw::camera

#endif  //__CAMERAMODEL_CAMERAUTILITIES_H__
