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


/// \file CameraUtilities.cc
///
/// This file contains miscellaneous functions for working with camera models.
///
#include <vw/Camera/CameraUtilities.h>

namespace vw {
namespace camera {

/// Load a pinhole camera model of any supported type
boost::shared_ptr<vw::camera::CameraModel>
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
boost::shared_ptr<vw::camera::CAHVModel>
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
    // The CAHV class is constructed from a Pinhole model.
    vw::camera::PinholeModel left_pin(camera_path);
    *(cahv.get()) = vw::camera::strip_lens_distortion(left_pin);

  } else {
    vw_throw(vw::ArgumentErr() << "load_cahv_pinhole_camera_model - unsupported camera file type.\n");
  }

  return cahv;
}

int auto_compute_sample_spacing(Vector2i const image_size) {
  int DEFAULT_SPACING = 50;
  int MIN_SIDE_POINTS = 100;
  int MAX_SIDE_POINTS = 600;
  
  // Use the default spacing unless it gives too many points, then use the min spacing
  //  that hits our point estimate.
  int spacing    = DEFAULT_SPACING;
  int side_total = image_size[0] + image_size[1];
  int num_pts    = side_total / DEFAULT_SPACING;
  
  if (num_pts < MIN_SIDE_POINTS)
    spacing = std::max(side_total / MIN_SIDE_POINTS, 1);
  
  if (num_pts > MAX_SIDE_POINTS)
    spacing = side_total / MAX_SIDE_POINTS;
  
  return spacing;
}

void update_rpc_undistortion(PinholeModel const& model){

  const vw::camera::LensDistortion* distortion = model.lens_distortion();
  std::string lens_name = distortion->name();
  if (lens_name != RPCLensDistortion::class_name())
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
    int sample_spacing = auto_compute_sample_spacing(rpc_dist->image_size());
    if (!rpc_dist->can_undistort()) 
      compute_undistortion<RPCLensDistortion>(*pin_ptr, rpc_dist->image_size(), sample_spacing);
  }

}

void 
resize_epipolar_cameras_to_fit(PinholeModel const& cam1,      PinholeModel const& cam2,
                               PinholeModel      & epi_cam1,  PinholeModel      & epi_cam2,
                               BBox2i       const& roi1,      BBox2i       const& roi2,
                               Vector2i          & epi_size1, Vector2i          & epi_size2) {
  // Get transforms from input images to epipolar images
  CameraTransform<PinholeModel, PinholeModel> in_to_epi1(cam1, epi_cam1);
  CameraTransform<PinholeModel, PinholeModel> in_to_epi2(cam2, epi_cam2);
  
  // Figure out the bbox needed to contain the transformed image
  // - This just uses the ROIs for the two input images so the cameras can be shifted to align with 0,0.
  BBox2 epi_bbox1 = compute_transformed_bbox_fast(roi1, in_to_epi1);
  BBox2 epi_bbox2 = compute_transformed_bbox_fast(roi2, in_to_epi2);

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

  // Apply the adjustments to the input epipolar cameras and recompute the bounding boxes  
  epi_cam1.set_point_offset(point_offset - center_adjust);
  epi_cam2.set_point_offset(point_offset - center_adjust);
  CameraTransform<PinholeModel, PinholeModel> in_to_epi1_new(cam1, epi_cam1);
  CameraTransform<PinholeModel, PinholeModel> in_to_epi2_new(cam2, epi_cam2);
  
  // Recompute the bounding boxes to be sure
  epi_bbox1 = compute_transformed_bbox_fast(roi1, in_to_epi1_new);
  epi_bbox2 = compute_transformed_bbox_fast(roi2, in_to_epi2_new);
 
  // Return the size required to contain all of the transformed image data.
  epi_size1 = epi_bbox1.max();
  epi_size2 = epi_bbox2.max();
} // End resize_epipolar_cameras_to_fit


}} // end namespace vw::camera

