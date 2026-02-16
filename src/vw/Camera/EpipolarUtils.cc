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

/// \file EpipolarUtils.cc
///
/// Utilities for epipolar rectification of stereo image pairs.

#include <vw/Camera/EpipolarUtils.h>
#include <vw/Camera/CameraUtilities.h>
#include <vw/Camera/CameraTransform.h>
#include <vw/Camera/CAHVModel.h>
#include <vw/Camera/CAHVORModel.h>
#include <vw/Camera/CAHVOREModel.h>
#include <vw/Camera/Extrinsics.h>
#include <vw/Image/Transform.h>
#include <vw/Image/Interpolation.h>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace vw {
namespace camera {

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
}

// Apply the epipolar alignment to images
void epipolar_transformed_images(
         std::string const& left_camera_file,
         std::string const& right_camera_file,
         boost::shared_ptr<CameraModel> const left_epi_cam,
         boost::shared_ptr<CameraModel> const right_epi_cam,
         ImageViewRef<PixelMask<float>> const& left_image_in,
         ImageViewRef<PixelMask<float>> const& right_image_in,
         BBox2i const& left_image_in_roi,
         BBox2i const& right_image_in_roi,
         Vector2i     & left_image_out_size,
         Vector2i     & right_image_out_size,
         ImageViewRef<PixelMask<float>> & left_image_out,
         ImageViewRef<PixelMask<float>> & right_image_out,
         ValueEdgeExtension<PixelMask<float>> const& edge_ext) {

  auto interp_func = BilinearInterpolation();

  // Load the original camera models with distortion
  boost::shared_ptr<PinholeModel> left_pin (new PinholeModel(left_camera_file));
  boost::shared_ptr<PinholeModel> right_pin(new PinholeModel(right_camera_file));

  // If necessary, replace the lens distortion models with a approximated models
  //  that will be much faster in the camera_transform calls below.
  Vector2i left_image_in_size (left_image_in.cols(),  left_image_in.rows());
  Vector2i right_image_in_size(right_image_in.cols(), right_image_in.rows());
  bool force_conversion = false;
  *left_pin.get() = fitPinholeModel(left_pin.get(), left_image_in_size,
                                    "TsaiLensDistortion", force_conversion);
  *right_pin.get() = fitPinholeModel(right_pin.get(), right_image_in_size,
                                     "TsaiLensDistortion", force_conversion);

  // Make sure there are no adjustments on the aligned camera models
  PinholeModel* left_epi_pin  = dynamic_cast<PinholeModel*>(unadjusted_model(&(*left_epi_cam)));
  PinholeModel* right_epi_pin = dynamic_cast<PinholeModel*>(unadjusted_model(&(*right_epi_cam)));

  // If so, create an adjusted version of the input files from disk and then transform,
  // otherwise just transform.
  Vector3 ZERO_TRANSLATION(0,0,0);
  Quat    ZERO_ROTATION(1,0,0,0);
  double  NO_SCALE = 1.0;
  if (left_image_in_roi.min() != Vector2i(0,0)) {
    // Get the input raw camera model with the crop applied.
    AdjustedCameraModel left_adj_cam(left_pin, ZERO_TRANSLATION, ZERO_ROTATION,
                                     left_image_in_roi.min(), NO_SCALE);

    // Convert from the cropped input to the epipolar-aligned camera model.
    left_image_out = camera_transform(left_image_in, left_adj_cam, *left_epi_pin,
                                      left_image_out_size, edge_ext, interp_func);
  } else {
    // Don't need to handle crops.
    left_image_out = camera_transform(left_image_in, *left_pin, *left_epi_pin,
                                      left_image_out_size, edge_ext, interp_func);
  }
  // Repeat for the right camera.
  if (right_image_in_roi.min() != Vector2i(0,0)) {
    AdjustedCameraModel right_adj_cam(right_pin, ZERO_TRANSLATION, ZERO_ROTATION,
                                      right_image_in_roi.min(), NO_SCALE);

    right_image_out = camera_transform(right_image_in, right_adj_cam, *right_epi_pin,
                                       right_image_out_size, edge_ext, interp_func);
  } else {
    // No crops
    right_image_out = camera_transform(right_image_in, *right_pin, *right_epi_pin,
                                       right_image_out_size, edge_ext, interp_func);
  }
}

// Given two epipolar rectified CAHV camera models and input images,
// get the aligned version of the images suitable for writing to disk and processing.
void epipolar_transformed_cahv_images(
        std::string const& left_camera_file,
        std::string const& right_camera_file,
        boost::shared_ptr<CameraModel> const left_cahv_camera,
        boost::shared_ptr<CameraModel> const right_cahv_camera,
        ImageViewRef<PixelMask<float>> const& left_image_in,
        ImageViewRef<PixelMask<float>> const& right_image_in,
        ImageViewRef<PixelMask<float>> & left_image_out,
        ImageViewRef<PixelMask<float>> & right_image_out,
        ValueEdgeExtension<PixelMask<float>> const& edge_ext) {

  auto interp_func = BilinearInterpolation();

  // In the epipolar alignment case, the "camera_models" function returns the CAHVModel type!
  CAHVModel* left_epipolar_cahv  = dynamic_cast<CAHVModel*>(unadjusted_model(&(*left_cahv_camera)));
  CAHVModel* right_epipolar_cahv = dynamic_cast<CAHVModel*>(unadjusted_model(&(*right_cahv_camera)));
  if (!left_epipolar_cahv || !right_epipolar_cahv)
    vw_throw(ArgumentErr() << "load_cahv_pinhole_camera_model: CAHVModel cast failed!\n");

  // Lens distortion is corrected in the transform from the input camera files
  //  to the input distortionless CAHV camera models.

  std::string lcase_file = boost::to_lower_copy(left_camera_file);

  // Remove lens distortion and create epipolar rectified images.
  if (boost::ends_with(lcase_file, ".cahvore")) {
    CAHVOREModel left_cahvore (left_camera_file);
    CAHVOREModel right_cahvore(right_camera_file);
    left_image_out  = camera_transform(left_image_in,  left_cahvore,  *left_epipolar_cahv,
                                       edge_ext, interp_func);
    right_image_out = camera_transform(right_image_in, right_cahvore, *right_epipolar_cahv,
                                       edge_ext, interp_func);

  } else if (boost::ends_with(lcase_file, ".cahvor") ||
             boost::ends_with(lcase_file, ".cmod")) {
    CAHVORModel left_cahvor (left_camera_file);
    CAHVORModel right_cahvor(right_camera_file);
    left_image_out  = camera_transform(left_image_in,  left_cahvor,  *left_epipolar_cahv,
                                       edge_ext, interp_func);
    right_image_out = camera_transform(right_image_in, right_cahvor, *right_epipolar_cahv,
                                       edge_ext, interp_func);

  } else if (boost::ends_with(lcase_file, ".cahv") ||
              boost::ends_with(lcase_file, ".pin")) {
    CAHVModel left_cahv (left_camera_file);
    CAHVModel right_cahv(right_camera_file);
    left_image_out  = camera_transform(left_image_in,  left_cahv,  *left_epipolar_cahv,
                                       edge_ext, interp_func);
    right_image_out = camera_transform(right_image_in, right_cahv, *right_epipolar_cahv,
                                       edge_ext, interp_func);

  } else if (boost::ends_with(lcase_file, ".pinhole") ||
              boost::ends_with(lcase_file, ".tsai")) {
    PinholeModel left_pin (left_camera_file);
    PinholeModel right_pin(right_camera_file);

    // If necessary, replace the lens distortion models with a approximated models
    //  that will be much faster in the camera_transform calls below.
    Vector2i left_image_size (left_image_in.cols(),  left_image_in.rows());
    Vector2i right_image_size(right_image_in.cols(), right_image_in.rows());
    bool force_conversion = false;
    left_pin = fitPinholeModel(&left_pin, left_image_size, "TsaiLensDistortion",
                               force_conversion);
    right_pin = fitPinholeModel(&right_pin, right_image_size, "TsaiLensDistortion",
                                force_conversion);

    left_image_out  = camera_transform(left_image_in,  left_pin,  *left_epipolar_cahv,
                                       edge_ext, interp_func);
    right_image_out = camera_transform(right_image_in, right_pin, *right_epipolar_cahv,
                                       edge_ext, interp_func);

  } else {
    vw_throw(ArgumentErr() << "epipolar_transformed_cahv_images: unsupported camera file type.\n");
  }
}

}} // end namespace vw::camera
