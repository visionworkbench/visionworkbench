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

// Apply the epipolar alignment to images
void epipolar_transformed_images(
         std::string const& left_camera_file,
         std::string const& right_camera_file,
         boost::shared_ptr<CameraModel> const left_epi_cam,
         boost::shared_ptr<CameraModel> const right_epi_cam,
         vw::ImageViewRef<vw::PixelMask<float>> const& left_image_in,
         vw::ImageViewRef<vw::PixelMask<float>> const& right_image_in,
         vw::BBox2i const& left_image_in_roi,
         vw::BBox2i const& right_image_in_roi,
         vw::Vector2i       & left_image_out_size,
         vw::Vector2i       & right_image_out_size,
         vw::ImageViewRef<vw::PixelMask<float>> & left_image_out,
         vw::ImageViewRef<vw::PixelMask<float>> & right_image_out,
         vw::ValueEdgeExtension<vw::PixelMask<float>> const& edge_ext);

// Given two epipolar rectified CAHV camera models and input images,
// get the aligned version of the images suitable for writing to disk and processing.
void epipolar_transformed_cahv_images(
        std::string const& left_camera_file,
        std::string const& right_camera_file,
        boost::shared_ptr<CameraModel> const left_cahv_camera,
        boost::shared_ptr<CameraModel> const right_cahv_camera,
        vw::ImageViewRef<vw::PixelMask<float>> const& left_image_in,
        vw::ImageViewRef<vw::PixelMask<float>> const& right_image_in,
        vw::ImageViewRef<vw::PixelMask<float>> & left_image_out,
        vw::ImageViewRef<vw::PixelMask<float>> & right_image_out,
        vw::ValueEdgeExtension<vw::PixelMask<float>> const& edge_ext);

}}      // namespace vw::camera

#endif  //__CAMERAMODEL_CAMERAUTILITIES_H__
