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

/// \file EpipolarUtils.h
///
/// Utilities for epipolar rectification of stereo image pairs.

#ifndef __VW_CAMERA_EPIPOLARUTILS_H__
#define __VW_CAMERA_EPIPOLARUTILS_H__

#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/EdgeExtension.h>

#include <boost/shared_ptr.hpp>

namespace vw {
namespace camera {

/// Adjust a pair of epipolar-aligned cameras so that the input images are fully
/// contained in the transformed images.
void resize_epipolar_cameras_to_fit(PinholeModel const& cam1,
                                    PinholeModel const& cam2,
                                    PinholeModel      & epi_cam1,
                                    PinholeModel      & epi_cam2,
                                    BBox2i       const& roi1,
                                    BBox2i       const& roi2,
                                    Vector2i          & epi_size1,
                                    Vector2i          & epi_size2);

/// Apply the epipolar alignment to images using pinhole camera models.
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
         ValueEdgeExtension<PixelMask<float>> const& edge_ext);

/// Given two epipolar rectified CAHV camera models and input images,
/// get the aligned version of the images suitable for writing to disk and processing.
void epipolar_transformed_cahv_images(
        std::string const& left_camera_file,
        std::string const& right_camera_file,
        boost::shared_ptr<CameraModel> const left_cahv_camera,
        boost::shared_ptr<CameraModel> const right_cahv_camera,
        ImageViewRef<PixelMask<float>> const& left_image_in,
        ImageViewRef<PixelMask<float>> const& right_image_in,
        ImageViewRef<PixelMask<float>> & left_image_out,
        ImageViewRef<PixelMask<float>> & right_image_out,
        ValueEdgeExtension<PixelMask<float>> const& edge_ext);

}} // namespace vw::camera

#endif // __VW_CAMERA_EPIPOLARUTILS_H__
