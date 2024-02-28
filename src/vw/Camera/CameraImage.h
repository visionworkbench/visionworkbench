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

#ifndef __VW_CAMERA_CAMERAIMAGE_H__
#define __VW_CAMERA_CAMERAIMAGE_H__

/// \file CameraIamge.h
/// Functionality for handling images and cameras 

#include <vw/Image/ImageViewRef.h>
#include <vw/Camera/CameraModel.h>

#include <boost/shared_ptr.hpp>

namespace vw { namespace camera {
                
// Estimate the GSD at the given pixel given an estimate of the ground point
// along the ray (or close to the ray). This can throw exceptions.
double estimatedGSD(vw::camera::CameraModel const* camera_model, 
                    vw::BBox2i const& image_bbox,
                    vw::Vector2 const& pixel,
                    vw::Vector3 const& position);                
                
}} // end namespace vw::camera

#endif // __VW_CAMERA_CAMERAIMAGE_H__

