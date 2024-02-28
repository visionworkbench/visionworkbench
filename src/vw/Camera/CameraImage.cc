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

/// \file CameraIamge.cc
/// Functionality for handling images and cameras 

#include <vw/Camera/CameraImage.h>

namespace vw { namespace camera {

// Estimate the GSD at the given pixel given an estimate of the ground point
// along the ray (or close to the ray). This can throw exceptions.
double estimatedGSD(vw::camera::CameraModel const* camera_model, 
                    vw::BBox2i const& image_bbox,
                    vw::Vector2 const& pixel,
                    vw::Vector3 const& position) {

  // The image bbox must be at least two pixels wide and tall
  if (image_bbox.width() < 2 || image_bbox.height() < 2)
    vw_throw(ArgumentErr() << "The image dimensions must be at least 2 x 2 pixels.\n");
  
  vw::Vector3 ctr = camera_model->camera_center(pixel);
  vw::Vector3 dir = camera_model->pixel_to_vector(pixel);
  dir = dir/norm_2(dir); // Likely not necessary
  
  // Distance from camera to ground
  double dist = norm_2(position - ctr);
  
  // Find the 3D point along the ray at the given distance
  vw::Vector3 ground = ctr + dist * dir;

  // Do the same for two neighboring pixels. Adjust them to 
  // ensure they are inside the image. Need to do it this way
  // because gsd along row and column may be different.
  vw::Vector2 pixel2 = pixel + vw::Vector2(1, 0);
  if (!image_bbox.contains(pixel2))
    pixel2 = pixel + vw::Vector2(-1, 0);
  vw::Vector2 pixel3 = pixel + vw::Vector2(0, 1);
  if (!image_bbox.contains(pixel3))
    pixel3 = pixel + vw::Vector2(0, -1);
  
  vw::Vector3 dir2 = camera_model->pixel_to_vector(pixel2);
  dir2 = dir2/norm_2(dir2); // Normalize the direction vector
  vw::Vector3 ground2 = camera_model->camera_center(pixel2) + dist * dir2;
  
  vw::Vector3 dir3 = camera_model->pixel_to_vector(pixel3);
  dir3 = dir3/norm_2(dir3); // Normalize the direction vector
  vw::Vector3 ground3 = camera_model->camera_center(pixel3) + dist * dir3;
  
  double gsd2 = norm_2(ground - ground2);
  double gsd3 = norm_2(ground - ground3);
  double gsd = (gsd2 + gsd3)/2.0;

  // Throw an exception if the GSD is not positive or nan
  if (gsd <= 0 || std::isnan(gsd))
    vw_throw(ArgumentErr() << "The estimated GSD is not positive or is nan.\n");
    
  return gsd;  
}    

}} // end namespace vw::camera
