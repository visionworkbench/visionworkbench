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

#include <vw/Camera/CameraModel.h>
#include <vw/Camera/LinescanErr.h>

#include <iomanip>

namespace vw {
namespace camera {

// Find two vectors that are perpendicular to each other and to the input unit
// vector.
void findPerpVecs(vw::Vector3 const& vec,
                  vw::Vector3 & perp1, 
                  vw::Vector3 & perp2) {

   // The input vec must have norm of 1, with tolerance
   if (std::abs(norm_2(vec) - 1.0) > 1e-5)
     vw_throw(ArgumentErr() << "findPerpVecs: Input vector must be a unit vector.\n");
     
   // Find the smallest coordinate in vec
   int min_i = 0;
   for (int i = 1; i < 3; i++) {
     if (std::abs(vec[i]) < std::abs(vec[min_i]))
       min_i = i;
   }
   
   // Find the other two indices
   int j = 0, k = 0;
   if (min_i == 0) {
     j = 1; k = 2;
   } else if (min_i == 1) {
     j = 0; k = 2;
   } else { // min_i == 2
     j = 0; k = 1;
   }
   
   // Find the vector that swaps j and k and keeps min_i at 0. This is a vector
   // that is perpendicular to vec.
   perp1 = Vector3(0, 0, 0);
   perp1[min_i] = 0; // Keep the smallest coordinate at 0
   perp1[j] = -vec[k]; // Set the other two coordinates to be perpendicular
   perp1[k] =  vec[j]; // This is the swap, so it is perpendicular to vec

   // Normalize 
   perp1 = normalize(perp1);
    
  // The second perpendicular vector is the cross product of the two
  perp2 = cross_prod(vec, perp1);
  // Normalize 
  perp2 = normalize(perp2);
}

LinescanErr::LinescanErr(const CameraModel* model, 
                         const vw::Vector3& pt, 
                         vw::Vector2 const& guess):
  m_model(model), m_point(pt) {
  
  // Find a direction from the camera to the ground
  vw::Vector3 cam_ctr = m_model->camera_center(guess);  
  vw::Vector3 ground_dir = normalize(m_point - cam_ctr);
  findPerpVecs(ground_dir, m_perp1, m_perp2);
}

// This must have the signature expected by Newton's method. Can throw
// exceptions. The math is described above.
vw::Vector2 LinescanErr::operator()(vw::Vector2 const& pix) const {
  
  vw::Vector3 cam_ctr = m_model->camera_center(pix);
  
  // Normalized direction from camera to ground point
  double dist_to_ground = norm_2(m_point - cam_ctr);
  vw::Vector3 ground_dir = (m_point - cam_ctr) / dist_to_ground;
  
  // Normalized direction from pixel to ground
  vw::Vector3 pix_dir = m_model->pixel_to_vector(pix);
  
  vw::Vector3 diff = pix_dir - ground_dir;

  // Find the components on this in a plane that is mostly perpendicular to vectors
  // from the cameras to the ground, so the diff vector is projected onto this plane.
  double dot1 = dot_prod(diff, m_perp1);
  double dot2 = dot_prod(diff, m_perp2);
  
  // Multiply by dist to ground so we can measure these at ground level. This avoids
  // numerical issues.
  dot1 *= dist_to_ground;
  dot2 *= dist_to_ground;

  return vw::Vector2(dot1, dot2);
}

}} // namespace asp::camera
