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

#include <vw/Camera/CameraTransform.h>

namespace vw {

SmartPtrCameraTransform::SmartPtrCameraTransform(boost::shared_ptr<CameraModel> const& src_camera,
                                                 boost::shared_ptr<CameraModel> const& dst_camera) :
  m_src_camera(src_camera), m_dst_camera(dst_camera) {}

Vector2 SmartPtrCameraTransform::reverse(const Vector2 &p) const {
  VW_ASSERT(m_src_camera->camera_center(p) == m_dst_camera->camera_center(p),
            LogicErr() << "CameraTransformFunctor: Camera transformation requires that "
                       << "the camera center is always the same for both cameras.");

  // (1) Call src PixelToVector to find the vector emanating from the camera center.
  Vector3 vec = m_dst_camera->pixel_to_vector(p);

  // (2) Take resulting vector and call dest camera's VectorToPixel on it
  return m_src_camera->point_to_pixel(vec + m_dst_camera->camera_center(p));
}

Vector2 SmartPtrCameraTransform::forward(const Vector2 &p) const {
  VW_ASSERT(m_src_camera->camera_center(p) == m_dst_camera->camera_center(p),
            LogicErr() << "CameraTransformFunctor: Camera transformation requires that "
                       << "the camera center is always the same for both cameras.");

  // (1) Call src PixelToVector to find the vector emanating from the camera center.
  Vector3 vec = m_src_camera->pixel_to_vector(p);

  // (2) Take resulting vector and call dest camera's VectorToPixel on it.
  return m_dst_camera->point_to_pixel(vec + m_src_camera->camera_center(p));
}

} // namespace vw
