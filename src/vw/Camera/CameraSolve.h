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


/// \file CameraSolve.h
///
#ifndef __CAMERA_CAMERA_SOLVE_H__
#define __CAMERA_CAMERA_SOLVE_H__

#include <vw/Camera/CameraModel.h>
#include <vw/Math/LevenbergMarquardt.h>

namespace vw {

namespace camera {

  // Levenberg Marquardt solver for line number (y) and column number (x) for
  // the given point in space (ie, point_to_pixel) The obtained solution pixel
  // (x, y) must be such that the vector from this camera pixel goes through the
  // given point. This version does uses only functions declared in CameraModel
  // class and should work for any camera. Some cameras may be able to use a
  // faster implementation that solves for y first and then for x separately.
  
  // See LinescanErr.h for an alternative 2D to 2D implementation which is about
  // 3x faster.
  
  class CameraGenericLMA: public math::LeastSquaresModelBaseFixed<CameraGenericLMA, 2, 3> {
    const CameraModel* m_model;
    Vector3 m_point;
  public:
    typedef Vector2 domain_type;     // 2D pixel, input to cost function vector
    typedef Vector3 result_type;     // 3D error, output of cost function vector
    typedef Matrix<double, 3, 2> jacobian_type;

    CameraGenericLMA( const CameraModel* model, const vw::Vector3& pt ) :
      m_model(model), m_point(pt) {}

    // Minimize the difference between the vector from the pixel and
    //  the vector from the point to the camera.
    inline result_type operator()( domain_type const& pix ) const {
      try {
        return m_model->pixel_to_vector(pix) - normalize(m_point - m_model->camera_center(pix));
      } catch(...) { // Handle camera misses by returning a very large error
        return Vector3(999999,999999,999999);
      }
    }
  }; // End class CameraGenericLMA
  
}} // end vw::camera

#endif//__CAMERA_CAMERA_SOLVE_H__
