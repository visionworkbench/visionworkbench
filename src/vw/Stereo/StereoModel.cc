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


#include <vw/Image/ImageView.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Math/Vector.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Core/Exception.h>

namespace vw {
namespace stereo {
namespace detail {
  class PointLMA : public math::LeastSquaresModelBase<PointLMA> {
    std::vector<const camera::CameraModel *> m_cameras;
    
  public:
    typedef Vector<double,4> result_type;
    typedef Vector<double,3> domain_type;
    typedef Matrix<double> jacobian_type;
    
    PointLMA(std::vector<const camera::CameraModel *>  const& cameras):
      m_cameras(cameras) {}
    
    inline result_type operator()( domain_type const& x ) const {
      Vector4 output;
      subvector(output,0,2) = m_cameras[0]->point_to_pixel( x );
      subvector(output,2,2) = m_cameras[1]->point_to_pixel( x );
      return output;
    }
  };
}
  
// Constructor with n cameras
StereoModel::StereoModel(std::vector<const camera::CameraModel *> const& cameras,
                         bool least_squares_refine ) :
  m_cameras(cameras),
  m_least_squares(least_squares_refine) {}
  
// Constructor with two cameras
StereoModel::StereoModel(camera::CameraModel const* camera_model1,
                         camera::CameraModel const* camera_model2,
                         bool least_squares_refine ){
  m_cameras.clear();
  m_cameras.push_back(camera_model1);
  m_cameras.push_back(camera_model2);
  m_least_squares = least_squares_refine;
}
  
ImageView<Vector3>
StereoModel::operator()(ImageView<PixelMask<Vector2f> > const& disparity_map,
                        ImageView<double> &error) const {
  
  // Error analysis
  double mean_error = 0.0;
  double max_error = 0.0;
  int32 point_count = 0;
  int32 divergent = 0;

  // Allocate xyz image and get pointer to buffer
  ImageView<Vector3> xyz(disparity_map.cols(), disparity_map.rows());
  error.set_size(disparity_map.cols(), disparity_map.rows());

  // Compute 3D position for each pixel in the disparity map
  vw_out() << "StereoModel: Applying camera models\n";
  for (int32 y = 0; y < disparity_map.rows(); y++) {
    if (y % 100 == 0) {
      printf("\tStereoModel computing points: %0.2f%% complete.\r", 100.0f*float(y)/disparity_map.rows());
      fflush(stdout);
    }
    for (int32 x = 0; x < disparity_map.cols(); x++) {
      if ( is_valid(disparity_map(x,y)) ) {
        xyz(x,y) = (*this)(Vector2( x, y),
                           Vector2( x+disparity_map(x,y)[0], y+disparity_map(x,y)[1]),
                           error(x,y) );

        if (error(x,y) >= 0) {
          // Keep track of error statistics
          if (error(x,y) > max_error)
            max_error = error(x,y);
          mean_error += error(x,y);
            ++point_count;
        } else {
          // rays diverge or are parallel
          xyz(x,y) = Vector3();
          divergent++;
        }
      } else {
        xyz(x,y) = Vector3();
        error(x,y) = 0;
      }
      }
  }

  if (divergent != 0)
    vw_out() << "WARNING in StereoModel: " << divergent
             << " rays diverged or were parallel!\n";

  vw_out() << "\tStereoModel computing points: Done.                  \n";
  vw_out() << "\tMean error = " << mean_error/double(point_count)
           << ",  Max error = " << max_error << std::endl;
    return xyz;
}

bool StereoModel::are_nearly_parallel(Vector3 const& vec1,
                                      Vector3 const& vec2) const{
  // If vec1 and vec2 are nearly parallel, there will be
  // very large numerical uncertainty about where to place the
  // point.  We set a threshold here to reject points that are
  // on nearly parallel rays.  The threshold of 1e-4 corresponds
  // to a convergence of less than theta = 0.81 degrees, so if
  // the two rays are within 0.81 degrees of being parallel, we
  // reject this point.
  //
  // This threshold was chosen empirically for now, but should
  // probably be revisited once a more rigorous analysis has
  // been completed. -mbroxton (11-MAR-07)
  return ( (1-dot_prod(vec1, vec2) < 1e-4 && !m_least_squares) ||
           (1-dot_prod(vec1, vec2) < 1e-5 && m_least_squares) );
}

Vector3 StereoModel::operator()(std::vector<Vector2> const& pixVec,
                                Vector3& errorVec) const {
  
  // Note: Class RPCStereoModel inherits from this class and
  // re-implements this function.

  VW_ASSERT(pixVec.size() == m_cameras.size(),
            vw::ArgumentErr() << "the number of rays must match "
            << "the number of cameras.\n");
            
  errorVec = Vector3();

  // Check for NaN and invalid pixels
  if (pixVec[0] != pixVec[0] || pixVec[1] != pixVec[1]) return Vector3();
  if (pixVec[0] == camera::CameraModel::invalid_pixel() ||
      pixVec[1] == camera::CameraModel::invalid_pixel()
      ) return Vector3();
  
  try {
    // Determine range by triangulation
    Vector3 vec1 = m_cameras[0]->pixel_to_vector(pixVec[0]);
    Vector3 vec2 = m_cameras[1]->pixel_to_vector(pixVec[1]);
    Vector3 origin1 = m_cameras[0]->camera_center(pixVec[0]);
    Vector3 origin2 = m_cameras[1]->camera_center(pixVec[1]);

    if (are_nearly_parallel(vec1, vec2)){
      return Vector3();
    }
    Vector3 result = triangulate_point(origin1, vec1,
                                       origin2, vec2,
                                       errorVec);
    
    if ( m_least_squares ){
      if (pixVec.size() == 2)
        refine_point(pixVec[0], pixVec[1], result);
      else
        vw::vw_throw(vw::NoImplErr() << "Least squares refinement is not "
                     << "implemented for multi-view stereo.");
    }
    
    // Reflect points that fall behind one of the two cameras
    if ( dot_prod(result - origin1, vec1) < 0 ||
         dot_prod(result - origin2, vec2) < 0 ) {
      result = -result + 2*origin1;
    }

    return result;

  } catch (const camera::PixelToRayErr& /*e*/) {
    return Vector3();
  }
}

Vector3 StereoModel::operator()(std::vector<Vector2> const& pixVec,
                                double& error) const {
  Vector3 errorVec;
  Vector3 result = operator()(pixVec, errorVec);
  error = norm_2(errorVec);
  return result;
}
  
Vector3 StereoModel::operator()(Vector2 const& pix1,
                                Vector2 const& pix2, Vector3& errorVec) const {
  std::vector<Vector2> pixVec;
  pixVec.push_back(pix1);
  pixVec.push_back(pix2);
  return operator()(pixVec, errorVec); 
}


Vector3 StereoModel::operator()(Vector2 const& pix1, Vector2 const& pix2,
                                double& error ) const {
  Vector3 errorVec;
  Vector3 result = operator()(pix1, pix2, errorVec);
  error = norm_2(errorVec);
  return result;
}

double StereoModel::convergence_angle(Vector2 const& pix1, Vector2 const& pix2) const {
  return acos(dot_prod(m_cameras[0]->pixel_to_vector(pix1),
                       m_cameras[1]->pixel_to_vector(pix2)));
}

Vector3 StereoModel::triangulate_point(Vector3 const& point1,
                                       Vector3 const& vec1,
                                       Vector3 const& point2,
                                       Vector3 const& vec2,
                                       Vector3& errorVec) const {

  // Triangulate the point by finding the midpoint of the segment
  // joining the closest points on the two rays emanating
  // from the camera.

  Vector3 v12 = cross_prod(vec1, vec2);
  Vector3 v1 = cross_prod(v12, vec1);
  Vector3 v2 = cross_prod(v12, vec2);

  Vector3 closestPoint1 = point1 + dot_prod(v2, point2-point1)/dot_prod(v2, vec1)*vec1;
  Vector3 closestPoint2 = point2 + dot_prod(v1, point1-point2)/dot_prod(v1, vec2)*vec2;

  errorVec = closestPoint1 - closestPoint2;

  return 0.5 * (closestPoint1 + closestPoint2);
}

void StereoModel::refine_point(Vector2 const& pix1,
                               Vector2 const& pix2,
                               Vector3& point) const {

  // Refine the point by minimizing the least squares error in pixel domain.

  detail::PointLMA model( m_cameras );
  Vector4 objective( pix1[0], pix1[1], pix2[0], pix2[1] );
  int status = 0;
  Vector3 npoint = levenberg_marquardt( model, point,
                                        objective, status, 1e-3, 1e-6, 10 );
  if ( status > 0 )
    point = npoint;
}

}} // vw::stereo
