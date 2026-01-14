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
using namespace std;

namespace vw { namespace stereo {

  // Two-ray triangulation. Triangulate the point by finding the
  // midpoint of the segment joining the closest points on the two
  // rays emanating from the camera.
  vw::Vector3 triangulate_pair(vw::Vector3 const& dir0, vw::Vector3 const& ctr0, 
                               vw::Vector3 const& dir1, vw::Vector3 const& ctr1, 
                               vw::Vector3& errorVec) {
    
    vw::Vector3 v12 = cross_prod(dir0, dir1);
    vw::Vector3 v1 = cross_prod(v12, dir0);
    vw::Vector3 v2 = cross_prod(v12, dir1);
    
    vw::Vector3 closestPoint1 = ctr0 + dot_prod(v2, ctr1-ctr0)/dot_prod(v2, dir0)*dir0;
    vw::Vector3 closestPoint2 = ctr1 + dot_prod(v1, ctr0-ctr1)/dot_prod(v1, dir1)*dir1;

    errorVec = closestPoint1 - closestPoint2;
    return 0.5 * (closestPoint1 + closestPoint2);
  }


  
// Constructor with n cameras
StereoModel::StereoModel(vector<const camera::CameraModel *> const& cameras,
                         double angle_tol):
  m_cameras(cameras),
  m_angle_tol(angle_tol) {}
  
// Constructor with two cameras
StereoModel::StereoModel(camera::CameraModel const* camera_model1,
                         camera::CameraModel const* camera_model2,
                         double angle_tol) {
  m_cameras.clear();
  m_cameras.push_back(camera_model1);
  m_cameras.push_back(camera_model2);
  m_angle_tol = angle_tol;
}

bool StereoModel::are_nearly_parallel(double angle_tol,
                                      std::vector<Vector3> const& camDirs) {

  // If the camera directions are nearly parallel, there will be very
  // large numerical uncertainty about where to place the point.  We
  // set a threshold here to reject points that are on nearly parallel
  // rays.  The threshold of 1e-4 corresponds to a convergence of less
  // than theta = 0.81 degrees, so if the two rays are within 0.81
  // degrees of being parallel, we reject this point.
  //
  // This threshold was chosen empirically for now, but should
  // probably be revisited once a more rigorous analysis has
  // been completed. -mbroxton (11-MAR-07)
  double tol = 1e-4;

  if (angle_tol > 0) tol = angle_tol; // can be over-ridden from the outside

  bool are_par = true;
  for (int p = 0; p < int(camDirs.size()) - 1; p++){
    if ( 1 - dot_prod(camDirs[p], camDirs[p+1]) >= tol )
      are_par = false;
  }
  return are_par;
}

// Compute the rays intersection.

// Note: Classes RPCStereoModel and BathyStereoModel inherit from this
// class and re-implement this function.
Vector3 StereoModel::operator()(vector<Vector2> const& pixVec,
                                Vector3& errorVec) const {
  
  // Initialize the outputs
  errorVec = Vector3();
  
  int num_cams = m_cameras.size();
  VW_ASSERT((int)pixVec.size() == num_cams,
            vw::ArgumentErr() << "the number of rays must match "
            << "the number of cameras.\n");
  
  try {

    vector<Vector3> camDirs(num_cams), camCtrs(num_cams);
    camDirs.clear(); camCtrs.clear();
    
    // Pick the valid rays
    for (int p = 0; p < num_cams; p++){
      
      Vector2 pix = pixVec[p];
      if (pix != pix || // i.e., NaN
          pix == camera::CameraModel::invalid_pixel() ) 
        continue;
      
      camDirs.push_back(m_cameras[p]->pixel_to_vector(pix));
      camCtrs.push_back(m_cameras[p]->camera_center(pix));
    }

    // Not enough valid rays
    if (camDirs.size() < 2) 
      return Vector3();

    if (are_nearly_parallel(m_angle_tol, camDirs)) 
      return Vector3();

    // Determine range by triangulation
    Vector3 result = triangulate_point(camDirs, camCtrs, errorVec);
    
    // Reflect points that fall behind one of the two cameras.
    bool reflect = false;
    for (int p = 0; p < (int)camCtrs.size(); p++)
      if (dot_prod(result - camCtrs[p], camDirs[p]) < 0 ) reflect = true;
    if (reflect)
      result = -result + 2*camCtrs[0];

    return result;
    
  } catch (const camera::PixelToRayErr& /*e*/) {
  }
  return Vector3();
}

Vector3 StereoModel::operator()(vector<Vector2> const& pixVec,
                                double& error) const {
  Vector3 errorVec;
  Vector3 result = operator()(pixVec, errorVec);
  error = norm_2(errorVec);
  return result;
}
  
Vector3 StereoModel::operator()(Vector2 const& pix1,
                                Vector2 const& pix2, Vector3& errorVec) const {
  vector<Vector2> pixVec;
  pixVec.push_back(pix1);
  pixVec.push_back(pix2);
  return operator()(pixVec, errorVec); 
}


Vector3 StereoModel::operator()(Vector2 const& pix1, Vector2 const& pix2,
                                double& error) const {
  Vector3 errorVec;
  Vector3 result = operator()(pix1, pix2, errorVec);
  error = norm_2(errorVec);
  return result;
}

double StereoModel::convergence_angle(Vector2 const& pix1, Vector2 const& pix2) const {
  return acos(dot_prod(m_cameras[0]->pixel_to_vector(pix1),
                       m_cameras[1]->pixel_to_vector(pix2)));
}

// The default 1-cos(x) function does badly when x is close to 0,
// which then leads to incorrect angle tolerance.
double StereoModel::robust_1_minus_cos(double x){
  if (std::abs(x) > 1e-7) 
    return 1-cos(x);
  return x*x/2;
}
  
Vector3 StereoModel::triangulate_point(vector<Vector3> const& camDirs,
                                       vector<Vector3> const& camCtrs,
                                       Vector3& errorVec){
  

  int num_cams = camDirs.size();
  if ( num_cams == 2 ){
    // Two-ray triangulation
    return triangulate_pair(camDirs[0], camCtrs[0], camDirs[1], camCtrs[1], errorVec);
  }

  // Multi-ray triangulation. Find the intersection of the rays in
  // least squares sense (the point from which the sum of square
  // distances to the rays is smallest).

  // For two rays, this will give the same result as the code above,
  // but it is a bit slower.
    
  // Based on:
  // Optimal Ray Intersection For Computing 3D Points
  // From N-View Correspondences
  // Greg Slabaugh, Ron Schafer, Mark Livingston
  // http://www.soi.city.ac.uk/~sbbh653/publications/opray.pdf

  Matrix<double,3,3> M(0.0);
  Vector3 R;
  for (int t = 0; t < num_cams; t++){
    Vector3 D = camDirs[t], C = camCtrs[t];
    double a = D[0], b = D[1], c = D[2], x = C[0], y = C[1], z = C[2];
    M(0, 0) = M(0, 0) + 1-a*a;
    M(1, 1) = M(1, 1) + 1-b*b;
    M(2, 2) = M(2, 2) + 1-c*c;
    M(0, 1) = M(0, 1) - a*b;
    M(0, 2) = M(0, 2) - a*c;
    M(1, 2) = M(1, 2) - b*c;
    R(0) = R(0) + (1-a*a)*x - a*b*y - a*c*z;
    R(1) = R(1) + (1-b*b)*y - a*b*x - b*c*z;
    R(2) = R(2) + (1-c*c)*z - a*c*x - b*c*y;
  }
  M(1, 0) = M(0, 1);
  M(2, 0) = M(0, 2);
  M(2, 1) = M(1, 2);
  Vector3 P = inverse(M)*R;
    
  // Find 2 times average distance from the intersection point to the
  // rays. For two rays, this will agree with the shortest distance
  // between rays.
  double err = 0.0;
  double px = P[0], py = P[1], pz = P[2];
  for (int t = 0; t < num_cams; t++){
    Vector3 D = camDirs[t], C = camCtrs[t];
    double a = D[0], b = D[1], c = D[2], x = C[0], y = C[1], z = C[2];
    double v = a*(px-x) + b*(py-y) + c*(pz-z);
    double dist = (px-x)*(px-x) + (py-y)*(py-y) + (pz-z)*(pz-z) - v*v;
    if (dist < 0) dist = 0; // if by some numerical fluke dist is negative
    dist = std::sqrt(dist);
    err += dist;
  }
  err = 2.0*err/num_cams;
    
  // For more than two cameras, the error vector between rays is not
  // meaningful.
  errorVec = Vector3(err, 0, 0);
    
  return P;
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
      printf("\tStereoModel computing points: %0.2f%% complete.\r",
             100.0f * float(y)/disparity_map.rows());
      fflush(stdout);
    }
    for (int32 x = 0; x < disparity_map.cols(); x++) {
      if ( is_valid(disparity_map(x,y)) ) {
        xyz(x,y) = (*this)(Vector2( x, y),
                           Vector2( x+disparity_map(x,y)[0],
                                    y+disparity_map(x,y)[1]),
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
           << ",  Max error = " << max_error << endl;
    return xyz;
}
  
}} // vw::stereo
