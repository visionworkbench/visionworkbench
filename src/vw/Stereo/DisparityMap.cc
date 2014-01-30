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


#include <vw/Stereo/DisparityMap.h>

namespace vw {
namespace stereo {

  StdDevImageFunc::StdDevImageFunc(int32 kernel_width, int32 kernel_height) :
    m_kernel_width(kernel_width), m_kernel_height(kernel_height) {
    VW_ASSERT(m_kernel_width > 0 && m_kernel_height > 0,
              ArgumentErr() << "StdDevImageFunc: kernel sizes must be non-zero.");
  }


  BBox2i StdDevImageFunc::work_area() const {
    return BBox2i(Vector2i(-m_kernel_width/2, -m_kernel_height/2),
                  Vector2i(m_kernel_width, m_kernel_height));
  }


  // Compute the plane that best fits a set of 3D points.
  // - The plane is described as z = ax + by + c
  //( the output vector contains [a, b, c]
  bool fitPlaneToPoints(const std::vector<Vector3> &points, Vector3 &planeDesc){
    const size_t X = 0; // Convenience constants
    const size_t Y = 1;
    const size_t Z = 2;
    const size_t numPoints = points.size();
    
    // Compute values in a matrix A and vector b
    //A: [xx  xy  x]     B: [xz]
    //   [xy  yy  y]        [yz]
    //   [ x   y  n]        [ z]
    Matrix3x3 matA(0, 0, 0, 0, 0, 0, 0, 0, 0); // A symmetric matrix, A' = A
    Vector3   vecB(0, 0, 0);
    for (size_t i=0; i<numPoints; ++i)
      {
        matA[0][0] += points[i][X] * points[i][X]; // sum xx
        matA[0][1] += points[i][X] * points[i][Y]; // sum xy
        matA[0][2] += points[i][X];                // sum x
        matA[1][0] += points[i][X] * points[i][Y]; // sum xy
        matA[1][1] += points[i][Y] * points[i][Y]; // sum yy
        matA[1][2] += points[i][Y];                // sum y
        matA[2][0] += points[i][X];                // sum x
        matA[2][1] += points[i][Y];                // sum y
      
        vecB[0]    += points[i][X] * points[i][Z]; // sum xz
        vecB[1]    += points[i][Y] * points[i][Z]; // sum yz
        vecB[2]    += points[i][Z];                // sum z
      }
    matA[2][2] = numPoints; // n
        
    // Now solve Ax = b (3x3)*(3x1) = (3x1)
    planeDesc = vw::math::solve(matA, vecB); // Throws!

    return true;
  }

  /// Computes the distance from a point to a plane
  double pointToPlaneDist(const Vector3 &point, const Vector3 &planeDesc){
    // Convert plane from format z = ax + by + c to format
    // 0 = ax + by + cz + d = 0
    double a = planeDesc[0];
    double b = planeDesc[1];
    double c = -1.0;
    double d = planeDesc[2];
   
    double numerator   = fabs(a*point[0] + b*point[1] + c*point[2] + d);
    double denominator = sqrt(a*a + b*b + c*c);
   
    //TODO: This will crash if an invalid plane is given! 
    return numerator / denominator;
  }

  /// Compute statistics for how well points fit a plane 
  bool checkPointToPlaneFit(const std::vector<Vector3> &points,
                            const Vector3 &planeDesc,
                            double &meanError, double &stdDevError){
    const size_t numPoints = points.size();
    std::vector<double> dists(numPoints);
    
    // Compute the mean
    meanError = 0;
    for (size_t i=0; i<numPoints; ++i)
      {
        dists[i] = pointToPlaneDist(points[i], planeDesc);
        meanError += dists[i];
      }
    meanError /= static_cast<double>(numPoints);
    
    // Compute the standard deviation
    double sumDiff = 0;
    for (size_t i=0; i<numPoints; ++i)
      {
        //double diffSq = (dists[i] - meanError) * (dists[i] - meanError);
        double diffSq = dists[i]*dists[i];
        sumDiff += diffSq;
      }
    stdDevError = sqrt(sumDiff / static_cast<double>(numPoints));
    
    return true;
  }
  
}}    // namespace vw::stereo
