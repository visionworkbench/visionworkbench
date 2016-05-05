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

#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/Geometry.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>

#include <math.h>

namespace vw{
namespace math{

double normalize_longitude(double lon, bool center_on_zero) {
  const double MULTIPLE = 360.0;
  double min, max;
  if (center_on_zero) {
    min=-180.0;
    max= 180.0;
  } else {
    min=  0.0;
    max=360.0;
  }
 
  if (lon < min){
    double factor = (min - lon) / MULTIPLE;
    lon += MULTIPLE*ceil(factor);
  }
  if (lon > max){
    double factor = (lon - max) / MULTIPLE;
    lon -= MULTIPLE*ceil(factor);
  }
   
  return lon;
}


double degree_diff(double d1, double d2) {
  double diff  = fabs(d1 - d2);
  double diff2 = fabs(d1 - (d2+360));
  double diff3 = fabs(d1 - (d2-360));
  if (diff2 < diff)
    diff = diff2;
  if (diff3 < diff)
    diff = diff3;
  return diff;
}


}} // End namespace vw::math

using namespace vw;
using namespace vw::math;


/// Solve for Normalization Similarity Matrix used for noise rej.
Matrix3x3
HomographyFittingFunctor::NormSimilarity( std::vector<Vector3> const& pts ) const {
  size_t num_points = pts.size();
  size_t dimension = 3;

  Vector2 translation;
  for ( size_t i = 0; i < num_points; i++ )
    translation+=subvector(pts[i],0,dimension-1);
  translation /= num_points;

  double scale = 0;
  for ( size_t i = 0; i < num_points; i++ )
    scale += norm_2( subvector(pts[i],0,dimension-1) - translation );
  scale = num_points*sqrt(2.)/scale;

  Matrix3x3 t;
  t(2,2) = 1;
  t(0,0) = scale;
  t(1,1) = scale;
  t(0,2) = -scale*translation[0];
  t(1,2) = -scale*translation[1];
  return t;
}

vw::Matrix3x3
HomographyFittingFunctor::BasicDLT( std::vector<Vector3 > const& input,
                                    std::vector<Vector3 > const& output )  const {
  VW_ASSERT( input.size() == 4 && output.size() == 4,
             vw::ArgumentErr() << "DLT in this implementation expects to have only 4 inputs." );
  VW_ASSERT( input[0][input[0].size()-1] == 1,
             vw::ArgumentErr() << "Input data doesn't seem to be normalized.");
  VW_ASSERT( output[0][output[0].size()-1] == 1,
             vw::ArgumentErr() << "Secondary input data doesn't seem to be normalized.");
  VW_ASSERT( input[0].size() == 3,
             vw::ArgumentErr() << "BasicDLT only supports homogeneous 2D vectors.");

  vw::Matrix<double,8,9> A;
  for ( uint8 i = 0; i < 4; i++ )
    for ( uint8 j = 0; j < 3; j++ ) {
      // Filling in -wi'*xi^T
      A(i,j+3) = -output[i][2]*input[i][j];
      // Filling in yi'*xi^T
      A(i,j+6) = output[i][1]*input[i][j];
      // Filling in wi'*xi^T
      A(i+4,j) = output[i][2]*input[i][j];
      // Filling in -xi'*xi^T
      A(i+4,j+6) = -output[i][0]*input[i][j];
    }

  Matrix<double> nullsp = nullspace(A);
  nullsp /= nullsp(8,0);
  Matrix3x3 H;
  for ( uint8 i = 0; i < 3; i++ )
    for ( uint8 j = 0; j < 3; j++ )
      H(i,j) = nullsp(i*3+j,0);
  return H;
}

