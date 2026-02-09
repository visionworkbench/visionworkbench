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


#include <vw/Camera/CameraGeometry.h>
#include <vw/InterestPoint/InterestPoint.h>

using namespace vw;
using namespace camera;

template <>
Vector3 camera::Convert2HomogenousVec2( ip::InterestPoint const& ip ) {
  return Vector3( ip.x, ip.y , 1 );
}

Matrix<double>
CameraMatrixFittingFunctor::BasicDLT( std::vector<Vector<double> > const& input,
                                      std::vector<Vector<double> > const& output ) const {
  VW_ASSERT( input[0].size() == 4 && output[0].size() == 3,
             vw::ArgumentErr() << "Camera Matrix requires Vector4 inputs and Vector3 outputs." );

  vw::Matrix<double,12,12> A;
  for ( unsigned i = 0; i < 6; i++ ) { // Measure iterator
    for ( unsigned j = 0; j < 4; j++ ) { // Measure's X elem iterator
      // Filling w*Xt
      A(2*i+1,j) = output[i][2]*input[i][j];
      // Filling -w*Xt
      A(2*i,j+4) = -output[i][2]*input[i][j];
      // Filling y*Xt
      A(2*i,j+8) = output[i][1]*input[i][j];
      // Filling -x*Xt
      A(2*i+1,j+8) = -output[i][0]*input[i][j];
    }
  }

  Matrix<double,3,4> p;
  // SVD for smallest singular value
  Matrix<double> U,VT;
  Vector<double> S;
  svd(A,U,S,VT);
  submatrix(p,0,0,1,4) = submatrix(VT,VT.rows()-1,0,1,4);
  submatrix(p,1,0,1,4) = submatrix(VT,VT.rows()-1,4,1,4);
  submatrix(p,2,0,1,4) = submatrix(VT,VT.rows()-1,8,1,4);
  return p;
}

CameraMatrixFittingFunctor::CameraMatrixModelLMA::CameraMatrixModelLMA
( std::vector<Vector<double> > const& input,
  std::vector<Vector<double> > const& output ) :
  m_world_input(input), m_image_output(output) {}

Vector<double>
CameraMatrixFittingFunctor::CameraMatrixModelLMA::operator()( Vector<double> const& x ) const {
  result_type output;
  output.set_size( m_world_input.size() );
  Matrix<double,3,4> P = unflatten(x);

  for ( uint32 i = 0; i < m_world_input.size(); i++ ) {
    Vector3 reproj = P*m_world_input[i];
    reproj /= reproj[2];
    output[i] = norm_2( subvector(m_image_output[i],0,2) - subvector(reproj,0,2) );
  }
  return output;
}


Vector<double>
CameraMatrixFittingFunctor::CameraMatrixModelLMA::flatten( Matrix<double> const& input ) const {
  Vector<double,12> output;
  for ( uint8 i = 0; i < 3; i++ ) {
    for ( uint8 j = 0; j < 4; j++ ) {
      output(4*i+j) = input(i,j);
    }
  }
  return output;
}

Matrix<double>
CameraMatrixFittingFunctor::CameraMatrixModelLMA::unflatten( Vector<double> const& input ) const {
  Matrix<double,3,4> output;
  for ( uint8 i = 0; i < 3; i++ ) {
    for ( uint8 j = 0; j < 4; j++ ) {
      output(i,j) = input(4*i+j);
    }
  }
  return output;
}

Vector<double>
FundamentalMatrixMLFittingFunctor::ProjectiveModelLMA::operator()( Vector<double> const& x )  const {
  unsigned number_of_measures = (x.size()-12)/3;
  Vector<double> output( number_of_measures*4 );

  // Create second camera
  Matrix<double> P2(3,4);
  select_row(P2,0) = subvector(x,0,4);
  select_row(P2,1) = subvector(x,4,4);
  select_row(P2,2) = subvector(x,8,4);

  // Calculate reprojection error
  for ( unsigned i = 0, j = 12; i < number_of_measures; i++, j+= 3 ) {
    Vector2 reprojection1(x[j]/x[j+2],x[j+1]/x[j+2]);
    subvector(output,4*i,2) = reprojection1;
    Vector3 reprojection2 = P2*Vector4(x[j],x[j+1],x[j+2],1);
    reprojection2 /= reprojection2[2];
    subvector(output,4*i+2,2) = subvector(reprojection2,0,2);
  }

  return output;
}
