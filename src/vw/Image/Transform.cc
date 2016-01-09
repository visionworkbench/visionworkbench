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


/// Transform.cc
///
///
// ImageView classes for transforming the domain of an image
// (image warping).

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/Image/Transform.h>

namespace vw{
  // Helper functions to pull and push a matrix to an affine transform

  Matrix3x3 affine2mat(AffineTransform const& transform){
    Vector2 O = transform.forward(Vector2(0, 0));
    Vector2 A = transform.forward(Vector2(1, 0));
    Vector2 B = transform.forward(Vector2(0, 1));

    A -= O;
    B -= O;
    Matrix3x3 T;
    T.set_identity();
    T(0, 0) = A[0]; T(1, 0) = A[1];
    T(0, 1) = B[0]; T(1, 1) = B[1];
    T(0, 2) = O[0]; T(1, 2) = O[1];
    return T;
  }

  AffineTransform mat2affine(Matrix3x3 const& T){
    return AffineTransform(submatrix(T, 0, 0, 2, 2), Vector2(T(0, 2), T(1, 2)));
  }
  
std::ostream& operator<<( std::ostream& os, AffineTransform const& trans ) {
  std::ostringstream oss; // To use custom precision
  oss.precision(10);
  oss << "AffineTransform: " << trans.m_a  << ", " << trans.m_b  << ", " 
                             << trans.m_c  << ", " << trans.m_d  << std::endl;
  oss << "               : " << trans.m_ai << ", " << trans.m_bi << ", " 
                             << trans.m_ci << ", " << trans.m_di << std::endl;
  oss << "               : " << trans.m_x  << ", " << trans.m_y  << std::endl;
  os << oss.str();
  return os;
}
  
  
std::ostream& operator<<( std::ostream& os, HomographyTransform const& trans ) {
  std::ostringstream oss; // To use custom precision
  oss.precision(10);
  oss << "HomograpyTransform: " << trans.m_H << std::endl;
  os << oss.str();
  return os;
}  
  
  
  
  
  
  
  
  
  
  
}
