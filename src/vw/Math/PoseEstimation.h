// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

//! \file PoseEstimation.h
//!
//! Routines for estimating relative poses in various geometrical
//! settings.  At the moment the only thing is relative_orientation(),
//! which computes the relative orientation of two sets of unit
//! 3-vectors.
//!
#ifndef __VW_MATH_POSEESTIMATION_H__
#define __VW_MATH_POSEESTIMATION_H__

#include <vw/Core/Exception.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Quaternion.h>
#include <vw/Math/LinearAlgebra.h>

namespace vw {
namespace math {

  // Compute the relative orientation of two sets of unit 3-vectors,
  // represented as 3xn matrices.  Specifically, computes the unit
  // quaternion q that minimizes sum(|v2-q.rotate(v1)|^2).  This can
  // be reduced to minimizing a quadratic form over unit 4-vectors,
  // which further reduces to an eigenvector problem.
  template <class Matrix1T, class Matrix2T>
  Quaternion<float64> relative_orientation( MatrixBase<Matrix1T> const& vectors1,
                                            MatrixBase<Matrix2T> const& vectors2 )
  {
    VW_ASSERT( vectors1.impl().rows()==3 && vectors2.impl().rows()==3 && 
               vectors1.impl().cols()==vectors2.impl().cols(),
               ArgumentErr( "relative_orientation(): Invalid or incompatible matrix dimensions" ) );
    Matrix4x4 M;
    for( int i=0; i<vectors1.impl().cols(); ++i ) {
      MatrixCol<const Matrix1T> v1(vectors1.impl(),i);
      MatrixCol<const Matrix2T> v2(vectors2.impl(),i);
      M(0,0) += -dot_prod(v1,v2);
      M(0,1) += v1(2)*v2(1) - v1(1)*v2(2);
      M(0,2) += v1(0)*v2(2) - v1(2)*v2(0);
      M(0,3) += v1(1)*v2(0) - v1(0)*v2(1);
      M(1,1) += - v1(0)*v2(0) + v1(1)*v2(1) + v1(2)*v2(2);
      M(2,2) +=   v1(0)*v2(0) - v1(1)*v2(1) + v1(2)*v2(2);
      M(3,3) +=   v1(0)*v2(0) + v1(1)*v2(1) - v1(2)*v2(2);
      M(1,2) -= v1(0)*v2(1) + v2(0)*v1(1);
      M(1,3) -= v1(2)*v2(0) + v2(2)*v1(0);
      M(2,3) -= v1(1)*v2(2) + v2(1)*v1(2);
    }
    M(1,0) = M(0,1);
    M(2,0) = M(0,2);
    M(2,1) = M(1,2);
    M(3,0) = M(0,3);
    M(3,1) = M(1,3);
    M(3,2) = M(2,3);
    Vector<std::complex<float64>,4> evals;
    Matrix<std::complex<float64>,4,4> evecs;
    eigen( M, evals, evecs );
    Quat result( real( select_col( evecs, 0 ) ) );
    return result;
  }

} // namespace math
} // namespace vw

#endif // __VW_MATH_POSEESTIMATION_H__
