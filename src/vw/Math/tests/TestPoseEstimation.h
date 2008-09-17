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

// TestPoseEstimation.h
#include <cxxtest/TestSuite.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/PoseEstimation.h>

using namespace vw;
using namespace vw::math;

class PoseEstimation : public CxxTest::TestSuite
{
public:
  void test_relative_orientation()
  {
    Matrix3x3 vectors1;
    select_col(vectors1,0) = normalize( Vector3(1,2,3) );
    select_col(vectors1,1) = normalize( Vector3(1,-2,3) );
    select_col(vectors1,2) = normalize( Vector3(1,2,0) );
    Quat q1 = normalize( Quat(1,2,3,4) );
    Matrix3x3 vectors2 = q1.rotation_matrix() * vectors1;
    Quat q2 = relative_orientation( vectors1, vectors2 );
    if( norm_2(q2+q1) < norm_2(q2-q1) ) q2 = -q2;
    TS_ASSERT_DELTA( q1.w(), q2.w(), 1e-6 );
    TS_ASSERT_DELTA( q1.x(), q2.x(), 1e-6 );
    TS_ASSERT_DELTA( q1.y(), q2.y(), 1e-6 );
    TS_ASSERT_DELTA( q1.z(), q2.z(), 1e-6 );
  }
};
