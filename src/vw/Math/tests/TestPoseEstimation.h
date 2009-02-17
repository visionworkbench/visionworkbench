// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
