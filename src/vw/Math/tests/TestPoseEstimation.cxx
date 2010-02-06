// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Quaternion.h>
#include <vw/Math/PoseEstimation.h>

using namespace vw;
using namespace vw::math;

TEST(PoseEstimation, RelativeOrientation)
{
  Matrix3x3 vectors1;
  select_col(vectors1,0) = normalize( Vector3(1,2,3) );
  select_col(vectors1,1) = normalize( Vector3(1,-2,3) );
  select_col(vectors1,2) = normalize( Vector3(1,2,0) );
  Quat q1 = normalize( Quat(1,2,3,4) );
  Matrix3x3 vectors2 = q1.rotation_matrix() * vectors1;
  Quat q2 = relative_orientation( vectors1, vectors2 );
  if( norm_2(q2+q1) < norm_2(q2-q1) ) q2 = -q2;

  EXPECT_NEAR( q2.w(), q1.w(), 1e-15 );
  EXPECT_NEAR( q2.x(), q1.x(), 1e-15 );
  EXPECT_NEAR( q2.y(), q1.y(), 1e-15 );
  EXPECT_NEAR( q2.z(), q1.z(), 1e-15 );
}
