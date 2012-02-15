// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <test/Helpers.h>
#include <vw/Camera/Extrinsics.h>

using namespace vw;
using namespace camera;

TEST( Extrinsics, HermitePositionInterpolation ) {
  std::vector<Vector3> m_p, m_v;

  {
    m_p.push_back( Vector3() ); m_p.push_back( Vector3(1,1,1) );
    m_v.push_back( Vector3(1,1,1) ); m_v.push_back( Vector3(1,1,1) );
    HermitePositionInterpolation interp( m_p, m_v, 0, 1 );
    EXPECT_VECTOR_NEAR( Vector3(1,1,1)/2, interp(0.5), 1e-6 );
  }

  {
    // Hermite expects a constant acceleration
    // a = c0
    // v = c0 * t + c1
    // p = c0 * c0 * t^2 / 2 + c1 * t + c2
    m_p.clear(); m_v.clear();
    m_p.push_back( Vector3(5,-3,1) );
    m_p.push_back( Vector3(7.02,-1.875,0.5) );
    m_v.push_back( Vector3(2,1,-1) );
    m_v.push_back( Vector3(1.8,1.5,0) );

    HermitePositionInterpolation interp( m_p, m_v, 0, 1.0 );
    EXPECT_VECTOR_NEAR( Vector3(6.005,-2.46875,.625),
                        interp( 0.5 ), 1e-1 );
    EXPECT_VECTOR_EQ( Vector3(5,-3,1), interp(0) );
    EXPECT_VECTOR_NEAR( Vector3(7.02,-1.875,0.5), interp(1.0 - 1e-7), 1e-6 );
    EXPECT_THROW( interp(-.5), ArgumentErr );
    EXPECT_THROW( interp(1.5), ArgumentErr );
  }
}

TEST( Extrinsics, LinearTime ) {
  LinearTimeInterpolation time( 5, 0.1 );

  EXPECT_NEAR( 5, time(0), 1e-12 );
  EXPECT_NEAR( 6, time(10), 1e-12 );
  EXPECT_NEAR( 4, time(-10), 1e-12 );
  EXPECT_NEAR( -1, time(-60), 1e-12 );
}

TEST( Extrinsics, TLCTime ) {
  std::vector<std::pair<double,double> > tlc(3);
  tlc[0].first = 1; // Line numbers
  tlc[1].first = 2;
  tlc[2].first = 4;
  tlc[0].second = 3; // Time
  tlc[1].second = 4;
  tlc[2].second = 4.5;

  TLCTimeInterpolation interp( tlc, -4 );

  // Assert data points are still the same
  EXPECT_NEAR( -1, interp(1), 1e-12 );
  EXPECT_NEAR( 0, interp(2), 1e-12 );
  EXPECT_NEAR( 0.5, interp(4), 1e-12 );

  // Check that interpolation is happening between points
  EXPECT_NEAR( -0.5, interp(1.5), 1e-12 );
  EXPECT_NEAR( 0.25, interp(3), 1e-12 );
  EXPECT_NEAR( -2, interp(0), 1e-12 );
  EXPECT_NEAR( 1, interp(6), 1e-12 );
}
