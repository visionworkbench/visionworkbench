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
    m_p.push_back( Vector3() ); m_p.push_back( Vector3(0.1,0.1,0.1) );
    m_v.push_back( Vector3(1,1,1) ); m_v.push_back( Vector3(1,1,1) );
    HermitePositionInterpolation interp( m_p, m_v, 2, 0.1 );
    EXPECT_VECTOR_NEAR( Vector3(0.05,0.05,0.05), interp(2.05), 1e-12 );
    EXPECT_VECTOR_NEAR( Vector3(), interp(2.00), 1e-12 );
    EXPECT_VECTOR_NEAR( Vector3(0.1,0.1,0.1), interp(2.1 - 1e-12), 1e-6 );
    EXPECT_VECTOR_NEAR( Vector3(0.01,0.01,0.01), interp(2.01), 1e-12 );
  }

  {
    // Hermite expects a constant acceleration
    // a = c0
    // v = c0 * t + c1
    // p = c0 * t^2 / 2 + c1 * t + c2
    Matrix<double,3,3> p_coeff;
    select_col(p_coeff,0) = Vector3(0.5,0.9,2.1); // Acceleration coeff
    select_col(p_coeff,1) = Vector3(-10,7.5,6.3);  // Velocity coeff
    select_col(p_coeff,2) = Vector3(500,0,-900);   // Position coeff

    Matrix<double,3,2> v_coeff;
    select_col(v_coeff,0) = select_col(p_coeff,0) * 2;
    select_col(v_coeff,1) = select_col(p_coeff,1);
    v_coeff(1,0) *= -1;

    m_p.clear(); m_v.clear();
    m_p.push_back( p_coeff * Vector3( 2 * 2, 2, 1 ) );
    m_p.push_back( p_coeff * Vector3( 2.1 * 2.1, 2.1, 1 ) );
    m_v.push_back( v_coeff * Vector2( 2, 1 ) );
    m_v.push_back( v_coeff * Vector2( 2.1, 1 ) );

    HermitePositionInterpolation interp( m_p, m_v, 2, 0.1 );
    EXPECT_VECTOR_NEAR( m_p[0], interp( 2 ), 1e-8 );
    EXPECT_VECTOR_NEAR( m_p[1], interp( 2.1 - 1e-12 ), 1e-8 );
    EXPECT_VECTOR_NEAR( p_coeff * Vector3( 2.05 * 2.05, 2.05, 1 ),
                        interp( 2.05 ), 1e-1 );
    EXPECT_VECTOR_NEAR( p_coeff * Vector3( 2.01 * 2.01, 2.01, 1 ),
                        interp( 2.01 ), 1e-1 );
    EXPECT_VECTOR_NEAR( p_coeff * Vector3( 2.09 * 2.09, 2.09, 1 ),
                        interp( 2.09 ), 1e-1 );
    EXPECT_THROW( interp(-.5), ArgumentErr );
    EXPECT_THROW( interp(1.5), ArgumentErr );
  }
}

TEST( Extrinsics, PiecewiseAPositionInterpolation ) {
  std::vector<Vector3> m_p, m_v;

  {
    m_p.push_back( Vector3() ); m_p.push_back( Vector3(0.1,0.1,0.1) );
    m_v.push_back( Vector3(1,1,1) ); m_v.push_back( Vector3(1,1,1) );
    PiecewiseAPositionInterpolation interp( m_p, m_v, 2, 0.1 );
    EXPECT_VECTOR_NEAR( Vector3(0.05,0.05,0.05), interp(2.05), 1e-12 );
    EXPECT_VECTOR_NEAR( Vector3(), interp(2.00), 1e-12 );
    EXPECT_VECTOR_NEAR( Vector3(0.1,0.1,0.1), interp(2.1 - 1e-12), 1e-6 );
    EXPECT_VECTOR_NEAR( Vector3(0.01,0.01,0.01), interp(2.01), 1e-12 );
  }

  {
    // Hermite expects a constant acceleration
    // a = c0
    // v = c0 * t + c1
    // p = c0 * t^2 / 2 + c1 * t + c2
    Matrix<double,3,3> p_coeff;
    select_col(p_coeff,0) = Vector3(0.5,0.9,2.1); // Acceleration coeff
    select_col(p_coeff,1) = Vector3(-10,7.5,6.3);  // Velocity coeff
    select_col(p_coeff,2) = Vector3(500,0,-900);   // Position coeff

    Matrix<double,3,2> v_coeff;
    select_col(v_coeff,0) = select_col(p_coeff,0) * 2;
    select_col(v_coeff,1) = select_col(p_coeff,1);
    //v_coeff(1,0) *= -1;

    m_p.clear(); m_v.clear();
    m_p.push_back( p_coeff * Vector3( 2 * 2, 2, 1 ) );
    m_p.push_back( p_coeff * Vector3( 2.1 * 2.1, 2.1, 1 ) );
    m_p.push_back( p_coeff * Vector3( 2.2 * 2.2, 2.2, 1 ) );
    m_v.push_back( v_coeff * Vector2( 2, 1 ) );
    m_v.push_back( v_coeff * Vector2( 2.1, 1 ) );
    m_v.push_back( v_coeff * Vector2( 2.2, 1 ) );

    PiecewiseAPositionInterpolation interp( m_p, m_v, 2, 0.1 );
    EXPECT_VECTOR_NEAR( m_p[0], interp( 2 ), 1e-12 );
    EXPECT_VECTOR_NEAR( m_p[1], interp( 2.1 - 1e-15), 1e-9 );
    EXPECT_VECTOR_NEAR( m_p[1], interp( 2.1 ), 1e-12 );
    EXPECT_VECTOR_NEAR( m_p[2], interp( 2.2 - 1e-15 ), 1e-9 );
    EXPECT_VECTOR_NEAR( p_coeff * Vector3( 2.05 * 2.05, 2.05, 1 ),
                        interp( 2.05 ), 1e-12 );
    EXPECT_VECTOR_NEAR( p_coeff * Vector3( 2.01 * 2.01, 2.01, 1 ),
                        interp( 2.01 ), 1e-12 );
    EXPECT_VECTOR_NEAR( p_coeff * Vector3( 2.09 * 2.09, 2.09, 1 ),
                        interp( 2.09 ), 1e-12 );
    EXPECT_THROW( interp(-.5), ArgumentErr );
    EXPECT_THROW( interp(1.5), ArgumentErr );
  }
}

TEST( Extrinsics, LinearPiecewisePositionInterpolation ) {
  std::vector<Vector3> m_p;
  m_p.push_back( Vector3() );
  m_p.push_back( Vector3(1,0,0) );
  m_p.push_back( Vector3(1,1,1) );

  LinearPiecewisePositionInterpolation acc( m_p, 4, 0.1 );
  EXPECT_VECTOR_NEAR( Vector3(), acc(4), 1e-12 );
  EXPECT_VECTOR_NEAR( Vector3(1,0,0), acc(4.1), 1e-12 );
  EXPECT_VECTOR_NEAR( Vector3(1,1,1), acc(4.2 - 1e-12), 1e-6 );

  EXPECT_VECTOR_NEAR( Vector3(0.5,0,0), acc(4.05), 1e-12 );
  EXPECT_VECTOR_NEAR( Vector3(1,0.5,0.5), acc(4.15), 1e-12 );
  EXPECT_VECTOR_NEAR( Vector3(1,0.1,0.1), acc(4.11), 1e-12 );
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
