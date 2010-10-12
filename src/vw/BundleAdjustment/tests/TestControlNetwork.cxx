// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <sstream>
#include <vw/BundleAdjustment/ControlNetwork.h>

#include <test/Helpers.h>

using namespace vw;
using namespace vw::ba;

TEST( ControlNetwork, Construction ) {

  ControlNetwork cnet( "TestCNET" );

  for ( uint32 i = 0; i < 4; i++ ) {
    ControlPoint cpoint;

    for ( uint32 j = 0; j < i+1; j++ ) {
      ControlMeasure cm( 100, 100, 1, 1, j );
      cpoint.add_measure( cm );
    }
    cnet.add_control_point( cpoint );
  }

  // Checking to see that sizes match
  ASSERT_EQ( cnet.size(), 4u );
  ASSERT_EQ( cnet[0].size(), 1u );
  ASSERT_EQ( cnet[1].size(), 2u );
  ASSERT_EQ( cnet[2].size(), 3u );
  ASSERT_EQ( cnet[3].size(), 4u );
  EXPECT_VECTOR_DOUBLE_EQ( cnet[0][0].position(),
                           Vector2(100,100) );
  EXPECT_VECTOR_DOUBLE_EQ( cnet[0][0].dominant(),
                           Vector2(100,100) );
  EXPECT_VECTOR_DOUBLE_EQ( cnet[0][0].sigma(),
                           Vector2(1,1) );
  EXPECT_NEAR( cnet[0][0].sigma_magnitude(),
               sqrt(2), 1e-6 );
  EXPECT_VECTOR_DOUBLE_EQ( cnet[0].position(),
                           Vector3() );
  EXPECT_EQ( cnet.num_tie_points(), 4u );
  EXPECT_EQ( cnet.num_ground_control_points(), 0u);

  // Checking that default settings were applied
  EXPECT_TRUE( cnet[0][0].is_pixels_dominant() );
  EXPECT_FALSE( cnet[0][0].ignore() );
  EXPECT_EQ( cnet.type(), ControlNetwork::ImageToImage );

  // Checking that ostreamstring still works
  std::ostringstream ostr;
  ostr << cnet[0][0];
  EXPECT_GT( ostr.str().size(), 3u );
  ostr.clear();
  ostr << cnet[0];
  EXPECT_GT( ostr.str().size(), 20u );
  ostr.clear();
  ostr << cnet;
  EXPECT_GT( ostr.str().size(), 100u );

  // Check clearing
  cnet[3].clear();
  ASSERT_EQ( cnet[3].size(), 0u );
  cnet.clear();
  ASSERT_EQ( cnet.size(), 0u );
}
