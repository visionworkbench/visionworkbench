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


#include <gtest/gtest_VW.h>
#include <test/Helpers.h>
#include <vw/InterestPoint/InterestPoint.h>

using namespace vw;
using namespace vw::ip;
using namespace vw::test;

TEST( InterestData, VWIP_IO_Loop ) {
  InterestPointList ip;
  for ( uint32 i = 0; i < 5; i++ ) {
    ip.push_back( InterestPoint( 2*i, 2*i+5, 1.0, -i, i, true, 5 ) );
    ip.back().descriptor = Vector3(5,6,i);
  }

  UnlinkName vwip_file( "monkey.vwip" );
  write_binary_ip_file( vwip_file, ip );
  std::vector<InterestPoint> result = read_binary_ip_file( vwip_file );

  InterestPointList::iterator ipiter = ip.begin();
  ASSERT_EQ( 5u, result.size() );
  for ( uint32 i = 0; i < 5; i++ ) {
    EXPECT_EQ( ipiter->x, result[i].x );
    EXPECT_EQ( ipiter->y, result[i].y );
    EXPECT_EQ( ipiter->scale, result[i].scale );
    EXPECT_EQ( ipiter->ix, result[i].ix );
    EXPECT_EQ( ipiter->iy, result[i].iy );
    EXPECT_EQ( ipiter->orientation, result[i].orientation );
    EXPECT_EQ( ipiter->interest, result[i].interest );
    EXPECT_EQ( ipiter->polarity, result[i].polarity );
    EXPECT_EQ( ipiter->octave, result[i].octave );
    EXPECT_EQ( ipiter->scale_lvl, result[i].scale_lvl );

    ASSERT_EQ( ipiter->size(), result[i].size() );
    EXPECT_VECTOR_FLOAT_EQ( ipiter->descriptor, result[i].descriptor );

    ipiter++;
  }
}

TEST( InterestData, MATCH_IO_Loop ) {
  std::vector<InterestPoint> ip1, ip2;
  ip1.reserve( 5 );
  ip2.reserve( 5 );
  for ( uint32 i = 0; i < 5; i++ ) {
    ip1.push_back( InterestPoint( 2*i, 2*i+5, 1.0, -i, i, true, 5 ) );
    ip1.back().descriptor = Vector3(5,6,i);
    ip2.push_back( InterestPoint( 20-2*i, i, 0.5, i, 5-i, false, 6 ) );
    ip2.back().descriptor = Vector3(7,i,2);
  }

  UnlinkName match_file( "monkey.match" );
  write_binary_match_file( match_file, ip1, ip2 );
  std::vector<InterestPoint> result1, result2;
  read_binary_match_file( match_file, result1, result2 );

  ASSERT_EQ( 5u, result1.size() );
  ASSERT_EQ( 5u, result2.size() );
  std::vector<InterestPoint>::iterator ip1iter = ip1.begin(),
    ip2iter = ip2.begin();
  for ( uint32 i = 0; i < 5; i++ ) {
    EXPECT_EQ( ip1iter->x, result1[i].x );
    EXPECT_EQ( ip1iter->y, result1[i].y );
    EXPECT_EQ( ip1iter->scale, result1[i].scale );
    EXPECT_EQ( ip1iter->ix, result1[i].ix );
    EXPECT_EQ( ip1iter->iy, result1[i].iy );
    EXPECT_EQ( ip1iter->orientation, result1[i].orientation );
    EXPECT_EQ( ip1iter->interest, result1[i].interest );
    EXPECT_EQ( ip1iter->polarity, result1[i].polarity );
    EXPECT_EQ( ip1iter->octave, result1[i].octave );
    EXPECT_EQ( ip1iter->scale_lvl, result1[i].scale_lvl );

    EXPECT_EQ( ip2iter->x, result2[i].x );
    EXPECT_EQ( ip2iter->y, result2[i].y );
    EXPECT_EQ( ip2iter->scale, result2[i].scale );
    EXPECT_EQ( ip2iter->ix, result2[i].ix );
    EXPECT_EQ( ip2iter->iy, result2[i].iy );
    EXPECT_EQ( ip2iter->orientation, result2[i].orientation );
    EXPECT_EQ( ip2iter->interest, result2[i].interest );
    EXPECT_EQ( ip2iter->polarity, result2[i].polarity );
    EXPECT_EQ( ip2iter->octave, result2[i].octave );
    EXPECT_EQ( ip2iter->scale_lvl, result2[i].scale_lvl );

    ASSERT_EQ( ip1iter->size(), result1[i].size() );
    EXPECT_VECTOR_FLOAT_EQ( ip1iter->descriptor, result1[i].descriptor );

    ASSERT_EQ( ip2iter->size(), result2[i].size() );
    EXPECT_VECTOR_FLOAT_EQ( ip2iter->descriptor, result2[i].descriptor );
    
    ip1iter++; ip2iter++;
  }
}
