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

  // Testing delete control point
  cnet.delete_control_point(2);
  ASSERT_EQ( cnet.size(), 3u );
  EXPECT_EQ( cnet[0].size(), 1u );
  EXPECT_EQ( cnet[1].size(), 2u );
  EXPECT_EQ( cnet[2].size(), 4u );

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
  
  // Do some basic read/write checks.
  std::string test_file = "temp_ControlNetwork.csv";

  EXPECT_NO_THROW(cnet.write_csv(test_file));

  ControlNetwork cnet_csv(test_file);
  //cnet_csv.read_csv(test_file);
  EXPECT_NO_THROW(cnet_csv.read_csv(test_file));

  EXPECT_EQ(cnet.num_tie_points(),            cnet_csv.num_tie_points());
  EXPECT_EQ(cnet.num_ground_control_points(), cnet_csv.num_ground_control_points());
  
  //EXPECT_NO_THROW(cnet_csv.write_in_gcp_format("cn_gcp_format.csv", cartography::Datum("WGS84")));
  cnet_csv.add_image_name("image1.tif");
  cnet_csv.add_image_name("image2.tif");
  cnet_csv.add_image_name("image3.tif");
  cnet_csv.add_image_name("image4.tif");
  cnet_csv.write_in_gcp_format("cn_gcp_format.csv", cartography::Datum("WGS84"));
  
  
  // Check clearing
  cnet[0].clear();
  ASSERT_EQ( cnet[0].size(), 0u );
  cnet.clear();
  ASSERT_EQ( cnet.size(), 0u );
}
