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


// TestBoxFilter.h
#include <gtest/gtest_VW.h>

#include <vw/InterestPoint/IntegralImage.h>
#include <vw/InterestPoint/BoxFilter.h>
#include <vw/Image/ImageView.h>

using namespace vw;
using namespace vw::ip;

TEST( BoxFilter, ApplyBoxFilter ) {
  ImageView<float> image(3,3);

  float count = 0;
  for ( uint32 i = 0; i < 3; i++ ) {
    for ( uint32 j = 0; j < 3; j++ ) {
      image(i,j) = count;
      count++;
    }
  }

  ImageView<float> integral = IntegralImage( image );
  EXPECT_EQ( 4, integral.cols() );
  EXPECT_EQ( 4, integral.rows() );

  BoxFilter filter;
  filter.resize(2);
  filter[0].start = Vector2i(0,0);
  filter[0].size = Vector2i(1,1);
  filter[0].weight = -9;
  filter[1].start = Vector2i(-1,-1);
  filter[1].size = Vector2i(3,3);
  filter[1].weight = 1;

  float response = apply_box_filter_at_point( integral.origin().advance(1,1),
                                              filter );
  EXPECT_NEAR( 0, response, 1e-5 );
}

TEST( BoxFilter, FilterView ) {
  ImageView<float> image(4,4);
  float count = 0;
  for ( uint8 i = 0; i < 4; i++ ) {
    for ( uint8 j = 0; j < 4; j++ ) {
      image(i,j) = count;
      count++;
    }
  }
  image(3,3) = 20;
  ImageView<float> integral = IntegralImage( image );
  EXPECT_EQ( 5, integral.cols() );
  EXPECT_EQ( 5, integral.rows() );

  BoxFilter filter;
  filter.resize(2);
  filter[0].start = Vector2i(0,0);
  filter[0].size = Vector2i(1,1);
  filter[0].weight = -9;
  filter[1].start = Vector2i(-1,-1);
  filter[1].size = Vector2i(3,3);
  filter[1].weight = 1;

  ImageView<float> applied = box_filter( integral, filter );
  EXPECT_EQ( 4, applied.cols() );
  EXPECT_EQ( 4, applied.rows() );

  EXPECT_NEAR( 0, applied(1,1), 1e-5 );
  EXPECT_NEAR( 5, applied(2,2), 1e-5 );
  EXPECT_NEAR( 0, applied(1,2), 1e-5 );
  EXPECT_NEAR( 0, applied(2,1), 1e-5 );
  EXPECT_NEAR( 0, applied(0,0), 1e-5 );
  EXPECT_NEAR( 0, applied(3,3), 1e-5 );
}
