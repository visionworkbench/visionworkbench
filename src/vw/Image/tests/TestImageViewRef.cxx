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


// TestImageViewRef.h
#include <gtest/gtest_VW.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Interpolation.h>

using namespace vw;

TEST( ImageViewRef, Construct ) {
  const int cols=3, rows=2;
  ImageView<float> image(cols,rows);
  for( int r=0; r<rows; ++r )
    for( int c=0; c<cols; ++c )
      image(c,r) = (float)(r*cols+c);

  ImageViewRef<float> ref = image;
  ASSERT_EQ( ref.cols(), image.cols() );
  ASSERT_EQ( ref.rows(), image.rows() );
  ASSERT_EQ( ref.planes(), image.planes() );

  // Test pixel indexing
  for( int r=0; r<rows; ++r )
    for( int c=0; c<cols; ++c )
      EXPECT_EQ( ref(c,r), (float)(r*cols+c) );

  // Test full rasterization: optimized case
  ImageView<float> im2 = ref;
  for( int r=0; r<rows; ++r )
    for( int c=0; c<cols; ++c )
      EXPECT_EQ( im2(c,r), ref(c,r) );

  // Test full rasterization: general case
  ImageView<double> im3 = ref;
  for( int r=0; r<rows; ++r )
    for( int c=0; c<cols; ++c )
      EXPECT_NEAR( im3(c,r), ref(c,r), 1e-8 );

  // Test accessor / generic rasterization
  ImageView<float> im4(cols,rows);
  vw::rasterize( ref, im4, BBox2i(0,0,cols,rows) );
  for( int r=0; r<rows; ++r )
    for( int c=0; c<cols; ++c )
      EXPECT_EQ( im4(c,r), ref(c,r) );

  // Test partial rasterization
  ImageView<float> im5(cols-1,rows-1);
  ref.rasterize( im5, BBox2i(1,1,cols-1,rows-1) );
  for( int r=0; r<rows-1; ++r )
    for( int c=0; c<cols-1; ++c )
      EXPECT_EQ( im5(c,r), ref(c+1,r+1) );

  // Test iterator
  int val=0;
  for( ImageViewRef<float>::iterator i=ref.begin(), end=ref.end();
       i!=end; ++i, ++val )
    EXPECT_EQ( *i, (float)(val) );
}

TEST( ImageViewRef, FloatAccess ) {
  ImageView<float> source(2,1); source(0,0) = 0; source(1,0) = 1;
  ImageViewRef<float> ref = interpolate( source );

  EXPECT_EQ( ref(double(0.5),double(0.0)), 0.5 );
  EXPECT_EQ( ref(double(0.5),double(0.0),0), 0.5 );
  EXPECT_EQ( ref(double(0.5),int32(0)), 0.5 );
  EXPECT_EQ( ref(double(0.5),int32(0),0), 0.5 );
  EXPECT_EQ( ref(double(0.5),uint32(0)), 0.5 );
  EXPECT_EQ( ref(double(0.5),uint16(0)), 0.5 );
  EXPECT_EQ( ref(double(0.5),uint64(0)), 0.5 );
  EXPECT_EQ( ref(int32(0),int32(0)), 0 );
  EXPECT_EQ( ref(int32(0),int32(0),0), 0 );
  EXPECT_EQ( ref(uint32(0),int32(0)), 0 );
  EXPECT_EQ( ref(uint32(0),int32(0),0), 0 );
  EXPECT_EQ( ref(char(0),int32(0)), 0 );
  EXPECT_EQ( ref(char(0),int32(0),0), 0 );
}
