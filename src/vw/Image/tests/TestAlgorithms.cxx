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


// Image/tests/TestAlgorithms.h
#include <gtest/gtest_VW.h>

#include <vw/Image/Algorithms.h>
#include <vw/Image/Grassfire.h>
#include <vw/Image/WindowAlgorithms.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/Manipulation.h>

#include <test/Helpers.h>

using namespace vw;

typedef PixelRGB<double> PixelT;

TEST( Algorithms, Fill ) {
  PixelT fillval(1,2,3);
  ImageView<PixelT> fl(2,1);
  fill( fl, fillval );

  EXPECT_PIXEL_EQ( fillval, fl(0,0) );
  EXPECT_PIXEL_EQ( fillval, fl(1,0) );
}

TEST( Algorithms, Clamp ) {
  ImageView<PixelT> im(2,1);
  im(0,0) = PixelT(0.5,0.9,1.5);
  im(1,0) = PixelT(-1,0,1);

  ImageView<PixelT> c1 = clamp(im,-0.2,0.2);
  EXPECT_PIXEL_EQ( PixelT(0.2, 0.2, 0.2), c1(0,0) );
  EXPECT_PIXEL_EQ( PixelT(-0.2, 0.0, 0.2), c1(1,0) );

  ImageView<PixelT> c2 = clamp(im,0.8);
  EXPECT_PIXEL_EQ( PixelT(0.5, 0.8, 0.8), c2(0,0) );
  EXPECT_PIXEL_EQ( PixelT(0.0, 0.0, 0.8), c2(1,0) );

  ImageView<PixelT> c3 = clamp(im);
  EXPECT_PIXEL_EQ( PixelT(0.5, 0.9, 1.0), c3(0,0) );
  EXPECT_PIXEL_EQ( PixelT(0.0, 0.0, 1.0), c3(1,0) );
}

TEST( Algorithms, Normalize ) {
  ImageView<PixelT > im(2,1);
  im(0,0) = PixelT(0.5,1.5,2.5);
  im(1,0) = PixelT(-1.5,0,1);

  ImageView<PixelT > n1 = normalize(im,-0.2,0.2);
  EXPECT_PIXEL_NEAR( PixelT(0.0, 0.1, 0.2), n1(0,0), 1e-7 );
  EXPECT_PIXEL_NEAR( PixelT(-0.2, -0.05, 0.05), n1(1,0), 1e-7 );

  ImageView<PixelT > n2 = normalize(im,0.8);
  EXPECT_PIXEL_NEAR( PixelT(0.4, 0.6, 0.8), n2(0,0), 1e-7 );
  EXPECT_PIXEL_NEAR( PixelT(0.0, 0.3, 0.5), n2(1,0), 1e-7 );

  ImageView<PixelT > n3 = normalize(im);
  EXPECT_PIXEL_NEAR( PixelT(0.5, 0.75, 1.0), n3(0,0), 1e-7 );
  EXPECT_PIXEL_NEAR( PixelT(0.0, 0.375, 0.625), n3(1,0), 1e-7 );

  ImageView<PixelT > n4 = normalize(im,0,1,-2,2);
  EXPECT_PIXEL_NEAR( PixelT(0, 4, 8), n4(0,0), 1e-7 );
  EXPECT_PIXEL_NEAR( PixelT(-8,-2,2), n4(1,0), 1e-7 );

  // Int16 version
  ImageView<int16> im2(4,1);
  im2(0,0) = -30000; im2(1,0) = -15000;
  im2(2,0) =  15000; im2(3,0) =  30000;
  ImageView<int16> out2 = normalize(im2,-30000,30000,0,30000);
  EXPECT_EQ( int16(    0), out2(0,0));
  EXPECT_EQ( int16( 7500), out2(1,0));
  EXPECT_EQ( int16(22500), out2(2,0));
  EXPECT_EQ( int16(30000), out2(3,0));
}

TEST( Algorithms, Threshold ) {
  ImageView<PixelT > im(2,1);
  im(0,0) = PixelT(0.5,1.5,2.5);
  im(1,0) = PixelT(-1.5,0,1);

  ImageView<PixelT > t1 = threshold(im,0.5,-0.2,0.2);
  EXPECT_PIXEL_EQ( PixelT( -0.2, 0.2, 0.2 ), t1(0,0) );
  EXPECT_PIXEL_EQ( PixelT( -0.2, -0.2, 0.2 ), t1(1,0) );

  ImageView<PixelT > t2 = threshold(im,0.6,0.8);
  EXPECT_PIXEL_EQ( PixelT( 0.0, 0.8, 0.8 ), t2(0,0) );
  EXPECT_PIXEL_EQ( PixelT( 0.0, 0.0, 0.8 ), t2(1,0) );

  ImageView<PixelT > t3 = threshold(im,0.6);
  EXPECT_PIXEL_EQ( PixelT( 0.0, 1.0, 1.0 ), t3(0,0) );
  EXPECT_PIXEL_EQ( PixelT( 0.0, 0.0, 1.0 ), t3(1,0) );

  ImageView<PixelT > t4 = threshold(im);
  EXPECT_PIXEL_EQ( PixelT( 1.0, 1.0, 1.0 ), t4(0,0) );
  EXPECT_PIXEL_EQ( PixelT( 0.0, 0.0, 1.0 ), t4(1,0) );
}

TEST( Algorithms, NormalizeRetainAlpha ) {
  typedef PixelGrayA<float> Px;
  ImageView<Px> im(2,1);
  im(0,0) = Px( 1.0, 0.2 );
  im(1,0) = Px( 2.0, 0.8 );

  ImageView<Px> n1 = normalize_retain_alpha(im,0,2,0,1);
  EXPECT_PIXEL_NEAR( Px( 0.5, 0.2), n1(0,0), 1e-7 );
  EXPECT_PIXEL_NEAR( Px( 1.0, 0.8), n1(1,0), 1e-7 );

  // Int16 Version
  ImageView<PixelGrayA<int16> > im2(4,1);
  im2(0,0) = PixelGrayA<int16>(-30000,-3200);
  im2(1,0) = PixelGrayA<int16>(-15000,100);
  im2(2,0) = PixelGrayA<int16>( 15000,300);
  im2(3,0) = PixelGrayA<int16>( 30000,32000);
  ImageView<PixelGrayA<int16> > out2 =
    normalize_retain_alpha(im2,-30000,30000,0,30000);

  EXPECT_PIXEL_EQ( PixelGrayA<int16>(0,    -3200), out2(0,0) );
  EXPECT_PIXEL_EQ( PixelGrayA<int16>(7500, 100),   out2(1,0) );
  EXPECT_PIXEL_EQ( PixelGrayA<int16>(22500,300),   out2(2,0) );
  EXPECT_PIXEL_EQ( PixelGrayA<int16>(30000,32000), out2(3,0) );
}

TEST( Algorithms, Grassfirealpha ) {
  ImageView<uint8> im(5,5);
  fill(crop(im,1,1,3,3), 255);
  ImageView<uint8> g = grassfire(im);
  EXPECT_EQ( 0, g(0,0) );
  EXPECT_EQ( 1, g(1,1) );
  EXPECT_EQ( 2, g(2,2) );
  EXPECT_EQ( 1, g(3,3) );
  EXPECT_EQ( 0, g(4,4) );
  EXPECT_EQ( 1, g(1,2) );

  // Just make sure this compiles
  ImageView<PixelGray<float> > im2(5,5);
  fill(im2, 100.0);
  ImageView<PixelGray<float> > g2;
  grassfire(im2, g2);
}

TEST( Algorithms, Grassfire_no_border ) {
  ImageView<uint16> im(5,5);
  fill(im, 255);
  im(2, 2) = 0; // Only the center pixel is zero
  ImageView<uint16> g = grassfire(im, true);

  for (int r=0; r<5; ++r) {
    for (int c=0; c<5; ++c) {
      uint16 ans = abs(2-r) + abs(2-c);
      EXPECT_EQ( ans, g(c,r) );
    }
  }
}

TEST( Algorithms, CenterlineWeights ) {
  ImageView<uint8> image(5,5);
  fill(crop(image,1,1,3,3), 255);
  ImageView<double> weights;
  centerline_weights(image, weights);
  
  const double EPS = 1e-7;
  EXPECT_NEAR(1e-16, weights(0,0), EPS);
  EXPECT_NEAR(0.25,  weights(1,1), EPS);
  EXPECT_NEAR(1,     weights(2,2), EPS);
  EXPECT_NEAR(0.25,  weights(3,3), EPS);
  EXPECT_NEAR(1e-16, weights(4,4), EPS);
  EXPECT_NEAR(0.5,   weights(1,2), EPS);
}

TEST( Algorithms, BoundingBox ) {
  ImageView<PixelT> im(2,5);
  BBox2i b1 = bounding_box(im);
  EXPECT_EQ( 2, b1.width() );
  EXPECT_EQ( 5, b1.height() );
  EXPECT_VECTOR_FLOAT_EQ( Vector2i(), b1.min() );
  EXPECT_VECTOR_FLOAT_EQ( Vector2i(2,5), b1.max() );

  im.set_size(10,66);
  BBox2i b2 = bounding_box(im);
  EXPECT_EQ( 10, b2.width() );
  EXPECT_EQ( 66, b2.height() );
  EXPECT_VECTOR_FLOAT_EQ( Vector2i(), b2.min() );
  EXPECT_VECTOR_FLOAT_EQ( Vector2i(10,66), b2.max() );
}

TEST( Algorithms, NonZeroBoundingBox ) {
  ImageView<uint8> im(10,10);
  fill(crop(im,4,2,4,8),255);
  BBox2i b1 = nonzero_data_bounding_box(im);
  EXPECT_EQ( 4, b1.width() );
  EXPECT_EQ( 8, b1.height() );
  EXPECT_VECTOR_FLOAT_EQ( Vector2i(4,2), b1.min() );
}

TEST( Algorithms, TransparentOpaque ) {
  typedef PixelGrayA<float> Px;
  ImageView<Px> im1(2,1);
  im1(1,0) = Px(0,.5);
  EXPECT_FALSE( is_transparent(im1) );
  EXPECT_FALSE( is_opaque(im1) );

  ImageView<Px> im2(2,1);
  im2(0,0) = Px(.5,0); im2(1,0) = Px(.8,0);
  EXPECT_TRUE( is_transparent(im2) );
  EXPECT_FALSE( is_opaque(im2) );

  ImageView<Px> im3(2,1);
  im3(0,0) = Px(.7,1); im3(1,0) = Px(0,1);
  EXPECT_FALSE( is_transparent(im3) );
  EXPECT_TRUE( is_opaque(im3) );

  ImageView<Px> im4(2,1);
  im4(0,0) = Px(.7,0.5); im4(1,0) = Px(1.0,0.5);
  EXPECT_FALSE( is_transparent(im4) );
  EXPECT_FALSE( is_opaque(im4) );

  EXPECT_FALSE( is_transparent(im4(0,0)) );
  EXPECT_FALSE( is_opaque(im4(0,0)) );

  EXPECT_FALSE( is_transparent(im4(1,0)) );
  EXPECT_FALSE( is_opaque(im4(1,0)) );
}

TEST( Algorithms, SubdivideBBox ) {
  ImageView<uint8> bin1(5,5);
  std::vector<BBox2i> bout1 = subdivide_bbox(bin1,2,2);
  EXPECT_EQ( 9u, bout1.size() );
  EXPECT_EQ( 2, bout1.front().width() );
  EXPECT_EQ( 2, bout1.front().height() );
  EXPECT_EQ( 1, bout1.back().width() );
  EXPECT_EQ( 1, bout1.back().height() );

  ImageView<uint8> bin2(9,10);
  std::vector<BBox2i> bout2 = subdivide_bbox(bin2,3,4);
  EXPECT_EQ( 9u, bout2.size() );
  EXPECT_EQ( 3, bout2.front().width() );
  EXPECT_EQ( 4, bout2.front().height() );
  EXPECT_EQ( 3, bout2.back().width() );
  EXPECT_EQ( 2, bout2.back().height() );
}

TEST( Algorithms, MeanFillTransparent ) {
  typedef PixelGrayA<uint8> Px8;
  typedef PixelRGBA<float32> PxF;
  typedef PixelGray<float32> NoP;

  {
    ImageView<Px8> im(2,2);
    fill( im, Px8() );
    im(1,1) = Px8(128);
    ImageView<Px8> r = mean_fill_transparent( im );
    EXPECT_FALSE( is_transparent(r(1,1)) );
    EXPECT_TRUE( is_transparent(r(0,0)) );
    EXPECT_TRUE( is_transparent(r(1,0)) );
    EXPECT_TRUE( is_transparent(r(0,1)) );
    EXPECT_EQ( 128, non_alpha_channels(r(0,0)) );
    EXPECT_EQ( 128, non_alpha_channels(r(1,0)) );
    EXPECT_EQ( 128, non_alpha_channels(r(0,1)) );
  }

  {
    ImageView<PxF> im(2,2);
    fill( im, PxF() );
    im(1,1) = PxF(PixelRGB<float32>(0.5,0.1,0.2));
    ImageView<PxF> r = mean_fill_transparent( im );
    EXPECT_FALSE( is_transparent(r(1,1)) );
    EXPECT_TRUE( is_transparent(r(0,0)) );
    EXPECT_TRUE( is_transparent(r(1,0)) );
    EXPECT_TRUE( is_transparent(r(0,1)) );
    EXPECT_VW_EQ( PixelRGB<float32>(0.5,0.1,0.2),
                  non_alpha_channels(r(0,0)) );
    EXPECT_VW_EQ( PixelRGB<float32>(0.5,0.1,0.2),
                  non_alpha_channels(r(0,0)) );
    EXPECT_VW_EQ( PixelRGB<float32>(0.5,0.1,0.2),
                  non_alpha_channels(r(0,0)) );
  }

  {
    ImageView<NoP> im(2,2);
    fill( im, NoP() );
    im(1,1) = NoP(0.7);
    ImageView<NoP> r = mean_fill_transparent( im );
    EXPECT_VW_EQ( NoP(), r(0,0) );
    EXPECT_VW_EQ( NoP(), r(1,0) );
    EXPECT_VW_EQ( NoP(), r(0,1) );
    EXPECT_VW_EQ( NoP(0.7), r(1,1) );
  }
}

TEST( Algorithms, TwoThresholdFill ) {

  const float low_thresh  = 2.0;
  const float high_thresh = 6.0;

  // Create a test image
  ImageView<float> image(6, 6);
  image(0,0) = 0.0;  image(1,0) = 0.0;  image(2,0) = 0.0;  image(3,0) = 7.0;  image(4,0) = 8.0;  image(5,0) = 0.0;
  image(0,1) = 0.0;  image(1,1) = 7.0;  image(2,1) = 0.0;  image(3,1) = 0.0;  image(4,1) = 0.0;  image(5,1) = 0.0;
  image(0,2) = 0.0;  image(1,2) = 8.0;  image(2,2) = 0.0;  image(3,2) = 0.0;  image(4,2) = 3.0;  image(5,2) = 5.0;
  image(0,3) = 5.0;  image(1,3) = 4.0;  image(2,3) = 0.0;  image(3,3) = 0.0;  image(4,3) = 0.0;  image(5,3) = 2.0;
  image(0,4) = 3.0;  image(1,4) = 0.0;  image(2,4) = 5.0;  image(3,4) = 0.0;  image(4,4) = 0.0;  image(5,4) = 0.0;
  image(0,5) = 0.0;  image(1,5) = 3.0;  image(2,5) = 0.0;  image(3,5) = 0.0;  image(4,5) = 0.0;  image(5,5) = 0.0;

  // Run the algorithm
  ImageView<uint8> output = two_threshold_fill(image, 256, low_thresh, high_thresh);

  // Check results
  const uint8 ON  = 1;
  const uint8 OFF = 0;
  EXPECT_EQ(ON,  output(1,1));
  EXPECT_EQ(ON,  output(0,4));
  EXPECT_EQ(ON,  output(1,5));
  EXPECT_EQ(ON,  output(4,0));
  EXPECT_EQ(ON,  output(2,4));
  EXPECT_EQ(OFF, output(0,0));
  EXPECT_EQ(OFF, output(4,2));
  EXPECT_EQ(OFF, output(5,3));
  EXPECT_EQ(OFF, output(1,4));
  EXPECT_EQ(OFF, output(4,1));
}

TEST( Algorithms, MedianFilter ) {

  // Create a test image
  ImageView<float> image(6, 6);
  image(0,0) = 10.0;  image(1,0) = 10.0;  image(2,0) = 10.0;  image(3,0) =  3.0;  image(4,0) =  2.0;  image(5,0) = 10.0;
  image(0,1) = 10.0;  image(1,1) =  3.0;  image(2,1) = 10.0;  image(3,1) = 10.0;  image(4,1) = 10.0;  image(5,1) = 10.0;
  image(0,2) = 10.0;  image(1,2) =  2.0;  image(2,2) = 10.0;  image(3,2) = 10.0;  image(4,2) =  7.0;  image(5,2) =  5.0;
  image(0,3) =  5.0;  image(1,3) =  6.0;  image(2,3) = 10.0;  image(3,3) = 10.0;  image(4,3) = 10.0;  image(5,3) =  8.0;
  image(0,4) =  7.0;  image(1,4) = 10.0;  image(2,4) =  5.0;  image(3,4) = 10.0;  image(4,4) = 10.0;  image(5,4) = 10.0;
  image(0,5) = 10.0;  image(1,5) =  7.0;  image(2,5) = 10.0;  image(3,5) = 10.0;  image(4,5) = 10.0;  image(5,5) = 10.0;

  // Run the algorithm
  ImageView<float> output = median_filter_view(image, Vector2i(3,3));

  // Check results
  EXPECT_EQ(10, output(1,1));
  EXPECT_EQ( 6, output(0,2));
  EXPECT_EQ( 8, output(5,2));
  EXPECT_EQ(10, output(2,4));
  EXPECT_EQ( 7, output(1,4));
}

TEST( Algorithms, StdDevFilter ) {

  // Create a test image
  ImageView<float> image(6, 6);
  image(0,0) = 10.0;  image(1,0) = 10.0;  image(2,0) = 10.0;  image(3,0) =  3.0;  image(4,0) =  2.0;  image(5,0) = 10.0;
  image(0,1) = 10.0;  image(1,1) =  3.0;  image(2,1) = 10.0;  image(3,1) = 10.0;  image(4,1) = 10.0;  image(5,1) = 10.0;
  image(0,2) = 10.0;  image(1,2) =  2.0;  image(2,2) = 10.0;  image(3,2) = 10.0;  image(4,2) =  7.0;  image(5,2) =  5.0;
  image(0,3) =  5.0;  image(1,3) =  6.0;  image(2,3) = 10.0;  image(3,3) = 10.0;  image(4,3) = 10.0;  image(5,3) =  8.0;
  image(0,4) =  7.0;  image(1,4) = 10.0;  image(2,4) =  5.0;  image(3,4) = 10.0;  image(4,4) = 10.0;  image(5,4) = 10.0;
  image(0,5) = 10.0;  image(1,5) =  7.0;  image(2,5) = 10.0;  image(3,5) = 10.0;  image(4,5) = 10.0;  image(5,5) = 10.0;

  // Run the algorithm
  ImageView<float> output = stddev_filter_view(image, Vector2i(3,3));

  // Check results
  const double eps = 0.01;
  EXPECT_NEAR(3.13, output(1,1), eps);
  EXPECT_NEAR( 3.08, output(0,2), eps);
  EXPECT_NEAR( 1.97, output(5,2), eps);
  EXPECT_NEAR(1.945, output(2,4), eps);
  EXPECT_NEAR( 2.096, output(1,4), eps);
}
