// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// Image/tests/TestAlgorithms.h
#include <gtest/gtest.h>

#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>

#include <test/Helpers.h>

using namespace vw;

typedef PixelRGB<double> PixelT;

TEST( Algorithms, Fill ) {
  PixelT fillval(1,2,3);
  ImageView<PixelT> fl(2,1);
  fill( fl, fillval );

  EXPECT_PIXEL_EQ( fl(0,0), fillval );
  EXPECT_PIXEL_EQ( fl(1,0), fillval );
}

TEST( Algorithms, Clamp ) {
  ImageView<PixelT> im(2,1);
  im(0,0) = PixelT(0.5,0.9,1.5);
  im(1,0) = PixelT(-1,0,1);

  ImageView<PixelT> c1 = clamp(im,-0.2,0.2);
  EXPECT_PIXEL_EQ( c1(0,0), PixelT(0.2, 0.2, 0.2) );
  EXPECT_PIXEL_EQ( c1(1,0), PixelT(-0.2, 0.0, 0.2) );

  ImageView<PixelT> c2 = clamp(im,0.8);
  EXPECT_PIXEL_EQ( c2(0,0), PixelT(0.5, 0.8, 0.8) );
  EXPECT_PIXEL_EQ( c2(1,0), PixelT(0.0, 0.0, 0.8) );

  ImageView<PixelT> c3 = clamp(im);
  EXPECT_PIXEL_EQ( c3(0,0), PixelT(0.5, 0.9, 1.0) );
  EXPECT_PIXEL_EQ( c3(1,0), PixelT(0.0, 0.0, 1.0) );
}

TEST( Algorithms, Normalize ) {
  ImageView<PixelT > im(2,1);
  im(0,0) = PixelT(0.5,1.5,2.5);
  im(1,0) = PixelT(-1.5,0,1);

  ImageView<PixelT > n1 = normalize(im,-0.2,0.2);
  EXPECT_PIXEL_NEAR( n1(0,0), PixelT(0.0, 0.1, 0.2), 1e-7 );
  EXPECT_PIXEL_NEAR( n1(1,0), PixelT(-0.2, -0.05, 0.05), 1e-7 );

  ImageView<PixelT > n2 = normalize(im,0.8);
  EXPECT_PIXEL_NEAR( n2(0,0), PixelT(0.4, 0.6, 0.8), 1e-7 );
  EXPECT_PIXEL_NEAR( n2(1,0), PixelT(0.0, 0.3, 0.5), 1e-7 );

  ImageView<PixelT > n3 = normalize(im);
  EXPECT_PIXEL_NEAR( n3(0,0), PixelT(0.5, 0.75, 1.0), 1e-7 );
  EXPECT_PIXEL_NEAR( n3(1,0), PixelT(0.0, 0.375, 0.625), 1e-7 );
}

TEST( Algorithms, Threshold ) {
  ImageView<PixelT > im(2,1);
  im(0,0) = PixelT(0.5,1.5,2.5);
  im(1,0) = PixelT(-1.5,0,1);

  ImageView<PixelT > t1 = threshold(im,0.5,-0.2,0.2);
  EXPECT_PIXEL_EQ( t1(0,0), PixelT( -0.2, 0.2, 0.2 ) );
  EXPECT_PIXEL_EQ( t1(1,0), PixelT( -0.2, -0.2, 0.2 ) );

  ImageView<PixelT > t2 = threshold(im,0.6,0.8);
  EXPECT_PIXEL_EQ( t2(0,0), PixelT( 0.0, 0.8, 0.8 ) );
  EXPECT_PIXEL_EQ( t2(1,0), PixelT( 0.0, 0.0, 0.8 ) );

  ImageView<PixelT > t3 = threshold(im,0.6);
  EXPECT_PIXEL_EQ( t3(0,0), PixelT( 0.0, 1.0, 1.0 ) );
  EXPECT_PIXEL_EQ( t3(1,0), PixelT( 0.0, 0.0, 1.0 ) );

  ImageView<PixelT > t4 = threshold(im);
  EXPECT_PIXEL_EQ( t4(0,0), PixelT( 1.0, 1.0, 1.0 ) );
  EXPECT_PIXEL_EQ( t4(1,0), PixelT( 0.0, 0.0, 1.0 ) );
}
