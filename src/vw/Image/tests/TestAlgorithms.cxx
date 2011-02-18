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

TEST( Algorithms, ImageBlocks ) {
  ImageView<uint8> bin1(5,5);
  std::vector<BBox2i> bout1 = image_blocks(bin1,2,2);
  EXPECT_EQ( 9u, bout1.size() );
  EXPECT_EQ( 2, bout1.front().width() );
  EXPECT_EQ( 2, bout1.front().height() );
  EXPECT_EQ( 1, bout1.back().width() );
  EXPECT_EQ( 1, bout1.back().height() );

  ImageView<uint8> bin2(9,10);
  std::vector<BBox2i> bout2 = image_blocks(bin2,3,4);
  EXPECT_EQ( 9u, bout2.size() );
  EXPECT_EQ( 3, bout2.front().width() );
  EXPECT_EQ( 4, bout2.front().height() );
  EXPECT_EQ( 3, bout2.back().width() );
  EXPECT_EQ( 2, bout2.back().height() );
}

TEST( Algorithms, BlobIndex ) {
  typedef PixelMask<uint8> MPx;
  ImageView<MPx> im(7,5);
  fill( crop(im,1,1,3,1), MPx(255) ); // Blob 1
  fill( crop(im,2,3,2,1), MPx(255) ); // Blob 2
  fill( crop(im,5,1,1,3), MPx(255) );
  im(6,2) = MPx(5); // Blob 3

  ImageView<uint8> idx =  blob_index( im );
  EXPECT_EQ( 3, int(max_channel_value(idx)) );
  EXPECT_EQ( idx(1,1), idx(2,1) ); // Verify blob 1
  EXPECT_NE( idx(0,0), idx(1,1) );
  EXPECT_EQ( idx(2,4), idx(3,4) ); // Verify blob 2
  EXPECT_NE( idx(1,1), idx(2,4) );
  EXPECT_EQ( idx(5,1), idx(6,2) ); // Verify blob 3
  EXPECT_NE( idx(6,2), idx(1,1) );
}
