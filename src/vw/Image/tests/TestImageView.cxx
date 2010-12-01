// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/ViewImageResource.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageIO.h>

using namespace vw;

TEST( ImageView, DefaultConstructor ) {
  ImageView<double> test_double;
  EXPECT_EQ(test_double.cols(), 0);
  EXPECT_EQ(test_double.rows(), 0);
  EXPECT_EQ(test_double.planes(), 0);
  ASSERT_EQ(test_double.data(), (double*)0);

  ImageView<PixelRGBA<vw::uint8> > test_rgba;
  EXPECT_EQ(test_rgba.cols(), 0);
  EXPECT_EQ(test_rgba.rows(), 0);
  EXPECT_EQ(test_rgba.planes(), 0);
  ASSERT_EQ(test_rgba.data(), (PixelRGBA<vw::uint8>*)0);
}

TEST( ImageView, ColsRowsConstructor ) {
  ImageView<double> test_double(3,4);
  ASSERT_TRUE( test_double );
  EXPECT_EQ(test_double.cols(), 3);
  EXPECT_EQ(test_double.rows(), 4);
  EXPECT_EQ(test_double.planes(), 1);
  ASSERT_NE(test_double.data(), (double*)0);

  ImageView<PixelRGBA<vw::uint8> > test_rgba(3,4);
  ASSERT_TRUE( test_rgba );
  EXPECT_EQ(test_rgba.cols(), 3);
  EXPECT_EQ(test_rgba.rows(), 4);
  EXPECT_EQ(test_rgba.planes(), 1);
  ASSERT_NE(test_rgba.data(), (PixelRGBA<vw::uint8>*)0);
}

TEST( ImageView, ColsRowsPlanesConstructor ) {
  ImageView<double> test_double(4,3,2);
  ASSERT_TRUE( test_double );
  EXPECT_EQ(test_double.cols(), 4);
  EXPECT_EQ(test_double.rows(), 3);
  EXPECT_EQ(test_double.planes(), 2);
  ASSERT_NE(test_double.data(), (double*)0);

  ImageView<PixelRGBA<vw::uint8> > test_rgba(4,3,2);
  ASSERT_TRUE( test_rgba );
  EXPECT_EQ(test_rgba.cols(), 4);
  EXPECT_EQ(test_rgba.rows(), 3);
  EXPECT_EQ(test_rgba.planes(), 2);
  ASSERT_NE(test_rgba.data(), (PixelRGBA<vw::uint8>*)0);
}

TEST( ImageView, CopyConstructor ) {
  ImageView<double> test_double(3,4);
  ImageView<double> test2_double( test_double );
  ASSERT_TRUE( test2_double );
  EXPECT_EQ(test2_double.cols(), 3);
  EXPECT_EQ(test2_double.rows(), 4);
  EXPECT_EQ(test2_double.planes(), 1);
  ASSERT_EQ(test2_double.data(),test_double.data());

  ImageView<PixelRGBA<vw::uint8> > test_rgba(3,4);
  ImageView<PixelRGBA<vw::uint8> > test2_rgba( test_rgba );
  ASSERT_TRUE( test2_rgba );
  EXPECT_EQ(test2_rgba.cols(), 3);
  EXPECT_EQ(test2_rgba.rows(), 4);
  EXPECT_EQ(test2_rgba.planes(), 1);
  ASSERT_EQ(test2_rgba.data(), test_rgba.data());
}

TEST( ImageView, SetSize ) {
  ImageView<double> test_double(3,4);
  ASSERT_TRUE( test_double );
  test_double.set_size(2,3);
  ASSERT_TRUE( test_double );
  EXPECT_EQ(test_double.cols(), 2);
  EXPECT_EQ(test_double.rows(), 3);
  EXPECT_EQ(test_double.planes(), 1);
  test_double.set_size(0,0);
  ASSERT_FALSE( test_double );
  EXPECT_EQ(test_double.cols(), 0);
  EXPECT_EQ(test_double.rows(), 0);
  ASSERT_EQ(test_double.data(), (double*)0);

  ImageView<PixelRGBA<vw::uint8> > test_rgba(3,4);
  ASSERT_TRUE( test_rgba );
  test_rgba.set_size(2,3);
  ASSERT_TRUE( test_rgba );
  EXPECT_EQ(test_rgba.cols(), 2);
  EXPECT_EQ(test_rgba.rows(), 3);
  EXPECT_EQ(test_rgba.planes(), 1);
  test_rgba.set_size(0,0);
  ASSERT_FALSE( test_rgba );
  EXPECT_EQ(test_rgba.cols(), 0);
  EXPECT_EQ(test_rgba.rows(), 0);
  ASSERT_EQ(test_rgba.data(), (PixelRGBA<vw::uint8>*)0);

  EXPECT_THROW(test_rgba.set_size(0,-1), ArgumentErr);
  EXPECT_THROW(test_rgba.set_size(0,std::numeric_limits<int32>::max()), ArgumentErr);
}

TEST( ImageView, Reset ) {
  ImageView<double> test_double(3,4);
  ASSERT_TRUE( test_double );
  test_double.reset();
  ASSERT_FALSE( test_double );
  EXPECT_EQ(test_double.cols(), 0);
  EXPECT_EQ(test_double.rows(), 0);
  EXPECT_EQ(test_double.planes(), 0);
  ASSERT_EQ(test_double.data(), (double*)0);

  ImageView<PixelRGBA<vw::uint8> > test_rgba(3,4);
  ASSERT_TRUE( test_rgba );
  test_rgba.reset();
  ASSERT_FALSE( test_rgba );
  EXPECT_EQ(test_rgba.cols(), 0);
  EXPECT_EQ(test_rgba.rows(), 0);
  EXPECT_EQ(test_rgba.planes(), 0);
  ASSERT_EQ(test_rgba.data(), (PixelRGBA<vw::uint8>*)0);
}

TEST( ImageView, Rasterization ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
  ImageView<double> im2(2,2);
  ASSERT_NO_THROW( im1.rasterize( im2, BBox2i(0,0,2,2) ) );
  EXPECT_EQ( im1(0,0), im2(0,0) );
  EXPECT_EQ( im1(1,0), im2(1,0) );
  EXPECT_EQ( im1(0,1), im2(0,1) );
  EXPECT_EQ( im1(1,1), im2(1,1) );
  EXPECT_NE( &im1(0,0), &im2(0,0) );
  EXPECT_NE( &im1(1,0), &im2(1,0) );
  EXPECT_NE( &im1(0,1), &im2(0,1) );
  EXPECT_NE( &im1(1,1), &im2(1,1) );
}

TEST( ImageView, DefaultRasterization1 ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
  ImageView<double> im2(2,2);
  ASSERT_NO_THROW( vw::rasterize( im1, im2, BBox2i(0,0,2,2) ) );
  EXPECT_EQ( im1(0,0), im2(0,0) );
  EXPECT_EQ( im1(1,0), im2(1,0) );
  EXPECT_EQ( im1(0,1), im2(0,1) );
  EXPECT_EQ( im1(1,1), im2(1,1) );
  EXPECT_NE( &im1(0,0), &im2(0,0) );
  EXPECT_NE( &im1(1,0), &im2(1,0) );
  EXPECT_NE( &im1(0,1), &im2(0,1) );
  EXPECT_NE( &im1(1,1), &im2(1,1) );
}

TEST( ImageView, DefaultRasterization2 ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
  ImageView<double> im2(2,2);
  ASSERT_NO_THROW( vw::rasterize( im1, im2 ) );
  EXPECT_EQ( im1(0,0), im2(0,0) );
  EXPECT_EQ( im1(1,0), im2(1,0) );
  EXPECT_EQ( im1(0,1), im2(0,1) );
  EXPECT_EQ( im1(1,1), im2(1,1) );
  EXPECT_NE( &im1(0,0), &im2(0,0) );
  EXPECT_NE( &im1(1,0), &im2(1,0) );
  EXPECT_NE( &im1(0,1), &im2(0,1) );
  EXPECT_NE( &im1(1,1), &im2(1,1) );
}

TEST( ImageView, DefaultRasterization3 ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
  ImageView<double> im2;
  ASSERT_NO_THROW( vw::rasterize( im1, im2, BBox2i(0,0,2,2) ) );
  EXPECT_EQ( im1(0,0), im2(0,0) );
  EXPECT_EQ( im1(1,0), im2(1,0) );
  EXPECT_EQ( im1(0,1), im2(0,1) );
  EXPECT_EQ( im1(1,1), im2(1,1) );
  EXPECT_NE( &im1(0,0), &im2(0,0) );
  EXPECT_NE( &im1(1,0), &im2(1,0) );
  EXPECT_NE( &im1(0,1), &im2(0,1) );
  EXPECT_NE( &im1(1,1), &im2(1,1) );
}

TEST( ImageView, DefaultRasterization4 ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
  ImageView<double> im2;
  ASSERT_NO_THROW( vw::rasterize( im1, im2 ) );
  EXPECT_EQ( im1(0,0), im2(0,0) );
  EXPECT_EQ( im1(1,0), im2(1,0) );
  EXPECT_EQ( im1(0,1), im2(0,1) );
  EXPECT_EQ( im1(1,1), im2(1,1) );
  EXPECT_NE( &im1(0,0), &im2(0,0) );
  EXPECT_NE( &im1(1,0), &im2(1,0) );
  EXPECT_NE( &im1(0,1), &im2(0,1) );
  EXPECT_NE( &im1(1,1), &im2(1,1) );
}

TEST( ImageView, PartialRasterization ) {
  ImageView<double> im1(2,2); im1(0,0)=1; im1(1,0)=2; im1(0,1)=3; im1(1,1)=4;
  ImageView<double> im2(1,2);
  ASSERT_NO_THROW( im1.rasterize( im2, BBox2i(1,0,1,2) ) );
  EXPECT_EQ( im1(1,0), im2(0,0) );
  EXPECT_EQ( im1(1,1), im2(0,1) );
  EXPECT_NE( &im1(1,0), &im2(0,0) );
  EXPECT_NE( &im1(1,1), &im2(0,1) );
  ImageView<double> im3(2,1);
  ASSERT_NO_THROW( im1.rasterize( im3, BBox2i(0,1,2,1) ) );
  EXPECT_EQ( im1(0,1), im3(0,0) );
  EXPECT_EQ( im1(1,1), im3(1,0) );
  EXPECT_NE( &im1(0,1), &im3(0,0) );
  EXPECT_NE( &im1(1,1), &im3(1,0) );
}

#define EXPECT_ITERATOR(I,C,R,P) \
  { EXPECT_GT( im.end() - I, 0 ); \
    EXPECT_EQ( &*I, &im(C,R,P) ); \
    EXPECT_EQ( I.col(), C ); \
    EXPECT_EQ( I.row(), R ); \
    EXPECT_EQ( I.plane(), P ); \
  }

TEST( ImageView, Iterator ) {
  ImageView<double> im(2,2,2);
  im(0,0,0)=1; im(1,0,0)=2; im(0,1,0)=3; im(1,1,0)=4;
  im(0,0,1)=5; im(1,0,1)=6; im(0,1,1)=7; im(1,1,1)=8;
  ImageView<double>::iterator i = im.begin();
  EXPECT_ITERATOR( i, 0, 0, 0 );
  ++i;
  EXPECT_ITERATOR( i, 1, 0, 0 );
  ++i;
  EXPECT_ITERATOR( i, 0, 1, 0 );
  ++i;
  EXPECT_ITERATOR( i, 1, 1, 0 );
  ++i;
  EXPECT_ITERATOR( i, 0, 0, 1 );
  ++i;
  EXPECT_ITERATOR( i, 1, 0, 1 );
  ++i;
  EXPECT_ITERATOR( i, 0, 1, 1 );
  ++i;
  EXPECT_ITERATOR( i, 1, 1, 1 );
  ++i;
  EXPECT_EQ( im.end() - i, 0 );
}

TEST( ImageView, BoundsCheck ) {
#if defined(VW_ENABLE_BOUNDS_CHECK) && (VW_ENABLE_BOUNDS_CHECK==1)
  ImageView<PixelRGB<vw::uint8> > test_img(2,2);

  // First, test to make sure that the operator() bounds checking is
  // working properly.
  PixelRGB<vw::uint8> test = test_img(0,0);
  ASSERT_NO_THROW ( test = test_img(1,0) );
  ASSERT_NO_THROW ( test = test_img(0,1) );
  ASSERT_NO_THROW ( test = test_img(1,1) );
  ASSERT_ANY_THROW ( test = test_img(0,2) );
  ASSERT_ANY_THROW ( test = test_img(2,0) );
  ASSERT_ANY_THROW ( test = test_img(-1,0) );
  ASSERT_ANY_THROW ( test = test_img(0,-1) );

  // Next, test to make sure that the MemoryStridingPixelAccessor
  // bounds checking is working properly.
  ImageView<PixelRGB<vw::uint8> >::pixel_accessor acc = test_img.origin();

  ASSERT_NO_THROW( test = *acc );

  // First, check to make sure that negative indices cause an error.
  acc.prev_col();
  ASSERT_ANY_THROW ( test = *acc );

  // Proceed through colums
  acc = test_img.origin();
  acc.next_col();
  ASSERT_NO_THROW ( test = *acc );
  acc.next_col();
  ASSERT_NO_THROW ( test = *acc );
  acc.next_col();
  ASSERT_NO_THROW ( test = *acc );
  acc.next_col();
  ASSERT_ANY_THROW ( test = *acc );

  // Proceed through rows
  acc = test_img.origin();
  acc.next_row();
  ASSERT_NO_THROW ( test = *acc );
  acc.next_row();
  ASSERT_ANY_THROW ( test = *acc );
#endif
}

TEST(ImageView, ViewImageResource) {
  typedef PixelRGB<uint8> Pixel;
  ImageView<Pixel> img1(2,2);
  img1(0,0) = Pixel(1,2,3);
  img1(0,1) = Pixel(5,6,7);
  img1(1,0) = Pixel(9,10,11);
  img1(1,1) = Pixel(13,14,15);

  ViewImageResource view(img1);

  ASSERT_EQ(img1.cols(),   view.cols());
  ASSERT_EQ(img1.rows(),   view.rows());
  ASSERT_EQ(img1.planes(), view.planes());

  ImageView<Pixel> img2;
  read_image(img2, view);

  ASSERT_EQ(img1.cols(),   img2.cols());
  ASSERT_EQ(img1.rows(),   img2.rows());
  ASSERT_EQ(img1.planes(), img2.planes());

  EXPECT_RANGE_EQ(img1.begin(), img1.end(), img2.begin(), img2.end());
}
