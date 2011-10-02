// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <gtest/gtest.h>
#include <test/Helpers.h>

#include <vw/Stereo/rewrite/SubpixelView.h>
#include <vw/Image.h>
#include <boost/foreach.hpp>

#include <boost/random/linear_congruential.hpp>

using namespace vw;
using namespace vw::stereo::rewrite;

template <int32 istretch>
class SubPixelCorrelateTest : public ::testing::Test {
  const int32 IMAGE_SIZE, HALF_IMAGE_SIZE;

public:
  SubPixelCorrelateTest() : IMAGE_SIZE(100), HALF_IMAGE_SIZE(50) {}

protected:

  void SetUp() {
    stretch = float32(istretch)/100;

    boost::rand48 gen(10);
    image1 = transform(channel_cast_rescale<uint8>(uniform_noise_view( gen, IMAGE_SIZE, IMAGE_SIZE )),
                       AffineTransform(Matrix2x2(3,0,0,3),Vector2()),
                       ZeroEdgeExtension(), BicubicInterpolation());
    translation = HALF_IMAGE_SIZE-HALF_IMAGE_SIZE*stretch;
    image2 = transform(image1, AffineTransform(Matrix2x2(stretch,0,0,1),
                                               Vector2(translation,0) ),
                       ZeroEdgeExtension(), BicubicInterpolation());

    starting_disp.set_size(IMAGE_SIZE,IMAGE_SIZE);
    for ( int32 i = 0; i < IMAGE_SIZE ; i++ ) {
      int32 disparity =
        boost::numeric_cast<int32>(stretch * i + translation - i);
      for ( int32 j = 0; j < IMAGE_SIZE; j++ )
        starting_disp(i,j) = disparity;
    }
  }

  template <class ViewT>
  double check_error( ImageViewBase<ViewT> const& input,
                      int32& invalid_count ) {
    ViewT const& disparity = input.impl();
    double error = 0;
    for ( int32 i = 0; i < IMAGE_SIZE; i++ ) {
      float expected = stretch * float(i) + translation - i;
      for ( int32 j = 0; j < IMAGE_SIZE; j++ ) {
        error += disparity(i,j)[1] + fabs(disparity(i,j)[0] - expected);
        if ( !is_valid(disparity(i,j)) )
          invalid_count++;
      }
    }
    return error / (double(IMAGE_SIZE)*double(IMAGE_SIZE));
  }

  float32 stretch, translation;
  ImageView<uint8> image1, image2;
  ImageView<PixelMask<Vector2f> > starting_disp;
};

TEST( ParabolaSubpixel, NullTest ) {
  ImageView<PixelMask<Vector2i> > disparity(5,5);
  fill( disparity, PixelMask<Vector2i>(Vector2i(1,1)) );
  ImageView<float> left(5,5), right(5,5);
  fill( left, 0.5 );
  fill( right, 0.6 );

  ImageView<PixelMask<Vector2f> > fdisparity =
    parabola_subpixel( disparity, left, right,
                       preprocessing<NULLOP>(),
                       Vector2i(3,3) );
  EXPECT_EQ( fdisparity.cols(), 5 );
  EXPECT_EQ( fdisparity.rows(), 5 );
  BOOST_FOREACH( PixelMask<Vector2f> const& fdisp, fdisparity ) {
    EXPECT_TRUE( is_valid( fdisp ) );
    EXPECT_VECTOR_NEAR( fdisp.child(), Vector2f(1,1), 0.1 );
  }

  fdisparity =
    parabola_subpixel( disparity, left, right,
                       preprocessing<LAPLACIAN_OF_GAUSSIAN>(1.4),
                       Vector2i(3,3) );
  EXPECT_EQ( fdisparity.cols(), 5 );
  EXPECT_EQ( fdisparity.rows(), 5 );
  BOOST_FOREACH( PixelMask<Vector2f> const& fdisp, fdisparity ) {
    EXPECT_TRUE( is_valid( fdisp ) );
    EXPECT_VECTOR_NEAR( fdisp.child(), Vector2f(1,1), 0.1 );
  }
}

typedef SubPixelCorrelateTest<95> SubPixelCorrelate95Test;
typedef SubPixelCorrelateTest<90> SubPixelCorrelate90Test;
typedef SubPixelCorrelateTest<80> SubPixelCorrelate80Test;
typedef SubPixelCorrelateTest<70> SubPixelCorrelate70Test;

TEST_F( SubPixelCorrelate95Test, Parabola ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    parabola_subpixel( starting_disp, image1, image2,
                       preprocessing<LAPLACIAN_OF_GAUSSIAN>(1.4),
                       Vector2i(7,7) );

  int32 invalid_count = 0;
  double error = check_error( disparity_map, invalid_count );
  EXPECT_LT(error, 0.6);
  EXPECT_LE(invalid_count, 0);
}
