// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <test/Helpers.h>

#include <vw/Stereo/rewrite/PreFilter.h>
#include <vw/Stereo/rewrite/CorrelationView.h>
#include <boost/random/linear_congruential.hpp>

#include <vw/Image.h>

using namespace vw;
using namespace vw::stereo::rewrite;

template <typename PixelT>
class PyramidViewTest : public ::testing::Test {
protected:
  PyramidViewTest(){}

  typedef PixelT pixel_type;
  typedef ImageView<PixelT> image_type;
  typedef typename PixelChannelType<PixelT>::type channel_type;
  image_type input1, input2;
  Vector2i kernel_size;
  BBox2i search_volume;

  Matrix2x2 scale;
  Vector2 translation;

  virtual void SetUp() {
    boost::rand48 gen(10);
    input1 = ChannelRange<channel_type>::max()*uniform_noise_view( gen, 300, 200 );
    scale = Matrix2x2(.9,0,0,.95);
    translation = Vector2(150*.1,100*.05);
    input2 = transform(input1, AffineTransform(scale,translation),
                       ConstantEdgeExtension(), BicubicInterpolation());

    search_volume = BBox2( -1.5 * translation, 1.5*translation );
    kernel_size = Vector2i(7,7);

    // These settings should invoke 4 pyramid levels in the Pyramid
    // Correation View.
  }

  template <class ViewT>
  void check_error( ImageViewBase<ViewT> const& input,
                    float correct = 0.9,
                    float attempted = 0.9,
                    std::string const& mesg = "" ) {
    ViewT const& disparity_map = input.impl();
    int64 count_correct = 0;
    int64 count_valid = 0;
    for ( int32 i = 0; i < disparity_map.cols(); ++i )
      for ( int32 j = 0; j < disparity_map.rows(); ++j )
        if ( is_valid( disparity_map(i,j) ) ) {
          count_valid++;
          Vector2 objectf = scale*Vector2(i,j)+translation - Vector2(i,j);
          Vector2i objective( round(objectf[0]), round(objectf[1]) );
          if ( disparity_map(i,j).child() == objective )
            count_correct++;
        }
    EXPECT_GT( float(count_correct)/float(count_valid), correct ) << mesg;
    EXPECT_GT( float(count_valid)/float(disparity_map.cols()*disparity_map.rows()), attempted ) << mesg;
  }
};

typedef PyramidViewTest<PixelGray<uint8> > PyramidViewGRAYU8;
typedef PyramidViewTest<PixelGray<int16> > PyramidViewGRAYI16;
typedef PyramidViewTest<PixelGray<float> > PyramidViewGRAYF32;
typedef PyramidViewTest<uint8>             PyramidViewU8;

TEST_F( PyramidViewGRAYU8, NullPreprocess ) {
  ImageView<PixelMask<Vector2i> > disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .999, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .990, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .999, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .999, "Cross Correlation" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       CROSS_CORRELATION, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Cross Correlation" );
}

TEST_F( PyramidViewGRAYI16, NullPreprocess ) {
  ImageView<PixelMask<Vector2i> > disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .999, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .990, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .999, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .87, .99, "Cross Correlation" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       CROSS_CORRELATION, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .87, .99, "Cross Correlation" );
}

TEST_F( PyramidViewGRAYF32, NullPreprocess ) {
  ImageView<PixelMask<Vector2i> > disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .999, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .990, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .999, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .999, "Cross Correlation" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       CROSS_CORRELATION, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Cross Correlation" );
}

TEST_F( PyramidViewU8, NullPreprocess ) {
  ImageView<PixelMask<Vector2i> > disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .999, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .990, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .999, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .999, "Cross Correlation" );

  disparity_map =
    pyramid_correlate( input1, input2, NullOperation(),
                       search_volume, kernel_size,
                       CROSS_CORRELATION, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Cross Correlation" );
}
