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


#include <test/Helpers.h>

#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/CorrelationView.h>
#include <boost/random/linear_congruential.hpp>


using namespace vw;
using namespace vw::stereo;

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
  int corr_timeout;
  double seconds_per_op;
  int max_levels;
  int filter_radius;
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
    corr_timeout = 0;
    seconds_per_op = 0;
    max_levels = 5;
    filter_radius = 5;

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
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .909, .9985, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .990, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .9985, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       CROSS_CORRELATION,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .998, "Cross Correlation" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       CROSS_CORRELATION,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Cross Correlation" );
}

TEST_F( PyramidViewGRAYI16, NullPreprocess ) {
  ImageView<PixelMask<Vector2i> > disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .909, .9985, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .990, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .9985, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       CROSS_CORRELATION,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .87, .99, "Cross Correlation" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       CROSS_CORRELATION,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .87, .99, "Cross Correlation" );
}

TEST_F( PyramidViewGRAYF32, NullPreprocess ) {
  ImageView<PixelMask<Vector2i> > disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .909, .9985, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .990, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .9985, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       CROSS_CORRELATION,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .998, "Cross Correlation" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       CROSS_CORRELATION,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Cross Correlation" );
}

TEST_F( PyramidViewU8, NullPreprocess ) {
  ImageView<PixelMask<Vector2i> > disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .909, .9985, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       ABSOLUTE_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .91, .990, "Absolute Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .9985, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       SQUARED_DIFFERENCE,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Squared Difference" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       CROSS_CORRELATION,
                       corr_timeout, seconds_per_op,
                       -1, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .998, "Cross Correlation" );

  disparity_map =
    pyramid_correlate( input1, input2,
                       constant_view(uint8(255), input1),
                       constant_view(uint8(255), input2),
                       PREFILTER_NONE, 0,
                       search_volume, kernel_size,
                       CROSS_CORRELATION,
                       corr_timeout, seconds_per_op,
                       2, 0, filter_radius, max_levels );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .90, .990, "Cross Correlation" );
}
