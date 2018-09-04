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
class CorrelationViewTest : public ::testing::Test {
protected:
  CorrelationViewTest(){}

  typedef PixelT            pixel_type;
  typedef ImageView<PixelT> image_type;
  typedef typename PixelChannelType<PixelT>::type channel_type;
  image_type input1, input2;
  Vector2i   kernel_size;
  BBox2i     search_volume;

  virtual void SetUp() {
    boost::rand48 gen(10);
    input1 = ChannelRange<channel_type>::max()*uniform_noise_view( gen, 19, 25 );
    input2 = transform(input1, TranslateTransform(1,1),
                       ConstantEdgeExtension(), NearestPixelInterpolation());
    search_volume = BBox2i(1,1,1,1);
    kernel_size = Vector2i(7,7);
  }

  template <class ViewT>
  void check_error( ImageViewBase<ViewT> const& input,
                    float correct = 0.9,
                    float attempted = 0.9,
                    std::string const& mesg = "") {
    ViewT const& disparity_map = input.impl();
    int64 count_correct = 0;
    int64 count_valid = 0;
    for ( int32 i = 0; i < disparity_map.cols(); ++i )
      for ( int32 j = 0; j < disparity_map.rows(); ++j )
        if ( is_valid( disparity_map(i,j) ) ) {
          count_valid++;
          if ( disparity_map(i,j).child() == Vector2f(1,1) )
            count_correct++;
        }
    EXPECT_GT( float(count_correct)/float(count_valid), correct ) << mesg;
    EXPECT_GT( float(count_valid)/float(disparity_map.cols()*disparity_map.rows()), attempted ) << mesg;
  }
};

typedef CorrelationViewTest<PixelGray<uint8> > CorrelationViewGRAYU8;
typedef CorrelationViewTest<PixelGray<int16> > CorrelationViewGRAYI16;
typedef CorrelationViewTest<PixelGray<float> > CorrelationViewGRAYF32;
typedef CorrelationViewTest<uint8>             CorrelationViewU8;

TEST_F( CorrelationViewGRAYU8, NullPreprocess ) {
  // Percentage correct should never go below 77.4%
  ImageView<PixelMask<Vector2i> > disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .983, .99, "Absolute Difference" );

  disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               ABSOLUTE_DIFFERENCE, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .983, .99, "Absolute Difference" );


  disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .966, .99, "Squared Difference" );

  disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               SQUARED_DIFFERENCE, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .966, .99, "Squared Difference" );

  disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .966, .99, "Cross Correlation" );

  disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               CROSS_CORRELATION, 2 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .966, .99, "Cross Correlation" );
}

TEST_F( CorrelationViewGRAYI16, NullPreprocess ) {
  // Percentage correct should never go below 77.4%
  ImageView<PixelMask<Vector2i> > disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .983, .99, "Absolute Difference" );

  disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .966, .99, "Squared Difference" );

  disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .964, .99, "Cross Correlation" );
}

TEST_F( CorrelationViewGRAYF32, NullPreprocess ) {
  // Percentage correct should never go below 77.4%
  ImageView<PixelMask<Vector2i> > disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .983, .99, "Absolute Difference" );

  disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .966, .99, "Squared Difference" );

  disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .959, .99, "Cross Correlation" );
}

TEST_F( CorrelationViewU8, NullPreprocess ) {
  // Percentage correct should never go below 77.4%
  ImageView<PixelMask<Vector2i> > disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .983, .99, "Absolute Difference" );

  disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .966, .99, "Squared Difference" );

  disparity_map =
    correlate( input1, input2, NullOperation(),
               search_volume, kernel_size,
               CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .966, .99, "Cross Correlation" );
}

TEST_F( CorrelationViewGRAYU8, LaplacianOfGaussian) {
  // Percentage correct should never go below 77.4%
  ImageView<PixelMask<Vector2i> > disparity_map =
    correlate( input1, input2, LaplacianOfGaussian(1.4),
               search_volume, kernel_size,
               ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .81, .99, "Absolute Difference" );

  disparity_map =
    correlate( input1, input2, LaplacianOfGaussian(1.4),
               search_volume, kernel_size,
               SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .8, .99, "Squared Difference" );

  disparity_map =
    correlate( input1, input2, LaplacianOfGaussian(1.4),
               search_volume, kernel_size,
               CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .8, .99, "Cross Correlation" );
}

TEST_F( CorrelationViewGRAYU8, SubtractedMean ) {
  // Percentage correct should never go below 77.4%
  ImageView<PixelMask<Vector2i> > disparity_map =
    correlate( input1, input2, SubtractedMean(5),
               search_volume, kernel_size,
               ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .93, .99, "Absolute Difference" );

  disparity_map =
    correlate( input1, input2, SubtractedMean(5),
               search_volume, kernel_size,
               SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .93, .99, "Squared Difference" );

  disparity_map =
    correlate( input1, input2, SubtractedMean(5),
               search_volume, kernel_size,
               CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .92, .99, "Cross Correlation" );
}

TEST_F( CorrelationViewGRAYF32, LaplacianOfGaussian) {
  ImageView<PixelMask<Vector2i> > disparity_map =
    correlate( input1, input2, LaplacianOfGaussian(1.4),
               search_volume, kernel_size,
               ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .7, .99, "Absolute Difference" );

  disparity_map =
    correlate( input1, input2, LaplacianOfGaussian(1.4),
               search_volume, kernel_size,
               SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .64, .99, "Squared Difference" );

  disparity_map =
    correlate( input1, input2, LaplacianOfGaussian(1.4),
               search_volume, kernel_size,
               CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .66, .99, "Cross Correlation" );
}

TEST_F( CorrelationViewGRAYF32, SubtractedMean ) {
  // Percentage correct should never go below 77.4%
  ImageView<PixelMask<Vector2i> > disparity_map =
    correlate( input1, input2, SubtractedMean(5),
               search_volume, kernel_size,
               ABSOLUTE_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .974, .99, "Absolute Difference" );

  disparity_map =
    correlate( input1, input2, SubtractedMean(5),
               search_volume, kernel_size,
               SQUARED_DIFFERENCE, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .966, .99, "Squared Difference" );

  disparity_map =
    correlate( input1, input2, SubtractedMean(5),
               search_volume, kernel_size,
               CROSS_CORRELATION, -1 );
  ASSERT_EQ( input1.cols(), disparity_map.cols() );
  ASSERT_EQ( input1.rows(), disparity_map.rows() );
  check_error( disparity_map, .966, .99, "Cross Correlation" );
}
