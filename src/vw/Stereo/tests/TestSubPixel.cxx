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


// TestCorrelator.h
#include <test/Helpers.h>

#include <vw/Image.h>
#include <vw/Stereo/SubpixelView.h>
#include <vw/FileIO.h>
#include <boost/foreach.hpp>
#include <boost/random/linear_congruential.hpp>

using namespace vw;
using namespace vw::stereo;

namespace vw {
  template<> struct PixelFormatID<Vector3>   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_3_CHANNEL; };
}

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
      int32 disparity = boost::numeric_cast<int32>(stretch * i + translation - i);
      for ( int32 j = 0; j < IMAGE_SIZE; j++ )
        starting_disp(i,j) = Vector2f(disparity, 0);
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

typedef SubPixelCorrelateTest<95> SubPixelCorrelate95Test;
typedef SubPixelCorrelateTest<90> SubPixelCorrelate90Test;
typedef SubPixelCorrelateTest<80> SubPixelCorrelate80Test;
typedef SubPixelCorrelateTest<70> SubPixelCorrelate70Test;

// Testing Parabola Subpixel View
//--------------------------------------------------------------
TEST( ParabolaSubpixel, NullTest ) {
  ImageView<PixelMask<Vector2i> > disparity(5,5);
  fill( disparity, PixelMask<Vector2i>(Vector2i(1,1)) );
  ImageView<float> left(5,5), right(5,5);
  fill( left, 0.5 );
  fill( right, 0.6 );

  ImageView<PixelMask<Vector2f> > fdisparity =
    parabola_subpixel( disparity, left, right,
                       PREFILTER_NONE, 1.4,
                       Vector2i(3,3) );
  EXPECT_EQ( fdisparity.cols(), 5 );
  EXPECT_EQ( fdisparity.rows(), 5 );
  BOOST_FOREACH( PixelMask<Vector2f> const& fdisp, fdisparity ) {
    EXPECT_TRUE( is_valid( fdisp ) );
    EXPECT_VECTOR_NEAR( fdisp.child(), Vector2f(1,1), 0.1 );
  }

  fdisparity =
    parabola_subpixel( disparity, left, right,
                       PREFILTER_LOG, 1.4,
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
                       PREFILTER_LOG, 1.4,
                       Vector2i(7,7) );

  int32 invalid_count = 0;
  double error = check_error( disparity_map, invalid_count );
  EXPECT_LT(error, 0.6);
  EXPECT_LE(invalid_count, 0);
}

// Testing Bayes EM SubPixel
//--------------------------------------------------------------
TEST_F( SubPixelCorrelate95Test, BayesEM95 ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    bayes_em_subpixel( starting_disp,
                       channel_cast_rescale<float>(image1),
                       channel_cast_rescale<float>(image2),
                       PREFILTER_LOG, 1.4,
                       Vector2i(7,7) );
  int32 invalid_count = 0;
  double error = check_error( disparity_map, invalid_count );
  //std::cout << "Err: " << error << " Cnt: " << invalid_count << "\n";
  EXPECT_LT(error, 0.35);
  EXPECT_LE(invalid_count, 12);
}

TEST_F( SubPixelCorrelate90Test, BayesEM90 ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    bayes_em_subpixel( starting_disp,
                       channel_cast_rescale<float>(image1),
                       channel_cast_rescale<float>(image2),
                       PREFILTER_LOG, 1.4,
                       Vector2i(7,7) );
  int32 invalid_count = 0;
  double error = check_error( disparity_map, invalid_count );
  //std::cout << "Err: " << error << " Cnt: " << invalid_count << "\n";
  EXPECT_LT(error, 0.36);
  EXPECT_LE(invalid_count, 12);
}

TEST_F( SubPixelCorrelate80Test, BayesEM80 ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    bayes_em_subpixel( starting_disp,
                       channel_cast_rescale<float>(image1),
                       channel_cast_rescale<float>(image2),
                       PREFILTER_LOG, 1.4,
                       Vector2i(7,7) );
  int32 invalid_count = 0;
  double error = check_error( disparity_map, invalid_count );
  //std::cout << "Err: " << error << " Cnt: " << invalid_count << "\n";
  EXPECT_LT(error, 0.52);
  EXPECT_LE(invalid_count, 22);
}

TEST_F( SubPixelCorrelate70Test, BayesEM70 ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    bayes_em_subpixel( starting_disp,
                       channel_cast_rescale<float>(image1),
                       channel_cast_rescale<float>(image2),
                       PREFILTER_LOG, 1.4,
                       Vector2i(7,7) );
  int32 invalid_count = 0;
  double error = check_error( disparity_map, invalid_count );
  //std::cout << "Err: " << error << " Cnt: " << invalid_count << "\n";
  EXPECT_LT(error, 0.9);
  EXPECT_LE(invalid_count, 48);
}
