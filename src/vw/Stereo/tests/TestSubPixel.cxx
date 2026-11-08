// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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


// TestSubPixel.h
#include <test/Helpers.h>

#include <vw/Image/Transform.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/UtilityViews.h>
#include <vw/Image/ImageChannels.h>
#include <vw/Stereo/ParabolaSubpixelView.h>
#include <vw/Stereo/PhaseSubpixelView.h>
#include <vw/Stereo/SubpixelView.h>
#include <boost/foreach.hpp>
#include <boost/random/linear_congruential.hpp>

using namespace vw;

// This setup creates two images, image1 and image2.  image1 contains a random noise
//  pattern and image2 contains a copy that has been "squished" towards the center
//  from both the left and the right sides.
template <int32 istretch>
class SubPixelCorrelateTest : public ::testing::Test {
  const int32 IMAGE_SIZE, HALF_IMAGE_SIZE;

public:
  SubPixelCorrelateTest() : IMAGE_SIZE(100), HALF_IMAGE_SIZE(50) {}

protected:

  void SetUp() {
    stretch = float32(istretch)/100;

    boost::rand48 gen(10);
    image1 = transform(channel_cast_rescale<uint8>(uniform_noise_view(gen, IMAGE_SIZE, IMAGE_SIZE)),
                       AffineTransform(Matrix2x2(3,0,0,3),Vector2()),
                       ZeroEdgeExtension(), BicubicInterpolation());
    translation = HALF_IMAGE_SIZE-HALF_IMAGE_SIZE*stretch;
    image2 = transform(image1, AffineTransform(Matrix2x2(stretch,0,0,1),
                                               Vector2(translation,0)),
                       ZeroEdgeExtension(), BicubicInterpolation());

    starting_disp.set_size(IMAGE_SIZE,IMAGE_SIZE);
    for (int32 i = 0; i < IMAGE_SIZE ; i++) {
      int32 disparity = boost::numeric_cast<int32>(stretch * i + translation - i);
      for (int32 j = 0; j < IMAGE_SIZE; j++)
        starting_disp(i,j) = Vector2f(disparity, 0);
    }
  }

  template <class ViewT>
  double check_error(ImageViewBase<ViewT> const& input,
                     int32& invalid_count) {
    ViewT const& disparity = input.impl();
    double error = 0;
    for (int32 i = 0; i < IMAGE_SIZE; i++) {
      float expected = stretch * float(i) + translation - i;
      for (int32 j = 0; j < IMAGE_SIZE; j++) {
        error += disparity(i,j)[1] + fabs(disparity(i,j)[0] - expected);
        if (!is_valid(disparity(i,j)))
          invalid_count++;
      }
    }
    return error / (double(IMAGE_SIZE)*double(IMAGE_SIZE));
  }

  float32 stretch, translation;
  ImageView<uint8> image1, image2;
  ImageView<PixelMask<Vector2f>> starting_disp;
};

typedef SubPixelCorrelateTest<95> SubPixelCorrelate95Test;
typedef SubPixelCorrelateTest<90> SubPixelCorrelate90Test;
typedef SubPixelCorrelateTest<80> SubPixelCorrelate80Test;
typedef SubPixelCorrelateTest<70> SubPixelCorrelate70Test;

// Testing Parabola Subpixel View
//--------------------------------------------------------------
TEST(ParabolaSubpixel, NullTest) {
  ImageView<PixelMask<Vector2i>> disparity(5,5);
  fill(disparity, PixelMask<Vector2i>(Vector2i(1,1)));
  ImageView<float> left(5,5), right(5,5);
  fill(left, 0.5);
  fill(right, 0.6);

  ImageView<PixelMask<Vector2f>> fdisparity =
    parabola_subpixel(pixel_cast<PixelMask<Vector2f>>(disparity),
                      pixel_cast<PixelGray<float>>(left),
                      pixel_cast<PixelGray<float>>(right),
                      PREFILTER_NONE, 1.4,
                      Vector2i(3,3));
  EXPECT_EQ(fdisparity.cols(), 5);
  EXPECT_EQ(fdisparity.rows(), 5);
  BOOST_FOREACH(PixelMask<Vector2f> const& fdisp, fdisparity) {
    EXPECT_TRUE(is_valid(fdisp));
    EXPECT_VECTOR_NEAR(fdisp.child(), Vector2f(1,1), 0.1);
  }

  fdisparity =
    parabola_subpixel(pixel_cast<PixelMask<Vector2f>>(disparity),
                      pixel_cast<PixelGray<float>>(left),
                      pixel_cast<PixelGray<float>>(right),
                      PREFILTER_LOG, 1.4,
                      Vector2i(3,3));
  EXPECT_EQ(fdisparity.cols(), 5);
  EXPECT_EQ(fdisparity.rows(), 5);
  BOOST_FOREACH(PixelMask<Vector2f> const& fdisp, fdisparity) {
    EXPECT_TRUE(is_valid(fdisp));
    EXPECT_VECTOR_NEAR(fdisp.child(), Vector2f(1,1), 0.1);
  }
}

TEST_F(SubPixelCorrelate95Test, Parabola) {
  ImageView<PixelMask<Vector2f>> disparity_map =
    parabola_subpixel(starting_disp,
                      pixel_cast<PixelGray<float>>(image1),
                      pixel_cast<PixelGray<float>>(image2),
                      PREFILTER_LOG, 1.4,
                      Vector2i(7,7));

  int32 invalid_count = 0;
  double error = check_error(disparity_map, invalid_count);
  EXPECT_LT(error, 0.6);
  EXPECT_LE(invalid_count, 0);
}

// Phase subpixel refinement on the synthetic squished pair.
TEST_F(SubPixelCorrelate80Test, Phase) {
  const BBox2 rect = bounding_box(image1);
  const int kern_size = 15; // Optimal size is products of 2's, 3's, and 5's

  ImageView<PixelMask<Vector2f>> disparity_map = crop(starting_disp, rect);
  ImageView<float> left_crop  = crop(image1,rect);
  ImageView<float> right_crop = crop(image2,rect);

  int max_pyramid_levels = 0;
  ImageView<PixelMask<Vector2f>> output_disp =
    phase_subpixel(disparity_map,
                   pixel_cast<PixelGray<float>>(left_crop),
                   pixel_cast<PixelGray<float>>(right_crop),
                   PREFILTER_LOG, 1.4,
                   Vector2i(kern_size, kern_size),
                   max_pyramid_levels);

  int32 invalid_count = 0;
  double error = check_error(output_disp, invalid_count);
  EXPECT_LT(error, 0.2);
  EXPECT_LE(invalid_count, 0);
}

// Testing Bayes EM SubPixel
//--------------------------------------------------------------
TEST_F(SubPixelCorrelate95Test, BayesEM95) {
  ImageView<PixelMask<Vector2f>> disparity_map =
    bayes_em_subpixel(starting_disp,
                      pixel_cast<PixelGray<float>>(channel_cast_rescale<float>(image1)),
                      pixel_cast<PixelGray<float>>(channel_cast_rescale<float>(image2)),
                      PREFILTER_LOG, 1.4,
                      Vector2i(7,7));

  int32 invalid_count = 0;
  double error = check_error(disparity_map, invalid_count);
  EXPECT_LT(error, 0.35);
  EXPECT_LE(invalid_count, 12);
}

TEST_F(SubPixelCorrelate90Test, BayesEM90) {
  ImageView<PixelMask<Vector2f>> disparity_map =
    bayes_em_subpixel(starting_disp,
                      pixel_cast<PixelGray<float>>(channel_cast_rescale<float>(image1)),
                      pixel_cast<PixelGray<float>>(channel_cast_rescale<float>(image2)),
                      PREFILTER_LOG, 1.4,
                      Vector2i(7,7));

  int32 invalid_count = 0;
  double error = check_error(disparity_map, invalid_count);
  EXPECT_LT(error, 0.36);
  EXPECT_LE(invalid_count, 12);
}

TEST_F(SubPixelCorrelate80Test, BayesEM80) {
  ImageView<PixelMask<Vector2f>> disparity_map =
    bayes_em_subpixel(starting_disp,
                      pixel_cast<PixelGray<float>>(channel_cast_rescale<float>(image1)),
                      pixel_cast<PixelGray<float>>(channel_cast_rescale<float>(image2)),
                      PREFILTER_LOG, 1.4,
                      Vector2i(7,7));

  int32 invalid_count = 0;
  double error = check_error(disparity_map, invalid_count);
  EXPECT_LT(error, 0.52);
  EXPECT_LE(invalid_count, 22);
}

TEST_F(SubPixelCorrelate70Test, BayesEM70) {
  ImageView<PixelMask<Vector2f>> disparity_map =
    bayes_em_subpixel(starting_disp,
                      pixel_cast<PixelGray<float>>(channel_cast_rescale<float>(image1)),
                      pixel_cast<PixelGray<float>>(channel_cast_rescale<float>(image2)),
                      PREFILTER_LOG, 1.4,
                      Vector2i(7,7));

  int32 invalid_count = 0;
  double error = check_error(disparity_map, invalid_count);
  EXPECT_LT(error, 0.9);
  EXPECT_LE(invalid_count, 48);
}
