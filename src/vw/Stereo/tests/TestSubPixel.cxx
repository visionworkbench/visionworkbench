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


// TestSubPixel.h
#include <test/Helpers.h>

#include <vw/Image/Transform.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/PixelMask.h>
#include <vw/Stereo/ParabolaSubpixelView.h>
#include <vw/Stereo/PhaseSubpixelView.h>
#include <vw/Stereo/SubpixelView.h>
#include <boost/foreach.hpp>
#include <boost/random/linear_congruential.hpp>

using namespace vw;
using namespace vw::stereo;

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

//typedef SubPixelCorrelateTest<95> SubPixelCorrelate95Test;
//typedef SubPixelCorrelateTest<90> SubPixelCorrelate90Test;
//typedef SubPixelCorrelateTest<80> SubPixelCorrelate80Test;
//typedef SubPixelCorrelateTest<70> SubPixelCorrelate70Test;

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


/*
/// Test the low level phase correlation code.
TEST_F( SubPixelCorrelate95Test, phase_simple) {

  Vector2f offset;
  BBox2 rect(25, 46, 15, 15);
  int subpixel_accuracy=50;
  ImageView<float> left_crop  = crop(image1,rect);
  ImageView<float> right_crop = crop(image2,rect);
  vw::stereo::phase_correlation_subpixel(left_crop, right_crop,
                                         offset, subpixel_accuracy, true); 
  
  EXPECT_NEAR(-100, offset[0], 0.01);
  EXPECT_NEAR(-100, offset[1], 0.01);
  EXPECT_TRUE(false);

}
*/

TEST_F( SubPixelCorrelate80Test, Phase) {
  
  //const BBox2 rect(2, 20, 90, 50);
  const BBox2 rect = bounding_box(image1);
  const int kern_size = 15; // Optimal size is products of 2's, 3's, and 5's
  const int subpixel_accuracy = 50;
  
  ImageView<PixelMask<Vector2f> > disparity_map = crop(starting_disp, rect);
  ImageView<float> left_crop  = crop(image1,rect);
  ImageView<float> right_crop = crop(image2,rect);
  /*
  ImageView<PixelMask<Vector2f> > disparity_map(rect.width(), rect.height());
  for(int r=0; r<rect.height(); ++r) {
    for(int c=0; c<rect.width(); ++c) {
      remove_mask(disparity_map(c,r)) = Vector2f(0,0);
      validate(disparity_map(c,r));
    }
  }
  */
  //write_image("/home/smcmich1/data/subpixel/crop1.tif", left_crop);
  //write_image("/home/smcmich1/data/subpixel/crop2.tif", right_crop);

  //write_image("/home/smcmich1/data/subpixel/starting_disp_x.tif",
  //            select_channel(disparity_map,0));
  //write_image("/home/smcmich1/data/subpixel/starting_disp_y.tif",
  //            select_channel(disparity_map,1));
  
  /*
  // The basic, single resolution call.
  subpixel_phase_2d(disparity_map,
                    left_crop, right_crop,
                    kern_size, kern_size,
                    bounding_box(left_crop),
                    subpixel_accuracy);
  write_image("/home/smcmich1/data/subpixel/disparity_map_out_x.tif",
              select_channel(disparity_map,0));
  */
  std::cout << "Starting test!\n";
  int max_pyramid_levels = 0;
  ImageView<PixelMask<Vector2f> > output_disp =
  phase_subpixel(disparity_map,
                 left_crop, right_crop,
                 PREFILTER_LOG, 1.4,
                 Vector2i(kern_size, kern_size),
                 max_pyramid_levels);
  
  //write_image("/home/smcmich1/data/subpixel/disparity_map_out_pyr_x.tif",
  //            select_channel(output_disp,0));
  
  
  int32 invalid_count = 0;
  //double error = check_error( disparity_map, invalid_count );
  double error = check_error( output_disp, invalid_count );
  //std::cout << "Err: " << error << " Cnt: " << invalid_count << "\n";
  EXPECT_LT(error, 0.2);
  EXPECT_LE(invalid_count, 0);
  
  //EXPECT_TRUE(false);
}





/*
TEST_F( SubPixelCorrelate90Test, affine_pyr) {
  
  //const BBox2 rect(2, 20, 90, 50);
  const BBox2 rect = bounding_box(image1);
  const int kern_size = 15;

  ImageView<float> left_crop  = crop(image1,rect);
  ImageView<float> right_crop = crop(image2,rect);

  ImageView<PixelMask<Vector2f> > disp_copy = crop(starting_disp, rect);
  ImageView<PixelMask<Vector2f> > disparity_map =
    affine_subpixel( disp_copy, left_crop, right_crop,
                       PREFILTER_LOG, 1.4,
                       Vector2i(kern_size,kern_size) );
  
  //write_image("/home/smcmich1/data/subpixel/affine90_pyr_out_x.tif",
  //            select_channel(disparity_map,0));
  
  int32 invalid_count = 0;
  double error = check_error( disparity_map, invalid_count );
  //std::cout << "Err: " << error << " Cnt: " << invalid_count << "\n";
  EXPECT_LT(error, 0.35);
  EXPECT_LE(invalid_count, 12);
  //EXPECT_TRUE(false);
}
*/
/*
TEST_F( SubPixelCorrelate80Test, affine_pyr) {
  
  //const BBox2 rect(2, 20, 90, 50);
  const BBox2 rect = bounding_box(image1);
  const int kern_size = 15;

  ImageView<float> left_crop  = crop(image1,rect);
  ImageView<float> right_crop = crop(image2,rect);

  ImageView<PixelMask<Vector2f> > disp_copy = crop(starting_disp, rect);
  ImageView<PixelMask<Vector2f> > disparity_map =
    affine_subpixel( disp_copy, left_crop, right_crop,
                       PREFILTER_LOG, 1.4,
                       Vector2i(kern_size,kern_size) );
  
  //write_image("/home/smcmich1/data/subpixel/affine80_pyr_out_x.tif",
  //            select_channel(disparity_map,0));

  int32 invalid_count = 0;
  double error = check_error( disparity_map, invalid_count );
  std::cout << "Err: " << error << " Cnt: " << invalid_count << "\n";
  EXPECT_LT(error, 0.35);
  EXPECT_LE(invalid_count, 12);
  EXPECT_TRUE(false);
}
*/
/*
TEST_F( SubPixelCorrelate95Test, parab) {
  
  const BBox2 rect(2, 20, 90, 50);
  //const BBox2 rect = bounding_box(image1);
  const int kern_size = 9;

  ImageView<float> left_crop  = crop(image1,rect);
  ImageView<float> right_crop = crop(image2,rect);

  ImageView<PixelMask<Vector2f> > disp_copy = crop(starting_disp, rect);
  ImageView<PixelMask<Vector2f> > disparity_map =
    parabola_subpixel( disp_copy, left_crop, right_crop,
                       PREFILTER_LOG, 1.4,
                       Vector2i(kern_size,kern_size) );
  
  write_image("/home/smcmich1/data/subpixel/parabola_disparity_map_out_x.tif",
              select_channel(disparity_map,0));
  
  EXPECT_TRUE(false);
}
*/
/*
TEST_F( SubPixelCorrelate80Test, Affine) {
  
  //const BBox2 rect(2, 20, 90, 50);
  const BBox2 rect = bounding_box(image1);
  const int kern_size = 15;

  ImageView<PixelMask<Vector2f> > disparity_map = crop(starting_disp, rect);
  //ImageView<PixelMask<Vector2f> > disparity_map = starting_disp;
  ImageView<float> left_crop  = crop(image1,rect);
  ImageView<float> right_crop = crop(image2,rect);

  subpixel_optimized_affine_2d(disparity_map,
                    left_crop, right_crop,
                    kern_size, kern_size,
                    bounding_box(left_crop), true, true, false);
  
  //write_image("/home/smcmich1/data/subpixel/affine80_out_x.tif",
  //            select_channel(disparity_map,0));
  
  EXPECT_TRUE(false);
}
*/

/*

// TODO: Add more phase tests
TEST_F( SubPixelCorrelate95Test, Phase ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    phase_subpixel( starting_disp, image1, image2,
                       PREFILTER_LOG, 1.4,
                       Vector2i(7,7) );

  int32 invalid_count = 0;
  double error = check_error( disparity_map, invalid_count );
  EXPECT_LT(error, 0.6);
  EXPECT_LE(invalid_count, 0);
}
*/

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
  //EXPECT_TRUE(false);
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
  //EXPECT_TRUE(false);
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
  //EXPECT_TRUE(false);
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
  //EXPECT_TRUE(false);
}
/*
TEST( SlipTest, Affine ) {
  
  ImageView<float> image1, image2;
  read_image(image1, "/home/smcmich1/data/subpixel/landsat/l_crop.tif");
  read_image(image2, "/home/smcmich1/data/subpixel/landsat/r_crop.tif");
  
  // Init the disparity image.
  ImageView<PixelMask<Vector2f> > starting_disp(image1.cols(), image1.rows());
  for(int r=0; r<starting_disp.rows(); ++r) {
    for(int c=0; c<starting_disp.cols(); ++c) {
      remove_mask(starting_disp(c,r)) = Vector2f(0,0);
      validate(starting_disp(c,r));
    }
  }

  ImageView<PixelMask<Vector2f> > disparity_map =
      affine_subpixel( starting_disp,
                       image1,
                       image2,
                       PREFILTER_LOG, 1.4, // Having a prefilter is important!
                       Vector2i(35,35) );

  //ImageView<PixelMask<Vector2f> > disparity_map = starting_disp;
  //subpixel_optimized_affine_2d(disparity_map,
  //                  image1, image2,
  //                  35, 35,
  //                  bounding_box(image1), true, true, false);

  write_image("/home/smcmich1/data/subpixel/landsat/slip_affine.tif",
              select_channel(disparity_map,0));
  EXPECT_TRUE(false);
}

*/
/*
TEST( SlipTest, EM2) {

  const int kern_size = 25;
  
  ImageView<float> image1, image2;
  read_image(image1, "/home/smcmich1/data/subpixel/landsat/l_crop.tif");
  read_image(image2, "/home/smcmich1/data/subpixel/landsat/r_crop.tif");
  
  // Init the disparity image.
  ImageView<PixelMask<Vector2f> > disparity_map(image1.cols(), image1.rows());
  for(int r=0; r<disparity_map.rows(); ++r) {
    for(int c=0; c<disparity_map.cols(); ++c) {
      remove_mask(disparity_map(c,r)) = Vector2f(0,0);
      validate(disparity_map(c,r));
    }
  }

  // The basic, single resolution call.
  subpixel_optimized_affine_2d_EM(disparity_map,
                    image1, image2,
                    kern_size, kern_size,
                    bounding_box(image1),
                    true, true, false);
  write_image("/home/smcmich1/data/subpixel/landsat/slip_EM2.tif",
              select_channel(disparity_map,0));

  
  EXPECT_TRUE(false);
}
*/
/*
TEST( SlipTest, Phase) {

  const int kern_size = 35; // Optimal size is products of 2's, 3's, and 5's
  const int subpixel_accuracy = 20;
  
  ImageView<float> image1, image2;
  read_image(image1, "/home/smcmich1/data/subpixel/landsat/l_crop.tif");
  read_image(image2, "/home/smcmich1/data/subpixel/landsat/r_crop.tif");
  
  // Init the disparity image.
  ImageView<PixelMask<Vector2f> > starting_disp(image1.cols(), image1.rows());
  for(int r=0; r<starting_disp.rows(); ++r) {
    for(int c=0; c<starting_disp.cols(); ++c) {
      remove_mask(starting_disp(c,r)) = Vector2f(0,0);
      validate(starting_disp(c,r));
    }
  }

  //// The basic, single resolution call.
  //subpixel_phase_2d(disparity_map,
  //                  image1, image2,
  //                  kern_size, kern_size,
  //                  bounding_box(image1),
  //                  subpixel_accuracy);
  //write_image("/home/smcmich1/data/subpixel/landsat/slip_phase.tif",
  //            select_channel(disparity_map,0));

  
  ImageView<PixelMask<Vector2f> > disparity_map =
      phase_subpixel( starting_disp,
                       image1,
                       image2,
                       PREFILTER_LOG, 1.4, // Having a prefilter is important!
                       Vector2i(35,35) );
      
  write_image("/home/smcmich1/data/subpixel/landsat/slip_phase_pyr.tif",
              select_channel(disparity_map,0));

  EXPECT_TRUE(false);
}
*/
