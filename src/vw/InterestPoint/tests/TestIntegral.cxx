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


// TestIntegral.cxx
#include <gtest/gtest_VW.h>

#include <vw/InterestPoint/IntegralImage.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Interpolation.h>
#include <vw/FileIO/DiskImageResource.h>

using namespace vw;
using namespace vw::ip;

// TODO(oalexan1): Why tests with png fail?
#if 0
TEST( Integral, IntegralSumming ) {
  ImageView<float> graffiti;
  read_image( graffiti, TEST_SRCDIR"/sub.png" );
  ImageView<float> integral;
  integral = IntegralImage( graffiti );

  for (unsigned size = 10; size <= 70; size+= 10) {
    float sum = IntegralBlock( integral,
                               Vector2i(0,0),
                               Vector2i(size,size) );
    float actual_sum = 0;
    for (unsigned i = 0; i < size; i++ ) {
      for (unsigned j = 0; j < size; j++ ) {
        actual_sum+=graffiti(i,j);
      }
    }

    EXPECT_NEAR(actual_sum, sum, 1e-2);
  }
}

TEST( Integral, HaarFilters ) {
  // Loading Image
  ImageView<float> graffiti, gradient;
  read_image( graffiti, TEST_SRCDIR"/sub.png" );
  read_image( gradient, TEST_SRCDIR"/noisy_gradient_60.png" );

  // Building Integrals
  ImageView<float> integral = IntegralImage( graffiti );

  float hand_response = 0;
  for ( unsigned i = 0; i < 2; i++ ) {
    for ( unsigned j = 0; j < 2; j++ ) {
      if ( i == 0 ) {
        hand_response -= graffiti(i,j);
      } else {
        hand_response += graffiti(i,j);
      }
    }
  }
  EXPECT_NEAR( hand_response,
               HHaarWavelet( integral,
                             1,1,2 ),
               1e-4 );

  hand_response = 0;
  for ( unsigned i = 0; i < 2; i++ ) {
    for ( unsigned j = 0; j < 2; j++ ) {
      if ( j == 0 ) {
        hand_response -= graffiti(i,j);
      } else {
        hand_response += graffiti(i,j);
      }
    }
  }
  EXPECT_NEAR( hand_response,
               VHaarWavelet( integral,
                             1, 1, 2 ),
               1e-4 );
}

TEST( Integral, HaarFilters2 ) {
  // Load Image
  ImageView<float> graffiti, gradient;
  read_image( graffiti, TEST_SRCDIR"/sub.png" );
  read_image( gradient, TEST_SRCDIR"/noisy_gradient_60.png" );

  // Rotating Image
  ImageView<float> graffiti_r, gradient_r;
  graffiti_r = rotate_180(graffiti);
  gradient_r = rotate_180(gradient);

  // Building Integrals
  ImageView<float> int_graf, int_graf_r, int_grad, int_grad_r;
  int_graf = IntegralImage( graffiti );
  int_graf_r = IntegralImage( graffiti_r);
  int_grad = IntegralImage( gradient );
  int_grad_r = IntegralImage( gradient_r );

  for ( unsigned i = 10; i < 90; i+=5 ) {
    for ( unsigned j = 10; j < 50; j+=5 ) {
      // Size 10
      EXPECT_NEAR( -HHaarWavelet( int_graf, i,j, 10),
                   HHaarWavelet( int_graf_r, int_graf.cols()-1-i,
                                 int_graf.rows()-1-j, 10 ),
                   1e-4 );
      EXPECT_NEAR( -VHaarWavelet( int_graf, i,j, 10),
                   VHaarWavelet( int_graf_r, int_graf.cols()-1-i,
                                 int_graf.rows()-1-j, 10 ),
                   1e-4 );
      EXPECT_NEAR( -HHaarWavelet( int_grad, i,j, 10),
                   HHaarWavelet( int_grad_r, int_grad.cols()-1-i,
                                 int_grad.rows()-1-j, 10 ),
                   1e-4 );
      EXPECT_NEAR( -VHaarWavelet( int_grad, i,j, 10),
                   VHaarWavelet( int_grad_r, int_grad.cols()-1-i,
                                 int_grad.rows()-1-j, 10 ),
                   1e-4 );

      // Size 20
      EXPECT_NEAR( -HHaarWavelet( int_graf, i,j, 20),
                   HHaarWavelet( int_graf_r, int_graf.cols()-1-i,
                                 int_graf.rows()-1-j, 20 ),
                   1e-4 );
      EXPECT_NEAR( -VHaarWavelet( int_graf, i,j, 20),
                   VHaarWavelet( int_graf_r, int_graf.cols()-1-i,
                                 int_graf.rows()-1-j, 20 ),
                   1e-4 );
      EXPECT_NEAR( -HHaarWavelet( int_grad, i,j, 20),
                   HHaarWavelet( int_grad_r, int_grad.cols()-1-i,
                                 int_grad.rows()-1-j, 20 ),
                   1e-4 );
      EXPECT_NEAR( -VHaarWavelet( int_grad, i,j, 20),
                   VHaarWavelet( int_grad_r, int_grad.cols()-1-i,
                                 int_grad.rows()-1-j, 20 ),
                   1e-4 );

      // Size 4
      EXPECT_NEAR( -HHaarWavelet( int_graf, i,j, 4),
                   HHaarWavelet( int_graf_r, int_graf.cols()-1-i,
                                 int_graf.rows()-1-j, 4 ),
                   1e-4 );
      EXPECT_NEAR( -VHaarWavelet( int_graf, i,j, 4),
                   VHaarWavelet( int_graf_r, int_graf.cols()-1-i,
                                 int_graf.rows()-1-j, 4 ),
                   1e-4 );
      EXPECT_NEAR( -HHaarWavelet( int_grad, i,j, 4),
                   HHaarWavelet( int_grad_r, int_grad.cols()-1-i,
                                 int_grad.rows()-1-j, 4 ),
                   1e-4 );
      EXPECT_NEAR( -VHaarWavelet( int_grad, i,j, 4),
                   VHaarWavelet( int_grad_r, int_grad.cols()-1-i,
                                 int_grad.rows()-1-j, 4 ),
                   1e-4 );

      // Size 16
      EXPECT_NEAR( -HHaarWavelet( int_graf, i,j, 16),
                   HHaarWavelet( int_graf_r, int_graf.cols()-1-i,
                                 int_graf.rows()-1-j, 16 ),
                   1e-4 );
      EXPECT_NEAR( -VHaarWavelet( int_graf, i,j, 16),
                   VHaarWavelet( int_graf_r, int_graf.cols()-1-i,
                                 int_graf.rows()-1-j, 16 ),
                   1e-4 );
      EXPECT_NEAR( -HHaarWavelet( int_grad, i,j, 16),
                   HHaarWavelet( int_grad_r, int_grad.cols()-1-i,
                                 int_grad.rows()-1-j, 16 ),
                   1e-4 );
      EXPECT_NEAR( -VHaarWavelet( int_grad, i,j, 16),
                   VHaarWavelet( int_grad_r, int_grad.cols()-1-i,
                                 int_grad.rows()-1-j, 16 ),
                   1e-4 );
    }
  }
}

TEST( Integral, DerivativeFilters ) {
  ImageView<float> graffiti;
  read_image( graffiti, TEST_SRCDIR"/sub.png" );
  ImageView<float> integral;
  integral = IntegralImage( graffiti );
  ImageView<float> rotated;
  rotated = rotate_180(graffiti);
  ImageView<float> r_integral;
  r_integral = IntegralImage( rotated );

  EXPECT_NEAR( XSecondDerivative( integral,
                                  49, 49, 51 ),
               XSecondDerivative( r_integral,
                                  50, 50, 51 ),
               1e-4 );

  EXPECT_NEAR( YSecondDerivative( integral,
                                  49, 49, 51 ),
               YSecondDerivative( r_integral,
                                      50, 50, 51 ),
               1e-4 );

  EXPECT_NEAR( XYDerivative( integral,
                             49, 49, 51 ),
               XYDerivative( r_integral,
                             50, 50, 51 ),
               1e-4 );

}

TEST( Integral, InterpolationProof ) {
  // Load Image
  ImageView<float> gradient;
  read_image( gradient, TEST_SRCDIR"/noisy_gradient_60.png" );

  // Building Integrals
  ImageView<float> integral;
  integral = IntegralImage( gradient );

  // Interpolation by hand from 2 filters
  float left, right;
  left = HHaarWavelet(  integral, int(10), int(10), 10 );
  right = HHaarWavelet( integral, int(11), int(10), 10 );

  // Wrapping Integral
  InterpolationView<EdgeExtensionView<ImageView<float>, ConstantEdgeExtension>, BilinearInterpolation> wrapped_integral = interpolate( integral, BilinearInterpolation() );

  // Calculating interpolated responses
  float int_response;
  int_response = -wrapped_integral(10.5-5,10-5);
  int_response += 2*wrapped_integral(10.5,10-5);
  int_response -= wrapped_integral(10.5+5,10-5);
  int_response += wrapped_integral(10.5-5, 10+5);
  int_response -= 2*wrapped_integral(10.5, 10+5);
  int_response += wrapped_integral(10.5+5, 10+5);

  EXPECT_NEAR( (left+right)/2, int_response, 1e-4 );

  // Using floating point function
  EXPECT_NEAR( (left+right)/2,
               HHaarWavelet( wrapped_integral,
                             10.5, 10.0, 10 ),
               1e-4 );
}

#endif
