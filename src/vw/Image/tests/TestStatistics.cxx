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


// Image/tests/TestStatistics.h
#include <gtest/gtest_VW.h>

#include <vw/Image/Statistics.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>

#include <test/Helpers.h>

using namespace vw;

ImageView<vw::uint8> im8(2,1), im8b(3,2);
ImageView<double> imf(2,1);
ImageView<PixelMask<vw::uint8> > im8ma(2,1);
ImageView<PixelMask<double> > imfma(2,1);
ImageView<PixelMask<vw::uint8> > im8mb(2,1);
ImageView<PixelMask<double> > imfmb(2,1);
ImageView<PixelRGB<vw::uint8> > imrgb8(2,1);
ImageView<PixelRGB<double> > imrgbf(2,1);
ImageView<PixelRGBA<vw::uint8> > imrgba8a(2,1);
ImageView<PixelRGBA<double> > imrgbafa(2,1);
ImageView<PixelRGBA<vw::uint8> > imrgba8b(2,1);
ImageView<PixelRGBA<double> > imrgbafb(2,1);
ImageView<PixelRGBA<vw::uint8> > imrgba8c(2,1);
ImageView<PixelRGBA<double> > imrgbafc(2,1);

void generate_data( void ) {
  im8(0,0) = 24;
  im8(1,0) = 119;

  im8b(0,0) = 24;
  im8b(0,1) = 255;
  im8b(1,0) = 119;
  im8b(1,1) = 119;
  im8b(2,0) = 0;
  im8b(2,1) = 0;

  imf(0,0) = 1.0;
  imf(1,0) = -1.0;

  im8ma(0,0) = 0;
  im8ma(1,0) = 119;

  imfma(0,0) = 0;
  imfma(1,0) = -1.0;

  imrgb8(0,0) = PixelRGB<vw::uint8>(24,67,89);
  imrgb8(1,0) = PixelRGB<vw::uint8>(119,76,228);

  imrgbf(0,0) = PixelRGB<double>(0.5,0.9,1.6);
  imrgbf(1,0) = PixelRGB<double>(-1,0,1);

  imrgba8a(0,0) = PixelRGBA<vw::uint8>(24,67,89,255);
  imrgba8a(1,0) = PixelRGBA<vw::uint8>(24,15,46,51);

  imrgbafa(0,0) = PixelRGBA<double>(0.5,0.9,1.6,1.0);
  imrgbafa(1,0) = PixelRGBA<double>(-0.5,0,0.5,0.5);

  imrgba8b(0,0) = PixelRGBA<vw::uint8>(24,67,89,255);
  imrgba8b(1,0) = PixelRGBA<vw::uint8>(20,15,100,0);

  imrgbafb(0,0) = PixelRGBA<double>(0.5,0.9,1.6,1.0);
  imrgbafb(1,0) = PixelRGBA<double>(-0.6,0,2.0,0);

  imrgba8c(0,0) = PixelRGBA<vw::uint8>(24,67,89,0);
  imrgba8c(1,0) = PixelRGBA<vw::uint8>(119,76,228,0);

  imrgbafc(0,0) = PixelRGBA<double>(0.5,0.9,1.6,0);
  imrgbafc(1,0) = PixelRGBA<double>(-1,0,1,0);
}

template <class T1, class T2>
static bool is_of_type( T2 ) {
  return boost::is_same<T1,T2>::value;
}

TEST( Statistics, MinPixel ) {
  generate_data();
  EXPECT_EQ( min_pixel_value(channels_to_planes(imrgb8)), 24 );
  ASSERT_TRUE( is_of_type<vw::uint8>( min_pixel_value(channels_to_planes(imrgb8)) ) );
  EXPECT_EQ( min_pixel_value(channels_to_planes(imrgbf)), -1.0 );
  ASSERT_TRUE( is_of_type<double>( min_pixel_value(channels_to_planes(imrgbf)) ) );
  EXPECT_EQ( min_pixel_value(im8), 24 );
  ASSERT_TRUE( is_of_type<vw::uint8>( min_pixel_value(im8) ) );
  EXPECT_EQ( min_pixel_value(imf), -1.0 );
  ASSERT_TRUE( is_of_type<double>( min_pixel_value(imf) ) );
  EXPECT_EQ( min_pixel_value(im8ma), 0 );
  ASSERT_TRUE( is_of_type<vw::uint8>( min_pixel_value(im8ma) ) );
  EXPECT_EQ( min_pixel_value(imfma), -1.0 );
  ASSERT_TRUE( is_of_type<double>( min_pixel_value(imfma) ) );
  ASSERT_THROW( min_pixel_value(im8mb), ArgumentErr );
  ASSERT_THROW( min_pixel_value(imfmb), ArgumentErr );
}

TEST( Statistics, MaxPixel ) {
  generate_data();
  EXPECT_EQ( max_pixel_value(channels_to_planes(imrgb8)), 228 );
  ASSERT_TRUE( is_of_type<vw::uint8>( max_pixel_value(channels_to_planes(imrgb8)) ) );
  EXPECT_EQ( max_pixel_value(channels_to_planes(imrgbf)), 1.6 );
  ASSERT_TRUE( is_of_type<double>( max_pixel_value(channels_to_planes(imrgbf)) ) );
  EXPECT_EQ( max_pixel_value(im8), 119 );
  ASSERT_TRUE( is_of_type<vw::uint8>( max_pixel_value(im8) ) );
  EXPECT_EQ( max_pixel_value(imf), 1.0 );
  ASSERT_TRUE( is_of_type<double>( max_pixel_value(imf) ) );
  EXPECT_EQ( max_pixel_value(im8ma), 119 );
  ASSERT_TRUE( is_of_type<vw::uint8>( max_pixel_value(im8ma) ) );
  EXPECT_EQ( max_pixel_value(imfma), 0.0 );
  ASSERT_TRUE( is_of_type<double>( max_pixel_value(imfma) ) );
  ASSERT_THROW( max_pixel_value(im8mb), ArgumentErr );
  ASSERT_THROW( max_pixel_value(imfmb), ArgumentErr );
}

TEST( Statistics, MinMaxPixel ) {
  generate_data();
  vw::uint8 mini, maxi;
  double minf, maxf;
  min_max_pixel_values( channels_to_planes(imrgb8), mini, maxi );
  EXPECT_EQ( mini, 24 );
  EXPECT_EQ( maxi, 228 );
  min_max_pixel_values( channels_to_planes(imrgbf), minf, maxf );
  EXPECT_EQ( minf, -1.0 );
  EXPECT_EQ( maxf, 1.6 );
  min_max_pixel_values( im8, mini, maxi );
  EXPECT_EQ( mini, 24 );
  EXPECT_EQ( maxi, 119 );
  min_max_pixel_values( imf, minf, maxf );
  EXPECT_EQ( minf, -1.0 );
  EXPECT_EQ( maxf, 1.0 );
  min_max_pixel_values( im8ma, mini, maxi );
  EXPECT_EQ( mini, 0 );
  EXPECT_EQ( maxi, 119 );
  min_max_pixel_values( imfma, minf, maxf );
  EXPECT_EQ( minf, -1.0 );
  EXPECT_EQ( maxf, 0 );
  ASSERT_THROW( min_max_pixel_values(im8mb, mini, maxi), ArgumentErr );
  ASSERT_THROW( min_max_pixel_values(imfmb, minf, maxf), ArgumentErr );
}

TEST( Statistics, MinChannel ) {
  generate_data();
  EXPECT_EQ( min_channel_value(im8), 24 );
  ASSERT_TRUE( is_of_type<vw::uint8>( min_channel_value(im8) ) );
  EXPECT_EQ( min_channel_value(imf), -1.0 );
  ASSERT_TRUE( is_of_type<double>( min_channel_value(imf) ) );
  EXPECT_EQ( min_channel_value(im8ma), 0 );
  ASSERT_TRUE( is_of_type<vw::uint8>( min_channel_value(im8ma) ) );
  EXPECT_EQ( min_channel_value(imfma), -1.0 );
  ASSERT_TRUE( is_of_type<double>( min_channel_value(imfma) ) );
  ASSERT_THROW( min_channel_value(im8mb), ArgumentErr );
  ASSERT_THROW( min_channel_value(imfmb), ArgumentErr );
  EXPECT_EQ( min_channel_value(imrgb8), 24 );
  ASSERT_TRUE( is_of_type<vw::uint8>( min_channel_value(imrgb8) ) );
  EXPECT_EQ( min_channel_value(imrgbf), -1 );
  ASSERT_TRUE( is_of_type<double>( min_channel_value(imrgbf) ) );
  EXPECT_EQ( min_channel_value(imrgba8a), 15 );
  ASSERT_TRUE( is_of_type<vw::uint8>( min_channel_value(imrgba8a) ) );
  EXPECT_EQ( min_channel_value(imrgbafa), -0.5 );
  ASSERT_TRUE( is_of_type<double>( min_channel_value(imrgbafa) ) );
  EXPECT_EQ( min_channel_value(imrgba8b), 0 );
  EXPECT_EQ( min_channel_value(imrgbafb), -0.6 );
  EXPECT_EQ( min_channel_value(imrgba8c), 0 );
  EXPECT_EQ( min_channel_value(imrgbafc), -1 );
}

TEST( Statistics, MaxChannel ) {
  generate_data();
  EXPECT_EQ( max_channel_value(im8), 119 );
  ASSERT_TRUE( is_of_type<vw::uint8>( max_channel_value(im8) ) );
  EXPECT_EQ( max_channel_value(imf), 1.0 );
  ASSERT_TRUE( is_of_type<double>( max_channel_value(imf) ) );
  EXPECT_EQ( max_channel_value(im8ma), 119 );
  ASSERT_TRUE( is_of_type<vw::uint8>( max_channel_value(im8ma) ) );
  EXPECT_EQ( max_channel_value(imfma), 0 );
  ASSERT_TRUE( is_of_type<double>( max_channel_value(imfma) ) );
  ASSERT_THROW( max_channel_value(im8mb), ArgumentErr );
  ASSERT_THROW( max_channel_value(imfmb), ArgumentErr );
  EXPECT_EQ( max_channel_value(imrgb8), 228 );
  ASSERT_TRUE( is_of_type<vw::uint8>( max_channel_value(imrgb8) ) );
  EXPECT_EQ( max_channel_value(imrgbf), 1.6 );
  ASSERT_TRUE( is_of_type<double>( max_channel_value(imrgbf) ) );
  EXPECT_EQ( max_channel_value(imrgba8a), 255 );
  ASSERT_TRUE( is_of_type<vw::uint8>( max_channel_value(imrgba8a) ) );
  EXPECT_EQ( max_channel_value(imrgbafa), 1.6 );
  ASSERT_TRUE( is_of_type<double>( max_channel_value(imrgbafa) ) );
  EXPECT_EQ( max_channel_value(imrgba8b), 255 );
  EXPECT_EQ( max_channel_value(imrgbafb), 2.0 );
  EXPECT_EQ( max_channel_value(imrgba8c), 228 );
  EXPECT_EQ( max_channel_value(imrgbafc), 1.6 );
}

TEST( Statistics, MinMaxChannel ) {
  generate_data();
  vw::uint8 mini, maxi;
  double minf, maxf;
  min_max_channel_values(im8,mini,maxi);
  EXPECT_EQ( mini, 24 );
  EXPECT_EQ( maxi, 119 );
  min_max_channel_values(imf,minf,maxf);
  EXPECT_EQ( minf, -1.0 );
  EXPECT_EQ( maxf, 1.0 );
  min_max_channel_values(im8ma,mini,maxi);
  EXPECT_EQ( mini, 0 );
  EXPECT_EQ( maxi, 119 );
  min_max_channel_values(imfma,minf,maxf);
  EXPECT_EQ( minf, -1.0 );
  EXPECT_EQ( maxf, 0.0 );
  ASSERT_THROW( min_max_channel_values(im8mb,mini,maxi), ArgumentErr );
  ASSERT_THROW( min_max_channel_values(imfmb,minf,maxf), ArgumentErr );
  min_max_channel_values(imrgb8,mini,maxi);
  EXPECT_EQ( mini, 24 );
  EXPECT_EQ( maxi, 228 );
  min_max_channel_values(imrgbf,minf,maxf);
  EXPECT_EQ( minf, -1 );
  EXPECT_EQ( maxf, 1.6 );
  min_max_channel_values(imrgba8a,mini,maxi);
  EXPECT_EQ( mini, 15 );
  EXPECT_EQ( maxi, 255 );
  min_max_channel_values(imrgbafa,minf,maxf);
  EXPECT_EQ( minf, -0.5 );
  EXPECT_EQ( maxf, 1.6 );
  min_max_channel_values(imrgba8b,mini,maxi);
  EXPECT_EQ( mini, 0 );
  EXPECT_EQ( maxi, 255 );
  min_max_channel_values(imrgbafb,minf,maxf);
  EXPECT_EQ( minf, -0.6 );
  EXPECT_EQ( maxf, 2.0 );
  min_max_channel_values(imrgba8c,mini,maxi);
  EXPECT_EQ( mini, 0 );
  EXPECT_EQ( maxi, 228 );
  min_max_channel_values(imrgbafc,minf,maxf);
  EXPECT_EQ( minf, -1.0 );
  EXPECT_EQ( maxf, 1.6 );
}

TEST( Statistics, SumPixel ) {
  generate_data();
  EXPECT_EQ( sum_of_pixel_values(im8), 143 );
  ASSERT_TRUE( is_of_type<vw::int32>( sum_of_pixel_values(im8) ) );
  EXPECT_EQ( sum_of_pixel_values(imf), 0.0 );
  ASSERT_TRUE( is_of_type<double>( sum_of_pixel_values(imf) ) );
  EXPECT_EQ( sum_of_pixel_values(im8ma), 119 );
  ASSERT_TRUE( is_of_type<vw::int32>( sum_of_pixel_values(im8ma) ) );
  EXPECT_EQ( sum_of_pixel_values(imfma), -1.0 );
  ASSERT_TRUE( is_of_type<double>( sum_of_pixel_values(imfma) ) );
  EXPECT_EQ( sum_of_pixel_values(im8mb), 0 );
  EXPECT_EQ( sum_of_pixel_values(imfmb), 0.0 );
  EXPECT_PIXEL_EQ( sum_of_pixel_values(imrgb8), PixelRGB<int>(143,143,317) );
  ASSERT_TRUE( is_of_type<PixelRGB<vw::int32> >( sum_of_pixel_values(imrgb8) ) );
  EXPECT_PIXEL_EQ( sum_of_pixel_values(imrgbf), PixelRGB<double>(-0.5,0.9,2.6) );
  ASSERT_TRUE( is_of_type<PixelRGB<double> >( sum_of_pixel_values(imrgbf) ) );
  EXPECT_PIXEL_EQ( sum_of_pixel_values(imrgba8a), PixelRGBA<int>(48,82,135,306) );
  ASSERT_TRUE( is_of_type<PixelRGBA<vw::int32> >( sum_of_pixel_values(imrgba8a) ) );
  EXPECT_PIXEL_EQ( sum_of_pixel_values(imrgbafa), PixelRGBA<double>(0.0,0.9,2.1,1.5) );
  ASSERT_TRUE( is_of_type<PixelRGBA<double> >( sum_of_pixel_values(imrgbafa) ) );
  EXPECT_PIXEL_EQ( sum_of_pixel_values(imrgba8b), PixelRGBA<int>(44,82,189,255) );
  EXPECT_PIXEL_EQ( sum_of_pixel_values(imrgbafb), PixelRGBA<double>((0.5+-0.6),0.9,3.6,1.0) );
  EXPECT_PIXEL_EQ( sum_of_pixel_values(imrgba8c), PixelRGBA<int>(143,143,317,0) );
  EXPECT_PIXEL_EQ( sum_of_pixel_values(imrgbafc), PixelRGBA<double>(-0.5,0.9,2.6,0) );
}

TEST( Statistics, SumChannel ) {
  generate_data();
  EXPECT_EQ( sum_of_channel_values(im8), 143 );
  ASSERT_TRUE( is_of_type<vw::int32>( sum_of_channel_values(im8) ) );
  EXPECT_EQ( sum_of_channel_values(imf), 0.0 );
  ASSERT_TRUE( is_of_type<double>( sum_of_channel_values(imf) ) );
  EXPECT_EQ( sum_of_channel_values(im8ma), 119 );
  ASSERT_TRUE( is_of_type<vw::int32>( sum_of_channel_values(im8ma) ) );
  EXPECT_EQ( sum_of_channel_values(imfma), -1.0 );
  ASSERT_TRUE( is_of_type<double>( sum_of_channel_values(imfma) ) );
  EXPECT_EQ( sum_of_channel_values(im8mb), 0 );
  EXPECT_EQ( sum_of_channel_values(imfmb), 0.0 );
  EXPECT_EQ( sum_of_channel_values(imrgb8), 603 );
  ASSERT_TRUE( is_of_type<vw::int32>( sum_of_channel_values(imrgb8) ) );
  EXPECT_EQ( sum_of_channel_values(imrgbf), 3.0 );
  ASSERT_TRUE( is_of_type<double>( sum_of_channel_values(imrgbf) ) );
  EXPECT_EQ( sum_of_channel_values(imrgba8a), 571 );
  ASSERT_TRUE( is_of_type<vw::int32>( sum_of_channel_values(imrgba8a) ) );
  EXPECT_EQ( sum_of_channel_values(imrgbafa), 4.5 );
  ASSERT_TRUE( is_of_type<double>( sum_of_channel_values(imrgbafa) ) );
  EXPECT_EQ( sum_of_channel_values(imrgba8b), 570 );
  EXPECT_EQ( sum_of_channel_values(imrgbafb), 5.4 );
  EXPECT_EQ( sum_of_channel_values(imrgba8c), 603 );
  EXPECT_EQ( sum_of_channel_values(imrgbafc), 3.0 );
}

TEST( Statistics, MeanPixel ) {
  generate_data();
  double val;
  PixelRGB<double> rgb;
  PixelRGBA<double> rgba;
  val = mean_pixel_value(im8);
  EXPECT_NEAR( val, 71.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_pixel_value(im8) ) );
  val = mean_pixel_value(imf);
  EXPECT_NEAR( val, 0, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_pixel_value(imf) ) );
  val = mean_pixel_value(im8ma);
  EXPECT_NEAR( val, 59.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_pixel_value(im8ma) ) );
  val = mean_pixel_value(imfma);
  EXPECT_NEAR( val, -0.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_pixel_value(imfma) ) );
  ASSERT_THROW( mean_pixel_value(im8mb), ArgumentErr );
  ASSERT_THROW( mean_pixel_value(imfmb), ArgumentErr );
  EXPECT_PIXEL_NEAR( mean_pixel_value(imrgb8),
                     PixelRGB<double>(71.5,71.5,158.5), 1e-6 );
  ASSERT_TRUE( is_of_type<PixelRGB<double> >( mean_pixel_value(imrgb8) ) );
  EXPECT_PIXEL_NEAR( mean_pixel_value(imrgbf),
                     PixelRGB<double>(-0.25,0.45,1.3), 1e-6 );
  EXPECT_PIXEL_NEAR( mean_pixel_value(imrgbf ),
                     PixelRGB<double>(-0.25,0.45,1.3), 1e-6 );
  ASSERT_TRUE( is_of_type<PixelRGB<double> >( mean_pixel_value(imrgbf) ) );
  EXPECT_PIXEL_NEAR( mean_pixel_value(imrgba8a),
                     PixelRGBA<double>(24,41,67.5,153), 1e-6 );
  ASSERT_TRUE( is_of_type<PixelRGBA<double> >( mean_pixel_value(imrgba8a) ) );
  EXPECT_PIXEL_NEAR( mean_pixel_value(imrgbafa),
                     PixelRGBA<double>(0.0,0.45,1.05,0.75), 1e-6 );
  ASSERT_TRUE( is_of_type<PixelRGBA<double> >( mean_pixel_value(imrgbafa) ) );
  EXPECT_PIXEL_NEAR( mean_pixel_value(imrgba8b),
                     PixelRGBA<double>(22,41,94.5,127.5), 1e-6 );
  EXPECT_PIXEL_NEAR( mean_pixel_value(imrgbafb),
                     PixelRGBA<double>(-0.05,0.45,1.8,0.5), 1e-6 );
  EXPECT_PIXEL_NEAR( mean_pixel_value(imrgba8c),
                     PixelRGBA<double>(71.5,71.5,158.5,0), 1e-6 );
  EXPECT_PIXEL_NEAR( mean_pixel_value(imrgbafc),
                     PixelRGBA<double>(-0.25,0.45,1.3,0), 1e-6 );
}

TEST( Statistics, WeightedMeanPixel ) {
  generate_data();
  double val;
  PixelRGB<double> rgb;
  PixelRGBA<double> rgba;
  val = weighted_mean_pixel_value(im8);
  EXPECT_NEAR( val, 71.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_pixel_value(im8) ) );
  val = weighted_mean_pixel_value(imf);
  EXPECT_NEAR( val, 0, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_pixel_value(imf) ) );
  val = weighted_mean_pixel_value(im8ma);
  EXPECT_NEAR( val, 59.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_pixel_value(im8ma) ) );
  val = weighted_mean_pixel_value(imfma);
  EXPECT_NEAR( val, -0.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_pixel_value(imfma) ) );
  ASSERT_THROW( weighted_mean_pixel_value(im8mb), ArgumentErr );
  ASSERT_THROW( weighted_mean_pixel_value(imfmb), ArgumentErr );
  EXPECT_PIXEL_NEAR( weighted_mean_pixel_value(imrgb8),
                     PixelRGB<double>(71.5,71.5,158.5), 1e-6 );
  ASSERT_TRUE( is_of_type<PixelRGB<double> >( weighted_mean_pixel_value(imrgb8) ) );
  EXPECT_PIXEL_NEAR( weighted_mean_pixel_value(imrgbf),
                     PixelRGB<double>(-0.25,0.45,1.3), 1e-6 );
  ASSERT_TRUE( is_of_type<PixelRGB<double> >( weighted_mean_pixel_value(imrgbf) ) );
  EXPECT_PIXEL_NEAR( weighted_mean_pixel_value(imrgba8a),
                     PixelRGB<double>(40,68.333333,112.5), 1e-6 );
  ASSERT_TRUE( is_of_type<PixelRGB<double> >( weighted_mean_pixel_value(imrgba8a) ) );
  EXPECT_PIXEL_NEAR( weighted_mean_pixel_value(imrgbafa),
                     PixelRGB<double>(0,0.6,1.4), 1e-6 );
  ASSERT_TRUE( is_of_type<PixelRGB<double> >( weighted_mean_pixel_value(imrgbafa) ) );
  // These next two are confusing because the function assumes
  // pre-multiplied images so the test images are technically invalid
  EXPECT_PIXEL_NEAR( weighted_mean_pixel_value(imrgba8b),
                     PixelRGB<double>(44,82,189), 1e-6 );
  EXPECT_PIXEL_NEAR( weighted_mean_pixel_value(imrgbafb),
                     PixelRGB<double>(-0.1,0.9,3.6), 1e-6 );
  ASSERT_THROW( weighted_mean_pixel_value(imrgba8c), ArgumentErr );
  ASSERT_THROW( weighted_mean_pixel_value(imrgbafc), ArgumentErr );
}

TEST( Statistics, MeanChannel ) {
  generate_data();
  EXPECT_NEAR( mean_channel_value(im8), 71.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_channel_value(im8) ) );
  EXPECT_NEAR( mean_channel_value(imf), 0.0, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_channel_value(imf) ) );
  EXPECT_NEAR( mean_channel_value(im8ma), 59.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_channel_value(im8ma) ) );
  EXPECT_NEAR( mean_channel_value(imfma), -0.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_channel_value(imfma) ) );
  ASSERT_THROW( mean_channel_value(im8mb), ArgumentErr );
  ASSERT_THROW( mean_channel_value(imfmb), ArgumentErr );
  EXPECT_NEAR( mean_channel_value(imrgb8), 100.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_channel_value(imrgb8) ) );
  EXPECT_NEAR( mean_channel_value(imrgbf), 0.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_channel_value(imrgbf) ) );
  EXPECT_NEAR( mean_channel_value(imrgba8a), 71.375, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_channel_value(imrgba8a) ) );
  EXPECT_NEAR( mean_channel_value(imrgbafa), 0.5625, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( mean_channel_value(imrgbafa) ) );
  EXPECT_NEAR( mean_channel_value(imrgba8b), 71.25, 1e-6 );
  EXPECT_NEAR( mean_channel_value(imrgbafb), 0.675, 1e-6 );
  EXPECT_NEAR( mean_channel_value(imrgba8c), 75.375, 1e-6 );
  EXPECT_NEAR( mean_channel_value(imrgbafc), 0.375, 1e-6 );
}

TEST( Statistics, WeightedMeanChannel ) {
  generate_data();
  EXPECT_NEAR( weighted_mean_channel_value(im8), 71.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_channel_value(im8) ) );
  EXPECT_NEAR( weighted_mean_channel_value(imf), 0.0, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_channel_value(imf) ) );
  EXPECT_NEAR( weighted_mean_channel_value(im8ma), 59.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_channel_value(im8ma) ) );
  EXPECT_NEAR( weighted_mean_channel_value(imfma), -0.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_channel_value(imfma) ) );
  ASSERT_THROW( weighted_mean_channel_value(im8mb), ArgumentErr );
  ASSERT_THROW( weighted_mean_channel_value(imfmb), ArgumentErr );
  EXPECT_NEAR( weighted_mean_channel_value(imrgb8), 100.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_channel_value(imrgb8) ) );
  EXPECT_NEAR( weighted_mean_channel_value(imrgbf), 0.5, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_channel_value(imrgbf) ) );
  EXPECT_NEAR( weighted_mean_channel_value(imrgba8a), 73.6111111, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_channel_value(imrgba8a) ) );
  EXPECT_NEAR( weighted_mean_channel_value(imrgbafa), 0.6666667, 1e-6 );
  ASSERT_TRUE( is_of_type<double>( weighted_mean_channel_value(imrgbafa) ) );
  EXPECT_NEAR( weighted_mean_channel_value(imrgba8b), 105.0, 1e-6 );
  EXPECT_NEAR( weighted_mean_channel_value(imrgbafb), 1.4666667, 1e-6 );
  ASSERT_THROW( weighted_mean_channel_value(imrgba8c), ArgumentErr );
  ASSERT_THROW( weighted_mean_channel_value(imrgbafc), ArgumentErr );
}

TEST( Statistics, StdDevPixel ) {
  generate_data();
  EXPECT_EQ(stddev_pixel_value(im8), 47.5);
  ASSERT_TRUE( is_of_type<double>( stddev_pixel_value(im8) ) );
  EXPECT_EQ(stddev_pixel_value(imf), 1.0);
  ASSERT_TRUE( is_of_type<double>( stddev_pixel_value(imf) ) );
}

TEST( Statistics, StdDevChannel ) {
  generate_data();
  EXPECT_EQ(stddev_channel_value(im8), 47.5);
  ASSERT_TRUE( is_of_type<double>( stddev_channel_value(im8) ) );
  EXPECT_EQ(stddev_channel_value(imf), 1.0);
  ASSERT_TRUE( is_of_type<double>( stddev_channel_value(imf) ) );
  EXPECT_EQ(stddev_channel_value(im8ma), 59.5);
  ASSERT_TRUE( is_of_type<double>( stddev_channel_value(im8ma) ) );
  EXPECT_EQ(stddev_channel_value(imfma), 0.5);
  ASSERT_TRUE( is_of_type<double>( stddev_channel_value(imfma) ) );
  ASSERT_THROW(stddev_channel_value(im8mb), ArgumentErr);
  ASSERT_THROW(stddev_channel_value(im8mb), ArgumentErr);
  EXPECT_NEAR(stddev_channel_value(imrgb8), 63.6468, 1e-4);
  ASSERT_TRUE( is_of_type<double>( stddev_channel_value(imrgb8) ) );
  EXPECT_NEAR(stddev_channel_value(imrgbf), 0.828654, 1e-6);
  ASSERT_TRUE( is_of_type<double>( stddev_channel_value(imrgbf) ) );
  EXPECT_NEAR(stddev_channel_value(imrgba8a), 73.1213, 1e-4);
  ASSERT_TRUE( is_of_type<double>( stddev_channel_value(imrgba8a) ) );
  EXPECT_NEAR(stddev_channel_value(imrgbafa), 0.5956, 1e-4);
  ASSERT_TRUE( is_of_type<double>( stddev_channel_value(imrgbafa) ) );
  EXPECT_NEAR(stddev_channel_value(imrgba8b), 77.4786, 1e-4);
  EXPECT_NEAR(stddev_channel_value(imrgbafb), 0.8166, 1e-4);
  EXPECT_NEAR(stddev_channel_value(imrgba8c), 70.2280, 1e-4);
  EXPECT_NEAR(stddev_channel_value(imrgbafc), 0.7495, 1e-4);
}

TEST( Statistics, WeightedStdDevChannel ) {
  generate_data();
  /*
    EXPECT_EQ(weighted_stddev_channel_value(im8), 47.5);
    EXPECT_EQ(weighted_stddev_channel_value(imf), 1.0);
    EXPECT_NEAR(weighted_stddev_channel_value(imrgb8), 63.6468, 1e-4);
    EXPECT_NEAR(weighted_stddev_channel_value(imrgbf), 0.828654, 1e-6);
    EXPECT_NEAR(weighted_stddev_channel_value(imrgba8a), 47.3288, 1e-4);
    EXPECT_NEAR(weighted_stddev_channel_value(imrgbafa), 0.763035, 1e-6);
    EXPECT_NEAR(weighted_stddev_channel_value(imrgba8b), 26.9938, 1e-4);
    EXPECT_NEAR(weighted_stddev_channel_value(imrgbafb), 0.454606, 1e-6);
    ASSERT_THROW(weighted_stddev_channel_value(imrgba8c), ArgumentErr);
    ASSERT_THROW(weighted_stddev_channel_value(imrgbafc), ArgumentErr);
  */
}

TEST( Statistics, MedianChannel ) {
  generate_data();
  ImageView<PixelRGBA<vw::uint8> > image1(4,1);
  image1(0,0) = PixelRGBA<vw::uint8>(8,7,6,5);
  image1(1,0) = PixelRGBA<vw::uint8>(4,3,2,1);
  image1(2,0) = PixelRGBA<vw::uint8>(9,8,7,0);
  image1(3,0) = PixelRGBA<vw::uint8>(12,11,10,0);
  EXPECT_EQ( median_channel_value(image1), 6 ); // 6.5 actually if float
  ASSERT_TRUE( is_of_type<vw::uint8>( median_channel_value(image1) ) );
  ImageView<PixelMask<PixelRGB<vw::uint8> > > image2(4,1);
  image2(0,0) = PixelMask<PixelRGB<vw::uint8> >(8,7,6);
  image2(1,0) = PixelMask<PixelRGB<vw::uint8> >(4,3,2);
  image2(2,0) = PixelMask<PixelRGB<vw::uint8> >(9,8,7);
  image2(2,0).invalidate();
  image2(3,0) = PixelMask<PixelRGB<vw::uint8> >(12,11,10);
  image2(3,0).invalidate();
  EXPECT_EQ( median_channel_value(image2), 5 );
  ASSERT_TRUE( is_of_type<vw::uint8>( median_channel_value(image2) ) );
}

TEST( Statistics, Histogram ) {
  
  vw::math::Histogram hist;
  int num_bins = 256;
  histogram(im8b, num_bins, hist);

  EXPECT_EQ(hist.get_bin_value(  0), 2);
  EXPECT_EQ(hist.get_bin_value( 24), 1);
  EXPECT_EQ(hist.get_bin_value(119), 2);
  EXPECT_EQ(hist.get_bin_value(255), 1);
}

TEST( Statistics, OptimalThreshold ) {
  double t  = optimal_threshold(im8b);
  double t0 = 0.27843137254902; // Computed in Matlab, t0 = graythresh(uint8_im)
  EXPECT_NEAR(t, t0, 1e-15);
}

TEST(BlockOperations, DISABLED_CDF) {

  const Vector2i block_size(128, 128);
  int sumsample_amount = 1;

  // Generate a test image.
  size_t real_count = 0;
  const int size = 512;
  ImageView<uint8> image(size,size);
  for (int i=0; i<size; ++i) {
    for (int j=0; j<size; ++j) {
      uint8 value = i % 10;
      image(i,j) = value;
    }
  }

  // Compute the CDF in using a single thread.
  ChannelAccumulator<vw::math::CDFAccumulator<float> > normal_cdf;
  for_each_pixel( subsample( edge_extend(image, ConstantEdgeExtension()),
                             sumsample_amount ),
                  normal_cdf );

  // Compute the CDF in parallel.
  vw::math::CDFAccumulator<float> parallel_cdf;
  block_cdf_computation(image, parallel_cdf, sumsample_amount, block_size);

  // Check results.
  const float EPS = 0.11; // Larger sample sizes reduce this error
  EXPECT_NEAR(normal_cdf.quantile(0),          parallel_cdf.quantile(0), EPS);
  EXPECT_NEAR(normal_cdf.quantile(1),          parallel_cdf.quantile(1), EPS);
  EXPECT_NEAR(normal_cdf.approximate_mean  (), parallel_cdf.approximate_mean  (), EPS);
  EXPECT_NEAR(normal_cdf.approximate_stddev(), parallel_cdf.approximate_stddev(), EPS);
  EXPECT_NEAR(normal_cdf.quantile(0.02),       parallel_cdf.quantile(0.02), EPS);
  EXPECT_NEAR(normal_cdf.quantile(0.98),       parallel_cdf.quantile(0.98), EPS);
}
