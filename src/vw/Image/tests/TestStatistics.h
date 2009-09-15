// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// Image/tests/TestStatistics.h
#include <cxxtest/TestSuite.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>

using namespace vw;

class TestImageStatistics : public CxxTest::TestSuite
{
  ImageView<vw::uint8> im8;
  ImageView<double> imf;
  ImageView<PixelMask<vw::uint8> > im8ma;
  ImageView<PixelMask<double> > imfma;
  ImageView<PixelMask<vw::uint8> > im8mb;
  ImageView<PixelMask<double> > imfmb;
  ImageView<PixelRGB<vw::uint8> > imrgb8;
  ImageView<PixelRGB<double> > imrgbf;
  ImageView<PixelRGBA<vw::uint8> > imrgba8a;
  ImageView<PixelRGBA<double> > imrgbafa;
  ImageView<PixelRGBA<vw::uint8> > imrgba8b;
  ImageView<PixelRGBA<double> > imrgbafb;
  ImageView<PixelRGBA<vw::uint8> > imrgba8c;
  ImageView<PixelRGBA<double> > imrgbafc;

  template <class T1, class T2>
  static bool is_of_type( T2 ) {
    return boost::is_same<T1,T2>::value;
  }

public:

  TestImageStatistics()
    : im8(2,1), imf(2,1), im8ma(2,1), imfma(2,1), im8mb(2,1), imfmb(2,1),
      imrgb8(2,1), imrgbf(2,1), imrgba8a(2,1), imrgbafa(2,1),
      imrgba8b(2,1), imrgbafb(2,1), imrgba8c(2,1), imrgbafc(2,1)
  {
    im8(0,0) = 24;
    im8(1,0) = 119;

    imf(0,0) = 1.0;
    imf(1,0) = -1.0;

    im8ma(1,0) = 119;

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

  void test_min_pixel_value()
  {
    TS_ASSERT_EQUALS( min_pixel_value(channels_to_planes(imrgb8)), 24 );
    TS_ASSERT( is_of_type<vw::uint8>( min_pixel_value(channels_to_planes(imrgb8)) ) );
    TS_ASSERT_EQUALS( min_pixel_value(channels_to_planes(imrgbf)), -1.0 );
    TS_ASSERT( is_of_type<double>( min_pixel_value(channels_to_planes(imrgbf)) ) );
    TS_ASSERT_EQUALS( min_pixel_value(im8), 24 );
    TS_ASSERT( is_of_type<vw::uint8>( min_pixel_value(im8) ) );
    TS_ASSERT_EQUALS( min_pixel_value(imf), -1.0 );
    TS_ASSERT( is_of_type<double>( min_pixel_value(imf) ) );
    TS_ASSERT_EQUALS( min_pixel_value(im8ma), 119 );
    TS_ASSERT( is_of_type<vw::uint8>( min_pixel_value(im8ma) ) );
    TS_ASSERT_EQUALS( min_pixel_value(imfma), -1.0 );
    TS_ASSERT( is_of_type<double>( min_pixel_value(imfma) ) );
    TS_ASSERT_THROWS( min_pixel_value(im8mb), ArgumentErr );
    TS_ASSERT_THROWS( min_pixel_value(imfmb), ArgumentErr );
  }

  void test_max_pixel_value()
  {
    TS_ASSERT_EQUALS( max_pixel_value(channels_to_planes(imrgb8)), 228 );
    TS_ASSERT( is_of_type<vw::uint8>( max_pixel_value(channels_to_planes(imrgb8)) ) );
    TS_ASSERT_EQUALS( max_pixel_value(channels_to_planes(imrgbf)), 1.6 );
    TS_ASSERT( is_of_type<double>( max_pixel_value(channels_to_planes(imrgbf)) ) );
    TS_ASSERT_EQUALS( max_pixel_value(im8), 119 );
    TS_ASSERT( is_of_type<vw::uint8>( max_pixel_value(im8) ) );
    TS_ASSERT_EQUALS( max_pixel_value(imf), 1.0 );
    TS_ASSERT( is_of_type<double>( max_pixel_value(imf) ) );
    TS_ASSERT_EQUALS( max_pixel_value(im8ma), 119 );
    TS_ASSERT( is_of_type<vw::uint8>( max_pixel_value(im8ma) ) );
    TS_ASSERT_EQUALS( max_pixel_value(imfma), -1.0 );
    TS_ASSERT( is_of_type<double>( max_pixel_value(imfma) ) );
    TS_ASSERT_THROWS( max_pixel_value(im8mb), ArgumentErr );
    TS_ASSERT_THROWS( max_pixel_value(imfmb), ArgumentErr );
  }

  void test_min_max_pixel_values()
  {
    vw::uint8 mini, maxi;
    double minf, maxf;
    min_max_pixel_values( channels_to_planes(imrgb8), mini, maxi );
    TS_ASSERT_EQUALS( mini, 24 );
    TS_ASSERT_EQUALS( maxi, 228 );
    min_max_pixel_values( channels_to_planes(imrgbf), minf, maxf );
    TS_ASSERT_EQUALS( minf, -1.0 );
    TS_ASSERT_EQUALS( maxf, 1.6 );
    min_max_pixel_values( im8, mini, maxi );
    TS_ASSERT_EQUALS( mini, 24 );
    TS_ASSERT_EQUALS( maxi, 119 );
    min_max_pixel_values( imf, minf, maxf );
    TS_ASSERT_EQUALS( minf, -1.0 );
    TS_ASSERT_EQUALS( maxf, 1.0 );
    min_max_pixel_values( im8ma, mini, maxi );
    TS_ASSERT_EQUALS( mini, 119 );
    TS_ASSERT_EQUALS( maxi, 119 );
    min_max_pixel_values( imfma, minf, maxf );
    TS_ASSERT_EQUALS( minf, -1.0 );
    TS_ASSERT_EQUALS( maxf, -1.0 );
    TS_ASSERT_THROWS( min_max_pixel_values(im8mb, mini, maxi), ArgumentErr );
    TS_ASSERT_THROWS( min_max_pixel_values(imfmb, minf, maxf), ArgumentErr );
  }

  void test_min_channel_value()
  {
    TS_ASSERT_EQUALS( min_channel_value(im8), 24 );
    TS_ASSERT( is_of_type<vw::uint8>( min_channel_value(im8) ) );
    TS_ASSERT_EQUALS( min_channel_value(imf), -1.0 );
    TS_ASSERT( is_of_type<double>( min_channel_value(imf) ) );
    TS_ASSERT_EQUALS( min_channel_value(im8ma), 119 );
    TS_ASSERT( is_of_type<vw::uint8>( min_channel_value(im8ma) ) );
    TS_ASSERT_EQUALS( min_channel_value(imfma), -1.0 );
    TS_ASSERT( is_of_type<double>( min_channel_value(imfma) ) );
    TS_ASSERT_THROWS( min_channel_value(im8mb), ArgumentErr );
    TS_ASSERT_THROWS( min_channel_value(imfmb), ArgumentErr );
    TS_ASSERT_EQUALS( min_channel_value(imrgb8), 24 );
    TS_ASSERT( is_of_type<vw::uint8>( min_channel_value(imrgb8) ) );
    TS_ASSERT_EQUALS( min_channel_value(imrgbf), -1 );
    TS_ASSERT( is_of_type<double>( min_channel_value(imrgbf) ) );
    TS_ASSERT_EQUALS( min_channel_value(imrgba8a), 15 );
    TS_ASSERT( is_of_type<vw::uint8>( min_channel_value(imrgba8a) ) );
    TS_ASSERT_EQUALS( min_channel_value(imrgbafa), -0.5 );
    TS_ASSERT( is_of_type<double>( min_channel_value(imrgbafa) ) );
    TS_ASSERT_EQUALS( min_channel_value(imrgba8b), 0 );
    TS_ASSERT_EQUALS( min_channel_value(imrgbafb), -0.6 );
    TS_ASSERT_EQUALS( min_channel_value(imrgba8c), 0 );
    TS_ASSERT_EQUALS( min_channel_value(imrgbafc), -1 );
  }

  void test_max_channel_value()
  {
    TS_ASSERT_EQUALS( max_channel_value(im8), 119 );
    TS_ASSERT( is_of_type<vw::uint8>( max_channel_value(im8) ) );
    TS_ASSERT_EQUALS( max_channel_value(imf), 1.0 );
    TS_ASSERT( is_of_type<double>( max_channel_value(imf) ) );
    TS_ASSERT_EQUALS( max_channel_value(im8ma), 119 );
    TS_ASSERT( is_of_type<vw::uint8>( max_channel_value(im8ma) ) );
    TS_ASSERT_EQUALS( max_channel_value(imfma), -1.0 );
    TS_ASSERT( is_of_type<double>( max_channel_value(imfma) ) );
    TS_ASSERT_THROWS( max_channel_value(im8mb), ArgumentErr );
    TS_ASSERT_THROWS( max_channel_value(imfmb), ArgumentErr );
    TS_ASSERT_EQUALS( max_channel_value(imrgb8), 228 );
    TS_ASSERT( is_of_type<vw::uint8>( max_channel_value(imrgb8) ) );
    TS_ASSERT_EQUALS( max_channel_value(imrgbf), 1.6 );
    TS_ASSERT( is_of_type<double>( max_channel_value(imrgbf) ) );
    TS_ASSERT_EQUALS( max_channel_value(imrgba8a), 255 );
    TS_ASSERT( is_of_type<vw::uint8>( max_channel_value(imrgba8a) ) );
    TS_ASSERT_EQUALS( max_channel_value(imrgbafa), 1.6 );
    TS_ASSERT( is_of_type<double>( max_channel_value(imrgbafa) ) );
    TS_ASSERT_EQUALS( max_channel_value(imrgba8b), 255 );
    TS_ASSERT_EQUALS( max_channel_value(imrgbafb), 2.0 );
    TS_ASSERT_EQUALS( max_channel_value(imrgba8c), 228 );
    TS_ASSERT_EQUALS( max_channel_value(imrgbafc), 1.6 );
  }

  void test_min_max_channel_values()
  {
    vw::uint8 mini, maxi;
    double minf, maxf;
    min_max_channel_values(im8,mini,maxi);
    TS_ASSERT_EQUALS( mini, 24 );
    TS_ASSERT_EQUALS( maxi, 119 );
    min_max_channel_values(imf,minf,maxf);
    TS_ASSERT_EQUALS( minf, -1.0 );
    TS_ASSERT_EQUALS( maxf, 1.0 );
    min_max_channel_values(im8ma,mini,maxi);
    TS_ASSERT_EQUALS( mini, 119 );
    TS_ASSERT_EQUALS( maxi, 119 );
    min_max_channel_values(imfma,minf,maxf);
    TS_ASSERT_EQUALS( minf, -1.0 );
    TS_ASSERT_EQUALS( maxf, -1.0 );
    TS_ASSERT_THROWS( min_max_channel_values(im8mb,mini,maxi), ArgumentErr );
    TS_ASSERT_THROWS( min_max_channel_values(imfmb,minf,maxf), ArgumentErr );
    min_max_channel_values(imrgb8,mini,maxi);
    TS_ASSERT_EQUALS( mini, 24 );
    TS_ASSERT_EQUALS( maxi, 228 );
    min_max_channel_values(imrgbf,minf,maxf);
    TS_ASSERT_EQUALS( minf, -1 );
    TS_ASSERT_EQUALS( maxf, 1.6 );
    min_max_channel_values(imrgba8a,mini,maxi);
    TS_ASSERT_EQUALS( mini, 15 );
    TS_ASSERT_EQUALS( maxi, 255 );
    min_max_channel_values(imrgbafa,minf,maxf);
    TS_ASSERT_EQUALS( minf, -0.5 );
    TS_ASSERT_EQUALS( maxf, 1.6 );
    min_max_channel_values(imrgba8b,mini,maxi);
    TS_ASSERT_EQUALS( mini, 0 );
    TS_ASSERT_EQUALS( maxi, 255 );
    min_max_channel_values(imrgbafb,minf,maxf);
    TS_ASSERT_EQUALS( minf, -0.6 );
    TS_ASSERT_EQUALS( maxf, 2.0 );
    min_max_channel_values(imrgba8c,mini,maxi);
    TS_ASSERT_EQUALS( mini, 0 );
    TS_ASSERT_EQUALS( maxi, 228 );
    min_max_channel_values(imrgbafc,minf,maxf);
    TS_ASSERT_EQUALS( minf, -1.0 );
    TS_ASSERT_EQUALS( maxf, 1.6 );
  }

  void test_sum_of_pixel_values() {
    TS_ASSERT_EQUALS( sum_of_pixel_values(im8), 143 );
    TS_ASSERT( is_of_type<vw::int32>( sum_of_pixel_values(im8) ) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(imf), 0.0 );
    TS_ASSERT( is_of_type<double>( sum_of_pixel_values(imf) ) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(im8ma), 119 );
    TS_ASSERT( is_of_type<vw::int32>( sum_of_pixel_values(im8ma) ) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(imfma), -1.0 );
    TS_ASSERT( is_of_type<double>( sum_of_pixel_values(imfma) ) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(im8mb), 0 );
    TS_ASSERT_EQUALS( sum_of_pixel_values(imfmb), 0.0 );
    TS_ASSERT_EQUALS( sum_of_pixel_values(imrgb8), PixelRGB<double>(143,143,317) );
    TS_ASSERT( is_of_type<PixelRGB<vw::int32> >( sum_of_pixel_values(imrgb8) ) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(imrgbf), PixelRGB<double>(-0.5,0.9,2.6) );
    TS_ASSERT( is_of_type<PixelRGB<double> >( sum_of_pixel_values(imrgbf) ) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(imrgba8a), PixelRGBA<double>(48,82,135,306) );
    TS_ASSERT( is_of_type<PixelRGBA<vw::int32> >( sum_of_pixel_values(imrgba8a) ) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(imrgbafa), PixelRGBA<double>(0.0,0.9,2.1,1.5) );
    TS_ASSERT( is_of_type<PixelRGBA<double> >( sum_of_pixel_values(imrgbafa) ) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(imrgba8b), PixelRGBA<double>(44,82,189,255) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(imrgbafb), PixelRGBA<double>((0.5+-0.6),0.9,3.6,1.0) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(imrgba8c), PixelRGBA<double>(143,143,317,0) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(imrgbafc), PixelRGBA<double>(-0.5,0.9,2.6,0) );
  }

  void test_sum_of_channel_values() {
    TS_ASSERT_EQUALS( sum_of_channel_values(im8), 143 );
    TS_ASSERT( is_of_type<vw::int32>( sum_of_channel_values(im8) ) );
    TS_ASSERT_EQUALS( sum_of_channel_values(imf), 0.0 );
    TS_ASSERT( is_of_type<double>( sum_of_channel_values(imf) ) );
    TS_ASSERT_EQUALS( sum_of_channel_values(im8ma), 119 );
    TS_ASSERT( is_of_type<vw::int32>( sum_of_channel_values(im8ma) ) );
    TS_ASSERT_EQUALS( sum_of_channel_values(imfma), -1.0 );
    TS_ASSERT( is_of_type<double>( sum_of_channel_values(imfma) ) );
    TS_ASSERT_EQUALS( sum_of_channel_values(im8mb), 0 );
    TS_ASSERT_EQUALS( sum_of_channel_values(imfmb), 0.0 );
    TS_ASSERT_EQUALS( sum_of_channel_values(imrgb8), 603 );
    TS_ASSERT( is_of_type<vw::int32>( sum_of_channel_values(imrgb8) ) );
    TS_ASSERT_EQUALS( sum_of_channel_values(imrgbf), 3.0 );
    TS_ASSERT( is_of_type<double>( sum_of_channel_values(imrgbf) ) );
    TS_ASSERT_EQUALS( sum_of_channel_values(imrgba8a), 571 );
    TS_ASSERT( is_of_type<vw::int32>( sum_of_channel_values(imrgba8a) ) );
    TS_ASSERT_EQUALS( sum_of_channel_values(imrgbafa), 4.5 );
    TS_ASSERT( is_of_type<double>( sum_of_channel_values(imrgbafa) ) );
    TS_ASSERT_EQUALS( sum_of_channel_values(imrgba8b), 570 );
    TS_ASSERT_EQUALS( sum_of_channel_values(imrgbafb), 5.4 );
    TS_ASSERT_EQUALS( sum_of_channel_values(imrgba8c), 603 );
    TS_ASSERT_EQUALS( sum_of_channel_values(imrgbafc), 3.0 );
  }

  void test_mean_pixel_value()
  {
    double val;
    PixelRGB<double> rgb;
    PixelRGBA<double> rgba;
    val = mean_pixel_value(im8);
    TS_ASSERT_DELTA( val, 71.5, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_pixel_value(im8) ) );
    val = mean_pixel_value(imf);
    TS_ASSERT_DELTA( val, 0, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_pixel_value(imf) ) );
    val = mean_pixel_value(im8ma);
    TS_ASSERT_DELTA( val, 119, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_pixel_value(im8ma) ) );
    val = mean_pixel_value(imfma);
    TS_ASSERT_DELTA( val, -1, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_pixel_value(imfma) ) );
    TS_ASSERT_THROWS( mean_pixel_value(im8mb), ArgumentErr );
    TS_ASSERT_THROWS( mean_pixel_value(imfmb), ArgumentErr );
    rgb = mean_pixel_value(imrgb8);
    TS_ASSERT_DELTA( rgb.r(), 71.5, 1e-6 );
    TS_ASSERT_DELTA( rgb.g(), 71.5, 1e-6 );
    TS_ASSERT_DELTA( rgb.b(), 158.5, 1e-6 );
    TS_ASSERT( is_of_type<PixelRGB<double> >( mean_pixel_value(imrgb8) ) );
    rgb = mean_pixel_value(imrgbf);
    TS_ASSERT_DELTA( rgb.r(), -0.25, 1e-6 );
    TS_ASSERT_DELTA( rgb.g(), 0.45, 1e-6 );
    TS_ASSERT_DELTA( rgb.b(), 1.3, 1e-6 );
    TS_ASSERT( is_of_type<PixelRGB<double> >( mean_pixel_value(imrgbf) ) );
    rgba = mean_pixel_value(imrgba8a);
    TS_ASSERT_DELTA( rgba.r(), 24, 1e-6 );
    TS_ASSERT_DELTA( rgba.g(), 41, 1e-6 );
    TS_ASSERT_DELTA( rgba.b(), 67.5, 1e-6 );
    TS_ASSERT_DELTA( rgba.a(), 153, 1e-6 );
    TS_ASSERT( is_of_type<PixelRGBA<double> >( mean_pixel_value(imrgba8a) ) );
    rgba = mean_pixel_value(imrgbafa);
    TS_ASSERT_DELTA( rgba.r(), 0.0, 1e-6 );
    TS_ASSERT_DELTA( rgba.g(), 0.45, 1e-6 );
    TS_ASSERT_DELTA( rgba.b(), 1.05, 1e-6 );
    TS_ASSERT_DELTA( rgba.a(), 0.75, 1e-6 );
    TS_ASSERT( is_of_type<PixelRGBA<double> >( mean_pixel_value(imrgbafa) ) );
    rgba = mean_pixel_value(imrgba8b);
    TS_ASSERT_DELTA( rgba.r(), 22, 1e-6 );
    TS_ASSERT_DELTA( rgba.g(), 41, 1e-6 );
    TS_ASSERT_DELTA( rgba.b(), 94.5, 1e-6 );
    TS_ASSERT_DELTA( rgba.a(), 127.5, 1e-6 );
    rgba = mean_pixel_value(imrgbafb);
    TS_ASSERT_DELTA( rgba.r(), -0.05, 1e-6 );
    TS_ASSERT_DELTA( rgba.g(), 0.45, 1e-6 );
    TS_ASSERT_DELTA( rgba.b(), 1.8, 1e-6 );
    TS_ASSERT_DELTA( rgba.a(), 0.5, 1e-6 );
    rgba = mean_pixel_value(imrgba8c);
    TS_ASSERT_DELTA( rgba.r(), 71.5, 1e-6 );
    TS_ASSERT_DELTA( rgba.g(), 71.5, 1e-6 );
    TS_ASSERT_DELTA( rgba.b(), 158.5, 1e-6 );
    TS_ASSERT_DELTA( rgba.a(), 0, 1e-6 );
    rgba = mean_pixel_value(imrgbafc);
    TS_ASSERT_DELTA( rgba.r(), -0.25, 1e-6 );
    TS_ASSERT_DELTA( rgba.g(), 0.45, 1e-6 );
    TS_ASSERT_DELTA( rgba.b(), 1.3, 1e-6 );
    TS_ASSERT_DELTA( rgba.a(), 0, 1e-6 );
  }

  void test_weighted_mean_pixel_value()
  {
    double val;
    PixelRGB<double> rgb;
    PixelRGBA<double> rgba;
    val = weighted_mean_pixel_value(im8);
    TS_ASSERT_DELTA( val, 71.5, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_pixel_value(im8) ) );
    val = weighted_mean_pixel_value(imf);
    TS_ASSERT_DELTA( val, 0, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_pixel_value(imf) ) );
    val = weighted_mean_pixel_value(im8ma);
    TS_ASSERT_DELTA( val, 119, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_pixel_value(im8ma) ) );
    val = weighted_mean_pixel_value(imfma);
    TS_ASSERT_DELTA( val, -1, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_pixel_value(imfma) ) );
    TS_ASSERT_THROWS( weighted_mean_pixel_value(im8mb), ArgumentErr );
    TS_ASSERT_THROWS( weighted_mean_pixel_value(imfmb), ArgumentErr );
    rgb = weighted_mean_pixel_value(imrgb8);
    TS_ASSERT_DELTA( rgb.r(), 71.5, 1e-6 );
    TS_ASSERT_DELTA( rgb.g(), 71.5, 1e-6 );
    TS_ASSERT_DELTA( rgb.b(), 158.5, 1e-6 );
    TS_ASSERT( is_of_type<PixelRGB<double> >( weighted_mean_pixel_value(imrgb8) ) );
    rgb = weighted_mean_pixel_value(imrgbf);
    TS_ASSERT_DELTA( rgb.r(), -0.25, 1e-6 );
    TS_ASSERT_DELTA( rgb.g(), 0.45, 1e-6 );
    TS_ASSERT_DELTA( rgb.b(), 1.3, 1e-6 );
    TS_ASSERT( is_of_type<PixelRGB<double> >( weighted_mean_pixel_value(imrgbf) ) );
    rgb = weighted_mean_pixel_value(imrgba8a);
    TS_ASSERT_DELTA( rgb.r(), 40, 1e-6 );
    TS_ASSERT_DELTA( rgb.g(), 68.3333333, 1e-6 );
    TS_ASSERT_DELTA( rgb.b(), 112.5, 1e-6 );
    TS_ASSERT( is_of_type<PixelRGB<double> >( weighted_mean_pixel_value(imrgba8a) ) );
    rgb = weighted_mean_pixel_value(imrgbafa);
    TS_ASSERT_DELTA( rgb.r(), 0, 1e-6 );
    TS_ASSERT_DELTA( rgb.g(), 0.6, 1e-6 );
    TS_ASSERT_DELTA( rgb.b(), 1.4, 1e-6 );
    TS_ASSERT( is_of_type<PixelRGB<double> >( weighted_mean_pixel_value(imrgbafa) ) );
    // These next two are confusing because the function assumes
    // pre-multiplied images so the test images are technically invalid
    rgb = weighted_mean_pixel_value(imrgba8b);
    TS_ASSERT_DELTA( rgb.r(), 44, 1e-6 );
    TS_ASSERT_DELTA( rgb.g(), 82, 1e-6 );
    TS_ASSERT_DELTA( rgb.b(), 189, 1e-6 );
    rgb = weighted_mean_pixel_value(imrgbafb);
    TS_ASSERT_DELTA( rgb.r(), -0.1, 1e-6 );
    TS_ASSERT_DELTA( rgb.g(), 0.9, 1e-6 );
    TS_ASSERT_DELTA( rgb.b(), 3.6, 1e-6 );
    TS_ASSERT_THROWS( weighted_mean_pixel_value(imrgba8c), ArgumentErr );
    TS_ASSERT_THROWS( weighted_mean_pixel_value(imrgbafc), ArgumentErr );
  }

  void test_mean_channel_value()
  {
    TS_ASSERT_DELTA( mean_channel_value(im8), 71.5, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_channel_value(im8) ) );
    TS_ASSERT_DELTA( mean_channel_value(imf), 0.0, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_channel_value(imf) ) );
    TS_ASSERT_DELTA( mean_channel_value(im8ma), 119, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_channel_value(im8ma) ) );
    TS_ASSERT_DELTA( mean_channel_value(imfma), -1.0, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_channel_value(imfma) ) );
    TS_ASSERT_THROWS( mean_channel_value(im8mb), ArgumentErr );
    TS_ASSERT_THROWS( mean_channel_value(imfmb), ArgumentErr );
    TS_ASSERT_DELTA( mean_channel_value(imrgb8), 100.5, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_channel_value(imrgb8) ) );
    TS_ASSERT_DELTA( mean_channel_value(imrgbf), 0.5, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_channel_value(imrgbf) ) );
    TS_ASSERT_DELTA( mean_channel_value(imrgba8a), 71.375, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_channel_value(imrgba8a) ) );
    TS_ASSERT_DELTA( mean_channel_value(imrgbafa), 0.5625, 1e-6 );
    TS_ASSERT( is_of_type<double>( mean_channel_value(imrgbafa) ) );
    TS_ASSERT_DELTA( mean_channel_value(imrgba8b), 71.25, 1e-6 );
    TS_ASSERT_DELTA( mean_channel_value(imrgbafb), 0.675, 1e-6 );
    TS_ASSERT_DELTA( mean_channel_value(imrgba8c), 75.375, 1e-6 );
    TS_ASSERT_DELTA( mean_channel_value(imrgbafc), 0.375, 1e-6 );
  }

  void test_weighted_mean_channel_value()
  {
    TS_ASSERT_DELTA( weighted_mean_channel_value(im8), 71.5, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_channel_value(im8) ) );
    TS_ASSERT_DELTA( weighted_mean_channel_value(imf), 0.0, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_channel_value(imf) ) );
    TS_ASSERT_DELTA( weighted_mean_channel_value(im8ma), 119, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_channel_value(im8ma) ) );
    TS_ASSERT_DELTA( weighted_mean_channel_value(imfma), -1.0, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_channel_value(imfma) ) );
    TS_ASSERT_THROWS( weighted_mean_channel_value(im8mb), ArgumentErr );
    TS_ASSERT_THROWS( weighted_mean_channel_value(imfmb), ArgumentErr );
    TS_ASSERT_DELTA( weighted_mean_channel_value(imrgb8), 100.5, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_channel_value(imrgb8) ) );
    TS_ASSERT_DELTA( weighted_mean_channel_value(imrgbf), 0.5, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_channel_value(imrgbf) ) );
    TS_ASSERT_DELTA( weighted_mean_channel_value(imrgba8a), 73.6111111, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_channel_value(imrgba8a) ) );
    TS_ASSERT_DELTA( weighted_mean_channel_value(imrgbafa), 0.6666667, 1e-6 );
    TS_ASSERT( is_of_type<double>( weighted_mean_channel_value(imrgbafa) ) );
    TS_ASSERT_DELTA( weighted_mean_channel_value(imrgba8b), 105.0, 1e-6 );
    TS_ASSERT_DELTA( weighted_mean_channel_value(imrgbafb), 1.4666667, 1e-6 );
    TS_ASSERT_THROWS( weighted_mean_channel_value(imrgba8c), ArgumentErr );
    TS_ASSERT_THROWS( weighted_mean_channel_value(imrgbafc), ArgumentErr );
  }

  void test_stddev_pixel_value()
  {
    TS_ASSERT_EQUALS(stddev_pixel_value(im8), 47.5);
    TS_ASSERT( is_of_type<double>( stddev_pixel_value(im8) ) );
    TS_ASSERT_EQUALS(stddev_pixel_value(imf), 1.0);
    TS_ASSERT( is_of_type<double>( stddev_pixel_value(imf) ) );
  }

  void test_stddev_channel_value()
  {
    TS_ASSERT_EQUALS(stddev_channel_value(im8), 47.5);
    TS_ASSERT( is_of_type<double>( stddev_channel_value(im8) ) );
    TS_ASSERT_EQUALS(stddev_channel_value(imf), 1.0);
    TS_ASSERT( is_of_type<double>( stddev_channel_value(imf) ) );
    TS_ASSERT_EQUALS(stddev_channel_value(im8ma), 0);
    TS_ASSERT( is_of_type<double>( stddev_channel_value(im8ma) ) );
    TS_ASSERT_EQUALS(stddev_channel_value(imfma), 0);
    TS_ASSERT( is_of_type<double>( stddev_channel_value(imfma) ) );
    TS_ASSERT_THROWS(stddev_channel_value(im8mb), ArgumentErr);
    TS_ASSERT_THROWS(stddev_channel_value(im8mb), ArgumentErr);
    TS_ASSERT_DELTA(stddev_channel_value(imrgb8), 63.6468, 1e-4);
    TS_ASSERT( is_of_type<double>( stddev_channel_value(imrgb8) ) );
    TS_ASSERT_DELTA(stddev_channel_value(imrgbf), 0.828654, 1e-6);
    TS_ASSERT( is_of_type<double>( stddev_channel_value(imrgbf) ) );
    TS_ASSERT_DELTA(stddev_channel_value(imrgba8a), 73.1213, 1e-4);
    TS_ASSERT( is_of_type<double>( stddev_channel_value(imrgba8a) ) );
    TS_ASSERT_DELTA(stddev_channel_value(imrgbafa), 0.5956, 1e-4);
    TS_ASSERT( is_of_type<double>( stddev_channel_value(imrgbafa) ) );
    TS_ASSERT_DELTA(stddev_channel_value(imrgba8b), 77.4786, 1e-4);
    TS_ASSERT_DELTA(stddev_channel_value(imrgbafb), 0.8166, 1e-4);
    TS_ASSERT_DELTA(stddev_channel_value(imrgba8c), 70.2280, 1e-4);
    TS_ASSERT_DELTA(stddev_channel_value(imrgbafc), 0.7495, 1e-4);
  }

  void test_weighted_stddev_channel_value()
  {
    /*
    TS_ASSERT_EQUALS(weighted_stddev_channel_value(im8), 47.5);
    TS_ASSERT_EQUALS(weighted_stddev_channel_value(imf), 1.0);
    TS_ASSERT_DELTA(weighted_stddev_channel_value(imrgb8), 63.6468, 1e-4);
    TS_ASSERT_DELTA(weighted_stddev_channel_value(imrgbf), 0.828654, 1e-6);
    TS_ASSERT_DELTA(weighted_stddev_channel_value(imrgba8a), 47.3288, 1e-4);
    TS_ASSERT_DELTA(weighted_stddev_channel_value(imrgbafa), 0.763035, 1e-6);
    TS_ASSERT_DELTA(weighted_stddev_channel_value(imrgba8b), 26.9938, 1e-4);
    TS_ASSERT_DELTA(weighted_stddev_channel_value(imrgbafb), 0.454606, 1e-6);
    TS_ASSERT_THROWS(weighted_stddev_channel_value(imrgba8c), ArgumentErr);
    TS_ASSERT_THROWS(weighted_stddev_channel_value(imrgbafc), ArgumentErr);
    */
  }

  void test_median_channel_value()
  {
    ImageView<PixelRGBA<vw::uint8> > image1(4,1);
    image1(0,0) = PixelRGBA<vw::uint8>(8,7,6,5);
    image1(1,0) = PixelRGBA<vw::uint8>(4,3,2,1);
    image1(2,0) = PixelRGBA<vw::uint8>(9,8,7,0);
    image1(3,0) = PixelRGBA<vw::uint8>(12,11,10,0);
    TS_ASSERT_EQUALS( median_channel_value(image1), 7 );
    TS_ASSERT( is_of_type<vw::uint8>( median_channel_value(image1) ) );
    ImageView<PixelMask<PixelRGB<vw::uint8> > > image2(4,1);
    image2(0,0) = PixelMask<PixelRGB<vw::uint8> >(8,7,6);
    image2(1,0) = PixelMask<PixelRGB<vw::uint8> >(4,3,2);
    image2(2,0) = PixelMask<PixelRGB<vw::uint8> >(9,8,7);
    image2(2,0).invalidate();
    image2(3,0) = PixelMask<PixelRGB<vw::uint8> >(12,11,10);
    image2(3,0).invalidate();
    TS_ASSERT_EQUALS( median_channel_value(image2), 6 );
    TS_ASSERT( is_of_type<vw::uint8>( median_channel_value(image2) ) );
  }

};


