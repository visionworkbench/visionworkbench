// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

// Image/tests/TestStatistics.h
#include <cxxtest/TestSuite.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>

using namespace std;
using namespace vw;

class TestImageStatistics : public CxxTest::TestSuite
{
  ImageView<PixelRGB<double> > im1;
  ImageView<PixelRGB<uint8> > im2;
  ImageView<PixelRGBA<double> > im3;
  ImageView<PixelRGBA<uint8> > im4;
  ImageView<PixelRGBA<double> > im5;
  ImageView<PixelRGBA<uint8> > im6;
  ImageView<PixelRGBA<double> > im7;
  ImageView<PixelRGBA<uint8> > im8;
  ImageView<PixelRGBA<double> > im9;
  ImageView<PixelRGBA<uint8> > im10;

public:

  TestImageStatistics()
    : im1(2,1), im2(2,1), im3(2,1), im4(2,1), im5(2,1),
      im6(2,1), im7(2,1), im8(2,1), im9(2,1), im10(2,1)
  {
    im1(0,0) = PixelRGB<double>(0.5,0.9,1.6);
    im1(1,0) = PixelRGB<double>(-1,0,1);

    im2(0,0) = PixelRGB<uint8>(24,67,89);
    im2(1,0) = PixelRGB<uint8>(119,76,228);

    im3(0,0) = PixelRGBA<double>(0.5,0.9,1.6,1);
    im3(1,0) = PixelRGBA<double>(-1,0,1,1);

    im4(0,0) = PixelRGBA<uint8>(24,67,89,255);
    im4(1,0) = PixelRGBA<uint8>(119,76,228,255);

    im5(0,0) = PixelRGBA<double>(0.5,0.9,1.6,1.0);
    im5(1,0) = PixelRGBA<double>(-0.5,0,0.5,0.5);

    im6(0,0) = PixelRGBA<uint8>(24,67,89,255);
    im6(1,0) = PixelRGBA<uint8>(24,15,46,51);

    im7(0,0) = PixelRGBA<double>(0.5,0.9,1.6,1.0);
    im7(1,0) = PixelRGBA<double>(-0.6,0,2.0,0);

    im8(0,0) = PixelRGBA<uint8>(24,67,89,255);
    im8(1,0) = PixelRGBA<uint8>(20,15,100,0);

    im9(0,0) = PixelRGBA<double>(0.5,0.9,1.6,0);
    im9(1,0) = PixelRGBA<double>(-1,0,1,0);

    im10(0,0) = PixelRGBA<uint8>(24,67,89,0);
    im10(1,0) = PixelRGBA<uint8>(119,76,228,0);
  }

  void test_min_pixel_value()
  {
    TS_ASSERT_EQUALS( min_pixel_value(channels_to_planes(im1)), -1.0 );
    TS_ASSERT_EQUALS( min_pixel_value(channels_to_planes(im2)), 24 );
  }

  void test_max_pixel_value()
  {
    TS_ASSERT_EQUALS( max_pixel_value(channels_to_planes(im1)), 1.6 );
    TS_ASSERT_EQUALS( max_pixel_value(channels_to_planes(im2)), 228 );
  }

  void test_min_max_pixel_values()
  {
    double minf, maxf;
    uint8 mini, maxi;
    min_max_pixel_values( channels_to_planes(im1), minf, maxf );
    TS_ASSERT_EQUALS( minf, -1.0 );
    TS_ASSERT_EQUALS( maxf, 1.6 );
    min_max_pixel_values( channels_to_planes(im2), mini, maxi );
    TS_ASSERT_EQUALS( mini, 24 );
    TS_ASSERT_EQUALS( maxi, 228 );
  }

  void test_min_channel_value()
  {
    TS_ASSERT_EQUALS( min_channel_value(im1), -1 );
    TS_ASSERT_EQUALS( min_channel_value(im2), 24 );
    TS_ASSERT_EQUALS( min_channel_value(im3), -1 );
    TS_ASSERT_EQUALS( min_channel_value(im4), 24 );
    TS_ASSERT_EQUALS( min_channel_value(im5), -0.5 );
    TS_ASSERT_EQUALS( min_channel_value(im6), 15 );
    TS_ASSERT_EQUALS( min_channel_value(im7), 0.5 );
    TS_ASSERT_EQUALS( min_channel_value(im8), 24 );
    TS_ASSERT_THROWS( min_channel_value(im9), ArgumentErr );
    TS_ASSERT_THROWS( min_channel_value(im10), ArgumentErr );
  }

  void test_max_channel_value()
  {
    TS_ASSERT_EQUALS( max_channel_value(im1), 1.6 );
    TS_ASSERT_EQUALS( max_channel_value(im2), 228 );
    TS_ASSERT_EQUALS( max_channel_value(im3), 1.6 );
    TS_ASSERT_EQUALS( max_channel_value(im4), 228 );
    TS_ASSERT_EQUALS( max_channel_value(im5), 1.6 );
    TS_ASSERT_EQUALS( max_channel_value(im6), 89 );
    TS_ASSERT_EQUALS( max_channel_value(im7), 1.6 );
    TS_ASSERT_EQUALS( max_channel_value(im8), 89 );
    TS_ASSERT_THROWS( max_channel_value(im9), ArgumentErr );
    TS_ASSERT_THROWS( max_channel_value(im10), ArgumentErr );
  }

  void test_min_max_channel_values()
  {
    double minf, maxf;
    uint8 mini, maxi;
    min_max_channel_values(im1,minf,maxf);
    TS_ASSERT_EQUALS( minf, -1 );
    TS_ASSERT_EQUALS( maxf, 1.6 );
    min_max_channel_values(im2,mini,maxi);
    TS_ASSERT_EQUALS( mini, 24 );
    TS_ASSERT_EQUALS( maxi, 228 );
    min_max_channel_values(im3,minf,maxf);
    TS_ASSERT_EQUALS( minf, -1 );
    TS_ASSERT_EQUALS( maxf, 1.6 );
    min_max_channel_values(im4,mini,maxi);
    TS_ASSERT_EQUALS( mini, 24 );
    TS_ASSERT_EQUALS( maxi, 228 );
    min_max_channel_values(im5,minf,maxf);
    TS_ASSERT_EQUALS( minf, -0.5 );
    TS_ASSERT_EQUALS( maxf, 1.6 );
    min_max_channel_values(im6,mini,maxi);
    TS_ASSERT_EQUALS( mini, 15 );
    TS_ASSERT_EQUALS( maxi, 89 );
    min_max_channel_values(im7,minf,maxf);
    TS_ASSERT_EQUALS( minf, 0.5 );
    TS_ASSERT_EQUALS( maxf, 1.6 );
    min_max_channel_values(im8,mini,maxi);
    TS_ASSERT_EQUALS( mini, 24 );
    TS_ASSERT_EQUALS( maxi, 89 );
    TS_ASSERT_THROWS( min_max_channel_values(im9,minf,maxf), ArgumentErr );
    TS_ASSERT_THROWS( min_max_channel_values(im10,mini,maxi), ArgumentErr );
  }

  void test_sum_of_pixel_values() {
    TS_ASSERT_EQUALS( sum_of_pixel_values(im1), PixelRGB<double>(-0.5,0.9,2.6) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(im2), PixelRGB<double>(143,143,317) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(im3), PixelRGBA<double>(-0.5,0.9,2.6,2.0) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(im4), PixelRGBA<double>(143,143,317,510) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(im5), PixelRGBA<double>(0.0,0.9,2.1,1.5) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(im6), PixelRGBA<double>(48,82,135,306) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(im7), PixelRGBA<double>(0.5,0.9,1.6,1.0) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(im8), PixelRGBA<double>(24,67,89,255) );
    TS_ASSERT_EQUALS( sum_of_pixel_values(im9), PixelRGBA<double>() );
    TS_ASSERT_EQUALS( sum_of_pixel_values(im10), PixelRGBA<double>() );
  }

  void test_sum_of_channel_values() {
    TS_ASSERT_EQUALS( sum_of_channel_values(im1), 3.0 );
    TS_ASSERT_EQUALS( sum_of_channel_values(im2), 603 );
    TS_ASSERT_EQUALS( sum_of_channel_values(im3), 3.0 );
    TS_ASSERT_EQUALS( sum_of_channel_values(im4), 603 );
    TS_ASSERT_EQUALS( sum_of_channel_values(im5), 3.0 );
    TS_ASSERT_EQUALS( sum_of_channel_values(im6), 265 );
    TS_ASSERT_EQUALS( sum_of_channel_values(im7), 3.0 );
    TS_ASSERT_EQUALS( sum_of_channel_values(im8), 180 );
    TS_ASSERT_EQUALS( sum_of_channel_values(im9), 0 );
    TS_ASSERT_EQUALS( sum_of_channel_values(im10), 0 );
  }

  void test_mean_pixel_value()
  {
    PixelRGB<double> val;
    PixelRGBA<double> vala;
    val = mean_pixel_value(im1);
    TS_ASSERT_DELTA( val.r(), -0.25, 1e-6 );
    TS_ASSERT_DELTA( val.g(), 0.45, 1e-6 );
    TS_ASSERT_DELTA( val.b(), 1.3, 1e-6 );
    val = mean_pixel_value(im2);
    TS_ASSERT_DELTA( val.r(), 71.5, 1e-4 );
    TS_ASSERT_DELTA( val.g(), 71.5, 1e-4 );
    TS_ASSERT_DELTA( val.b(), 158.5, 1e-4 );
    vala = mean_pixel_value(im3);
    TS_ASSERT_DELTA( vala.r(), -0.25, 1e-6 );
    TS_ASSERT_DELTA( vala.g(), 0.45, 1e-6 );
    TS_ASSERT_DELTA( vala.b(), 1.3, 1e-6 );
    TS_ASSERT_DELTA( vala.a(), 1.0, 1e-6 );
    vala = mean_pixel_value(im4);
    TS_ASSERT_DELTA( vala.r(), 71.5, 1e-4 );
    TS_ASSERT_DELTA( vala.g(), 71.5, 1e-4 );
    TS_ASSERT_DELTA( vala.b(), 158.5, 1e-4 );
    TS_ASSERT_DELTA( vala.a(), 255, 1e-4 );
    vala = mean_pixel_value(im5);
    TS_ASSERT_DELTA( vala.r(), 0.0, 1e-6 );
    TS_ASSERT_DELTA( vala.g(), 0.45, 1e-6 );
    TS_ASSERT_DELTA( vala.b(), 1.05, 1e-6 );
    TS_ASSERT_DELTA( vala.a(), 0.75, 1e-6 );
    vala = mean_pixel_value(im6);
    TS_ASSERT_DELTA( vala.r(), 24, 1e-4 );
    TS_ASSERT_DELTA( vala.g(), 41, 1e-4 );
    TS_ASSERT_DELTA( vala.b(), 67.5, 1e-4 );
    TS_ASSERT_DELTA( vala.a(), 153, 1e-4 );
    vala = mean_pixel_value(im7);
    TS_ASSERT_DELTA( vala.r(), 0.25, 1e-6 );
    TS_ASSERT_DELTA( vala.g(), 0.45, 1e-6 );
    TS_ASSERT_DELTA( vala.b(), 0.8, 1e-6 );
    TS_ASSERT_DELTA( vala.a(), 0.5, 1e-6 );
    vala = mean_pixel_value(im8);
    TS_ASSERT_DELTA( vala.r(), 12, 1e-4 );
    TS_ASSERT_DELTA( vala.g(), 33.5, 1e-4 );
    TS_ASSERT_DELTA( vala.b(), 44.5, 1e-4 );
    TS_ASSERT_DELTA( vala.a(), 127.5, 1e-4 );
    vala = mean_pixel_value(im9);
    TS_ASSERT_EQUALS( vala.r(), 0 );
    TS_ASSERT_EQUALS( vala.g(), 0 );
    TS_ASSERT_EQUALS( vala.b(), 0 );
    TS_ASSERT_EQUALS( vala.a(), 0 );
    vala = mean_pixel_value(im10);
    TS_ASSERT_EQUALS( vala.r(), 0 );
    TS_ASSERT_EQUALS( vala.g(), 0 );
    TS_ASSERT_EQUALS( vala.b(), 0 );
    TS_ASSERT_EQUALS( vala.a(), 0 );
  }

  void test_weighted_mean_pixel_value()
  {
    PixelRGB<double> val;
    PixelRGBA<double> vala;
    val = weighted_mean_pixel_value(im1);
    TS_ASSERT_DELTA( val.r(), -0.25, 1e-6 );
    TS_ASSERT_DELTA( val.g(), 0.45, 1e-6 );
    TS_ASSERT_DELTA( val.b(), 1.3, 1e-6 );
    val = weighted_mean_pixel_value(im2);
    TS_ASSERT_DELTA( val.r(), 71.5, 1e-4 );
    TS_ASSERT_DELTA( val.g(), 71.5, 1e-4 );
    TS_ASSERT_DELTA( val.b(), 158.5, 1e-4 );
    vala = weighted_mean_pixel_value(im3);
    TS_ASSERT_DELTA( vala.r(), -0.25, 1e-6 );
    TS_ASSERT_DELTA( vala.g(), 0.45, 1e-6 );
    TS_ASSERT_DELTA( vala.b(), 1.3, 1e-6 );
    TS_ASSERT_DELTA( vala.a(), 1.0, 1e-6 );
    vala = weighted_mean_pixel_value(im4);
    TS_ASSERT_DELTA( vala.r(), 71.5, 1e-4 );
    TS_ASSERT_DELTA( vala.g(), 71.5, 1e-4 );
    TS_ASSERT_DELTA( vala.b(), 158.5, 1e-4 );
    TS_ASSERT_DELTA( vala.a(), 255, 1e-4 );
    vala = weighted_mean_pixel_value(im5);
    TS_ASSERT_DELTA( vala.r(), 0.0, 1e-6 );
    TS_ASSERT_DELTA( vala.g(), 0.6, 1e-6 );
    TS_ASSERT_DELTA( vala.b(), 1.4, 1e-6 );
    TS_ASSERT_DELTA( vala.a(), 1.0, 1e-6 );
    vala = weighted_mean_pixel_value(im6);
    TS_ASSERT_DELTA( vala.r(), 40, 1e-4 );
    TS_ASSERT_DELTA( vala.g(), 68.3333, 1e-4 );
    TS_ASSERT_DELTA( vala.b(), 112.5, 1e-4 );
    TS_ASSERT_DELTA( vala.a(), 255, 1e-4 );
    vala = weighted_mean_pixel_value(im7);
    TS_ASSERT_DELTA( vala.r(), 0.5, 1e-6 );
    TS_ASSERT_DELTA( vala.g(), 0.9, 1e-6 );
    TS_ASSERT_DELTA( vala.b(), 1.6, 1e-6 );
    TS_ASSERT_DELTA( vala.a(), 1.0, 1e-6 );
    vala = weighted_mean_pixel_value(im8);
    TS_ASSERT_DELTA( vala.r(), 24, 1e-4 );
    TS_ASSERT_DELTA( vala.g(), 67, 1e-4 );
    TS_ASSERT_DELTA( vala.b(), 89, 1e-4 );
    TS_ASSERT_DELTA( vala.a(), 255, 1e-4 );
    TS_ASSERT_THROWS(weighted_mean_pixel_value(im9), ArgumentErr);
    TS_ASSERT_THROWS(weighted_mean_pixel_value(im10), ArgumentErr);
  }

  void test_mean_channel_value()
  {
    TS_ASSERT_DELTA(mean_channel_value(im1), 0.5, 1e-6);
    TS_ASSERT_DELTA(mean_channel_value(im2), 100.5, 1e-4);
    TS_ASSERT_DELTA(mean_channel_value(im3), 0.5, 1e-6);
    TS_ASSERT_DELTA(mean_channel_value(im4), 100.5, 1e-4);
    TS_ASSERT_DELTA(mean_channel_value(im5), 0.5, 1e-6);
    TS_ASSERT_DELTA(mean_channel_value(im6), 44.1667, 1e-4);
    TS_ASSERT_DELTA(mean_channel_value(im7), 0.5, 1e-6);
    TS_ASSERT_DELTA(mean_channel_value(im8), 30.0, 1e-4);
    TS_ASSERT_EQUALS(mean_channel_value(im9), 0);
    TS_ASSERT_EQUALS(mean_channel_value(im10), 0);
  }

  void test_weighted_mean_channel_value()
  {
    TS_ASSERT_DELTA(weighted_mean_channel_value(im1), 0.5, 1e-6);
    TS_ASSERT_DELTA(weighted_mean_channel_value(im2), 100.5, 1e-4);
    TS_ASSERT_DELTA(weighted_mean_channel_value(im3), 0.5, 1e-6);
    TS_ASSERT_DELTA(weighted_mean_channel_value(im4), 100.5, 1e-4);
    TS_ASSERT_DELTA(weighted_mean_channel_value(im5), 0.666667, 1e-6);
    TS_ASSERT_DELTA(weighted_mean_channel_value(im6), 73.6111, 1e-4);
    TS_ASSERT_DELTA(weighted_mean_channel_value(im7), 1.0, 1e-6);
    TS_ASSERT_DELTA(weighted_mean_channel_value(im8), 60.0, 1e-4);
    TS_ASSERT_THROWS(weighted_mean_channel_value(im9), ArgumentErr);
    TS_ASSERT_THROWS(weighted_mean_channel_value(im10), ArgumentErr);
  }

  void test_stddev_channel_value()
  {
    TS_ASSERT_DELTA(stddev_channel_value(im1), 0.828654, 1e-6);
    TS_ASSERT_DELTA(stddev_channel_value(im2), 63.6468, 1e-4);
    TS_ASSERT_DELTA(stddev_channel_value(im3), 0.828654, 1e-6);
    TS_ASSERT_DELTA(stddev_channel_value(im4), 63.6468, 1e-4);
    TS_ASSERT_DELTA(stddev_channel_value(im5), 0.763035, 1e-6);
    TS_ASSERT_DELTA(stddev_channel_value(im6), 47.3288, 1e-4);
    TS_ASSERT_DELTA(stddev_channel_value(im7), 0.454606, 1e-6);
    TS_ASSERT_DELTA(stddev_channel_value(im8), 26.9938, 1e-4);
    TS_ASSERT_THROWS(stddev_channel_value(im9), ArgumentErr);
    TS_ASSERT_THROWS(stddev_channel_value(im10), ArgumentErr);
  }

  void test_median_channel_value()
  {
    ImageView<PixelRGBA<uint8> > image(4,1);
    image(0,0) = PixelRGBA<uint8>(8,7,6,5);
    image(1,0) = PixelRGBA<uint8>(4,3,2,1);
    image(2,0) = PixelRGBA<uint8>(9,8,7,0);
    image(3,0) = PixelRGBA<uint8>(12,11,10,0);
    TS_ASSERT_EQUALS( median_channel_value(image), 6 );
  }

};


