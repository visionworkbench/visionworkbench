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

#include <vw/Image/Statistics.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>

using namespace std;
using namespace vw;

class TestImageStatistics : public CxxTest::TestSuite
{
public:

  void test_min_channel_value()
  {
    ImageView<PixelRGB<double> > im(2,1);
    im(0,0) = PixelRGB<double>(0.5,0.9,1.5);
    im(1,0) = PixelRGB<double>(-1,0,1);

    double out_min = min_channel_value(im);
    
    TS_ASSERT_EQUALS(out_min, -1);
  }

  void test_max_channel_value()
  {
    ImageView<PixelRGB<double> > im(2,1);
    im(0,0) = PixelRGB<double>(0.5,0.9,1.5);
    im(1,0) = PixelRGB<double>(-1,0,1);

    double out_max = max_channel_value(im);
    
    TS_ASSERT_EQUALS(out_max, 1.5);
  }

  void test_min_max_channel_values()
  {
    ImageView<PixelRGB<double> > im(2,1);
    im(0,0) = PixelRGB<double>(0.5,0.9,1.5);
    im(1,0) = PixelRGB<double>(-1,0,1);

    double out_min, out_max;
    min_max_channel_values(im, out_min, out_max);
    
    TS_ASSERT_EQUALS(out_max, 1.5);
    TS_ASSERT_EQUALS(out_min, -1);
  }

  void test_mean_channel_value()
  {
    ImageView<PixelRGB<double> > im(2,1);
    im(0,0) = PixelRGB<double>(0.5,0.9,1.6);
    im(1,0) = PixelRGB<double>(-1,0,1);

    double out_mean = mean_channel_value(im);

    TS_ASSERT_EQUALS(out_mean, 0.5);
  }

};


