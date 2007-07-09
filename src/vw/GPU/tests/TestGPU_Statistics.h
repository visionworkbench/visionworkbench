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

#include <vw/GPU.h>

using namespace std;
using namespace vw;
using namespace GPU;

#define DELTA_PRECISION 1e-5

class TestImageStatistics : public CxxTest::TestSuite
{
public:
	
	void test_min_channel_value()
	{
		gpu_init();
		{
			GPUImage<PixelRGB<float> > im(2,1);
			im.pixel(0,0) = PixelRGB<float>(0.5,0.9,1.5);
			im.pixel(1,0) = PixelRGB<float>(-1,0,1);
			float out_min = min_channel_value(im);
			TS_ASSERT_DELTA(out_min, -1, DELTA_PRECISION);
		}
		{
			GPUImage<float> im(5, 5);
			int total = im.cols() * im.rows();
			int i = total-1;
			for(int x=0; x < im.cols(); x++) {
				for(int y=0; y< im.rows(); y++) {
					im.pixel(x, y) = 0.2 + 0.5*(i-- / (float) total);
				}
			}
			float out_min = min_channel_value(im);
			TS_ASSERT_DELTA(out_min, 0.2, DELTA_PRECISION);
		}
	}

void test_max_channel_value()
{
	gpu_init();
	{
		GPUImage<PixelRGB<float> > im(2,1);
		im.pixel(0,0) = PixelRGB<float>(0.5,0.9,1.5);
		im.pixel(1,0) = PixelRGB<float>(-1,0,1);
		float out_max = max_channel_value(im);
		TS_ASSERT_DELTA(out_max, 1.5, DELTA_PRECISION);
    }
	{
		GPUImage<float> im(5, 5);
		int total = im.cols() * im.rows();
		int i = total-1;
		for(int x=0; x < im.cols(); x++) {
			for(int y=0; y< im.rows(); y++) {
				im.pixel(x, y) = 0.7 - 0.5*(i-- / (float) total);
			}
		}
		float out_min = max_channel_value(im);
		TS_ASSERT_DELTA(out_min, 0.7, DELTA_PRECISION);
	}
	
}

void test_min_max_channel_values()
{
	gpu_init();
    GPUImage<PixelRGB<float> > im(2,1);
    im.pixel(0,0) = PixelRGB<float>(0.5,0.9,1.5);
    im.pixel(1,0) = PixelRGB<float>(-1,0,1);
	
    float out_min, out_max;
    min_max_channel_values(im, out_min, out_max);
    
    TS_ASSERT_DELTA(out_max, 1.5, DELTA_PRECISION);
    TS_ASSERT_DELTA(out_min, -1, DELTA_PRECISION);
}

void test_mean_channel_value()
{
	gpu_init();
	{
		GPUImage<PixelRGB<float> > im(2,1);
		im.pixel(0,0) = PixelRGB<float>(0.5,0.9,1.6);
		im.pixel(1,0) = PixelRGB<float>(-1,0,1);
		float out_mean = mean_channel_value(im);
		TS_ASSERT_DELTA(out_mean, 0.5, DELTA_PRECISION);
	}
	{
		GPUImage<float> im(5, 5);
		int total = 0;
		int i = 0;
		for(int x=0; x < im.cols(); x++) {
			for(int y=0; y< im.rows(); y++) {
				im.pixel(x, y) = i;
				total += i;
				i++;
			}
		}
		float out_mean = mean_channel_value(im);
		TS_ASSERT_DELTA(out_mean, total / (float) i, DELTA_PRECISION);
	}
}

};


