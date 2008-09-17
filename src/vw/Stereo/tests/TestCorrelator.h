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

// TestCorrelator.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/UtilityViews.h>
#include <vw/Stereo/Correlator.h>
#include <vw/Image/Transform.h>

using namespace std;
using namespace vw;
using namespace vw::stereo;

class TestCorrelator : public CxxTest::TestSuite
{
public:
//   void test_image_sum_functions() {
//     ImageView<uint8> image = ConstantView<uint8>(1, 100, 100);
//     int kern_width = 10;
//     int kern_height = 10;

//     // Computes the complete sum at every pixel
//     for (int j = kern_height/2; j < image.rows()-kern_height/2; ++j)
//       for (int i = kern_width/2; i < image.cols()-kern_width/2; ++i) {
//         int16 sum = compute_sum(image, BBox2i(i-kern_width/2,j-kern_height/2,kern_width,kern_height));
//         TS_ASSERT_EQUALS(sum,kern_width*kern_height);
//       }

//     // Computes the complete sum at the beginning of the scanline, and
//     // updates as it moves across the image.
//     for (int j = kern_height/2; j < image.rows()-kern_height/2; ++j) {
//       bool need_to_compute_full_sum = true;
//       for (int i = kern_width/2; i < image.cols()-kern_width/2; ++i) {
//         int16 sum;
//         if (need_to_compute_full_sum) {
//           sum = compute_sum(image, BBox2i(i-kern_width/2,j-kern_height/2,kern_width,kern_height));
//           need_to_compute_full_sum = false;
//         } else {
//           sum = update_sum(image, sum, BBox2i(i-kern_width/2,j-kern_height/2,kern_width,kern_height));
//         }
//         TS_ASSERT_EQUALS(sum,kern_width*kern_height);
//       }
//     }
//   }


  void test_basic_correlation() {
    BBox2i search_range(-5,-5,10,10);
    Correlator corr(search_range, 5,5,1.5,1);

    ImageView<uint8> image1(50,50);
    for (int j = 0; j < image1.rows(); ++j) {
      for (int i = 0; i < image1.cols(); ++i) {
        image1(i,j) = j*image1.cols()+i;
      }
    }
    ImageView<uint8> image2 = transform(image1, TranslateTransform(3,3), ZeroEdgeExtension(), NearestPixelInterpolation());

    ImageView<PixelDisparity<float> > disparity_map;
    disparity_map = corr(image1,image2);

    for (int j = 0; j < disparity_map.rows(); ++j)
      for (int i = 0; i < disparity_map.cols(); ++i)
        TS_TRACE(stringify(disparity_map(i,j)));

  }


};
