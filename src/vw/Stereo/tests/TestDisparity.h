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

// TestDisparity.h
#include <cxxtest/TestSuite.h>

#include <vw/Math.h>
#include <vw/Stereo.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Transform.h>

using namespace std;
using namespace vw;
using namespace vw::stereo;

class TestDisparity : public CxxTest::TestSuite
{
 public:
  void test_disparity_transform_1() {

    vw::Matrix<double,3,3> align_matrix;
    align_matrix.set_identity();
    // Adding just a translation.
    align_matrix(0,2) = 45;
    align_matrix(1,2) = -30;
    std::cout << "\nMatrix: " << align_matrix;

    // Building disparity map
    ImageView<PixelDisparity<float> > map(5,5);
    for (unsigned i = 0; i < 5; i++)
      for (unsigned j = 0; j < 5; j++) {
	map(i,j)[0] = i*5+j;
	map(i,j)[1] = j*7+i;
	map(i,j)[2] = 0; // non-missing pixel
      }

    // Applying the inverse of the align matrix
    ImageViewRef<PixelDisparity<float> > result;
    result = disparity::transform_disparities(map, HomographyTransform(align_matrix)); 
    
    // Comparing results
    for (unsigned i = 0; i < 5; i++)
      for (unsigned j = 0; j < 5; j++) {
	Vector3 t_disparity(result(i,j)[0],
			    result(i,j)[1], 1);
	Vector3 location(i,j,1);
	Vector3 check = t_disparity + (align_matrix*location - location);
	// Disparity map maps the relative change between
	// images. Relative change from a transform is
	// (align_matrix*loc - loc)
	
	TS_ASSERT_DELTA( check[0], map(i,j)[0], .1);
	TS_ASSERT_DELTA( check[1], map(i,j)[1], .1);
      }
  }

  void test_disparity_transform_2() {

    vw::Matrix<double,3,3> align_matrix;
    align_matrix.set_identity();
    align_matrix(0,0) = 1.00679;
    align_matrix(0,1) = -0.0125401;
    align_matrix(0,2) = 116.812;
    align_matrix(1,0) = 0.00788373;
    align_matrix(1,1) = .996033;
    align_matrix(1,2) = -1.93039; //Homography affine
     std::cout << "\nMatrix: " << align_matrix;

    // Building disparity map
    ImageView<PixelDisparity<float> > map(5,5);
    for (unsigned i = 0; i < 5; i++)
      for (unsigned j = 0; j < 5; j++) {
	map(i,j)[0] = i*5+j;
	map(i,j)[1] = j*7+i;
	map(i,j)[2] = 0; // non-missing pixel
      }

    // Applying the inverse of the align matrix
    ImageViewRef<PixelDisparity<float> > result;
    result = disparity::transform_disparities(map, HomographyTransform(align_matrix)); 
    
    // Comparing results
    for (unsigned i = 0; i < 5; i++)
      for (unsigned j = 0; j < 5; j++) {
	Vector3 t_disparity(result(i,j)[0],
			    result(i,j)[1], 1);
	Vector3 location(i,j,1);
	Vector3 check = t_disparity + (align_matrix*location-location);
	TS_ASSERT_DELTA( check[0], map(i,j)[0], .1 );
	TS_ASSERT_DELTA( check[1], map(i,j)[1], .1 );
      }
  }

};
