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

// TestLinearAlgebra.h
#include <cxxtest/TestSuite.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Geometry.h>

using namespace vw;
using namespace vw::math;

class TestGeometry : public CxxTest::TestSuite
{
public:

  void test_homography_fitting_functor() 
  {
    
    static double A_data[] = {  
      0.0153, 0.9318, 0.8462, 0.6721, 0.6813, 0.5028, 0.3046, 0.6822, 0.1509, 0.8600,
      0.7468, 0.4660, 0.5252, 0.8381, 0.3795, 0.7095, 0.1897, 0.3028, 0.6979, 0.8537,
      0.4451, 0.4186, 0.2026, 0.0196, 0.8318, 0.4289, 0.1934, 0.5417, 0.3784, 0.5936 };
    
    vw::MatrixProxy<double,3,10> A(A_data);
    
    static double B_data[] = {  
      0.4966, 0.6449, 0.3420, 0.5341, 0.8385, 0.7027, 0.6946, 0.9568, 0.1730, 0.2523,
      0.8998, 0.8180, 0.2897, 0.7271, 0.5681, 0.5466, 0.6213, 0.5226, 0.9797, 0.8757,
      0.8216, 0.6602, 0.3412, 0.3093, 0.3704, 0.4449, 0.7948, 0.8801, 0.2714, 0.7373 };
    
    vw::MatrixProxy<double,3,10> B(B_data);

    std::vector<Vector3> p1, p2;
    
    for (int i = 0; i < A.cols(); ++i) {
      p1.push_back(Vector3(A(0,i), A(1,i), A(2,i)));
      p2.push_back(Vector3(B(0,i), B(1,i), B(2,i)));
    }

    HomographyFittingFunctor hom_fit(true);
    vw::Matrix<double> H = hom_fit(p1,p2);

    TS_ASSERT_DELTA(H(0,0), 0.6254, 0.001);
    TS_ASSERT_DELTA(H(0,1), 0.0385, 0.001);
    TS_ASSERT_DELTA(H(0,2), 1.1318, 0.001);
    TS_ASSERT_DELTA(H(1,0), -0.1214, 0.001);
    TS_ASSERT_DELTA(H(1,1), 1.4143, 0.001);
    TS_ASSERT_DELTA(H(1,2), 0.7983, 0.001);
    TS_ASSERT_DELTA(H(2,0), 0.2912, 0.001);
    TS_ASSERT_DELTA(H(2,1), 0.4783, 0.001);
    TS_ASSERT_DELTA(H(2,2), 1.0, 0.001);
  }    

  void test_homography_fitting_functor2() 
  {
    
    static double A_data[] = {  
      0.0153, 0.9318, 0.8462, 0.6721, 0.6813, 0.5028, 0.3046, 0.6822, 0.1509, 0.8600,
      0.7468, 0.4660, 0.5252, 0.8381, 0.3795, 0.7095, 0.1897, 0.3028, 0.6979, 0.8537,
      0.4451, 0.4186, 0.2026, 0.0196, 0.8318, 0.4289, 0.1934, 0.5417, 0.3784, 0.5936 };
    
    vw::MatrixProxy<double,3,10> A(A_data);
    
    static double H_data[] = {
      1.1567,  0.5916,  0.5557,
      0.2814,  1.0851,  0.0225,
      0.7388,  0.9278,  1.0000};

    vw::MatrixProxy<double,3,3> H(H_data);
    
    vw::Matrix<double> B = H*A;
    
    std::vector<Vector3> p1, p2;    
    for (int i = 0; i < A.cols(); ++i) {
      p1.push_back(Vector3(A(0,i), A(1,i), A(2,i)));
      p2.push_back(Vector3(B(0,i), B(1,i), B(2,i)));
    }

    HomographyFittingFunctor hom_fit(true);
    vw::Matrix<double> newH = hom_fit(p1,p2);

    TS_ASSERT_DELTA(H(0,0), newH(0,0), 0.001);
    TS_ASSERT_DELTA(H(0,1), newH(0,1), 0.001);
    TS_ASSERT_DELTA(H(0,2), newH(0,2), 0.001);
    TS_ASSERT_DELTA(H(1,0), newH(1,0), 0.001);
    TS_ASSERT_DELTA(H(1,1), newH(1,1), 0.001);
    TS_ASSERT_DELTA(H(1,2), newH(1,2), 0.001);
    TS_ASSERT_DELTA(H(2,0), newH(2,0), 0.001);
    TS_ASSERT_DELTA(H(2,1), newH(2,1), 0.001);
    TS_ASSERT_DELTA(H(2,2), newH(2,2), 0.001);
    
    HomographyFittingFunctor hom_fit2;
    newH = hom_fit2(p1,p2);

    TS_ASSERT_DELTA(H(0,0), newH(0,0), 0.001);
    TS_ASSERT_DELTA(H(0,1), newH(0,1), 0.001);
    TS_ASSERT_DELTA(H(0,2), newH(0,2), 0.001);
    TS_ASSERT_DELTA(H(1,0), newH(1,0), 0.001);
    TS_ASSERT_DELTA(H(1,1), newH(1,1), 0.001);
    TS_ASSERT_DELTA(H(1,2), newH(1,2), 0.001);
    TS_ASSERT_DELTA(H(2,0), newH(2,0), 0.001);
    TS_ASSERT_DELTA(H(2,1), newH(2,1), 0.001);
    TS_ASSERT_DELTA(H(2,2), newH(2,2), 0.001);
  }    

  void test_similarity_fitting_functor() 
  {
    
    static double A_data[] = {  
      0.0153, 0.9318, 0.8462, 0.6721, 0.6813, 0.5028, 0.3046, 0.6822, 0.1509, 0.8600,
      0.7468, 0.4660, 0.5252, 0.8381, 0.3795, 0.7095, 0.1897, 0.3028, 0.6979, 0.8537,
      0.4451, 0.4186, 0.2026, 0.0196, 0.8318, 0.4289, 0.1934, 0.5417, 0.3784, 0.5936 };
    
    vw::MatrixProxy<double,3,10> A(A_data);
    
    static double H_data[] = {
      1.1567,  0.5916,  0.5557,
      0.2814,  1.0851,  0.0225,
      0,  0,  1.0000};

    vw::MatrixProxy<double,3,3> H(H_data);
    
    vw::Matrix<double> B = H*A;
    
    std::vector<Vector3> p1, p2;    
    for (int i = 0; i < A.cols(); ++i) {
      p1.push_back(Vector3(A(0,i), A(1,i), A(2,i)));
      p2.push_back(Vector3(B(0,i), B(1,i), B(2,i)));
    }

    SimilarityFittingFunctor sim_fit(true);
    vw::Matrix<double> newH = sim_fit(p1,p2);


    TS_ASSERT_DELTA(H(0,0), newH(0,0), 0.001);
    TS_ASSERT_DELTA(H(0,1), newH(0,1), 0.001);
    TS_ASSERT_DELTA(H(0,2), newH(0,2), 0.001);
    TS_ASSERT_DELTA(H(1,0), newH(1,0), 0.001);
    TS_ASSERT_DELTA(H(1,1), newH(1,1), 0.001);
    TS_ASSERT_DELTA(H(1,2), newH(1,2), 0.001);
    TS_ASSERT_DELTA(H(2,0), newH(2,0), 0.001);
    TS_ASSERT_DELTA(H(2,1), newH(2,1), 0.001);
    TS_ASSERT_DELTA(H(2,2), newH(2,2), 0.001);
    
    SimilarityFittingFunctor sim_fit2;
    newH = sim_fit2(p1,p2);

    TS_ASSERT_DELTA(H(0,0), newH(0,0), 0.001);
    TS_ASSERT_DELTA(H(0,1), newH(0,1), 0.001);
    TS_ASSERT_DELTA(H(0,2), newH(0,2), 0.001);
    TS_ASSERT_DELTA(H(1,0), newH(1,0), 0.001);
    TS_ASSERT_DELTA(H(1,1), newH(1,1), 0.001);
    TS_ASSERT_DELTA(H(1,2), newH(1,2), 0.001);
    TS_ASSERT_DELTA(H(2,0), newH(2,0), 0.001);
    TS_ASSERT_DELTA(H(2,1), newH(2,1), 0.001);
    TS_ASSERT_DELTA(H(2,2), newH(2,2), 0.001);
  }    

}; // class TestGeometry
