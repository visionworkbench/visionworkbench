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
  void test_homography_fitting_functor2() 
  {    
    static double A_data[] = {  
      0.0153, 0.9318, 1,
      0.6721, 0.6813, 1,
      0.3046, 0.6822, 1,
      0.8600, 0.7468, 1,
      0.5252, 0.8381, 1,
      0.7095, 0.1897, 1,
      0.6979, 0.8537, 1,
      0.4186, 0.2026, 1,
      0.8318, 0.4289, 1,
      0.5417, 0.3784, 1 };
    
    vw::MatrixProxy<double,10,3> A(A_data);
    
    static double H_data[] = {
      1.1567,  0.5916,  0.5557,
      0.2814,  1.0851,  0.0225,
      0.7388,  0.9278,  1.0000};

    vw::MatrixProxy<double,3,3> H(H_data);
    
    // and apply it to some points
    vw::Matrix<double> B = transpose(H*transpose(A));
    std::vector<Vector3> p1, p2;    
    for (unsigned i = 0; i < A.rows(); ++i) {
      p1.push_back(select_row(A,i));
      p2.push_back(select_row(B,i));
    }

    HomographyFittingFunctor hom_fit;
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

void test_affinity_fitting_functor() 
  {
    
    static double A_data[] = {  
      0.0153, 0.9318, 1,
      0.6721, 0.6813, 1,
      0.3046, 0.6822, 1,
      0.8600, 0.7468, 1,
      0.5252, 0.8381, 1,
      0.7095, 0.1897, 1,
      0.6979, 0.8537, 1,
      0.4186, 0.2026, 1,
      0.8318, 0.4289, 1,
      0.5417, 0.3784, 1 };    
    vw::MatrixProxy<double,10,3> A(A_data);
    
    static double S_data[] = {
      1.1567,  0.5916,  0.5557,
      0.2814,  1.0851,  0.0225,
      0,  0,  1.0000};

    vw::MatrixProxy<double,3,3> S(S_data);
    
    vw::Matrix<double> B = transpose(S*transpose(A));
    std::vector<Vector3> p1, p2;    
    for (unsigned i = 0; i < A.rows(); ++i) {
      p1.push_back(select_row(A,i));
      p2.push_back(select_row(B,i));
    }
    
    AffineFittingFunctor sim_fit;
    vw::Matrix<double> newS = sim_fit(p1,p2);

    TS_ASSERT_DELTA(S(0,0), newS(0,0), 0.001);
    TS_ASSERT_DELTA(S(0,1), newS(0,1), 0.001);
    TS_ASSERT_DELTA(S(0,2), newS(0,2), 0.001);
    TS_ASSERT_DELTA(S(1,0), newS(1,0), 0.001);
    TS_ASSERT_DELTA(S(1,1), newS(1,1), 0.001);
    TS_ASSERT_DELTA(S(1,2), newS(1,2), 0.001);
    TS_ASSERT_DELTA(S(2,0), newS(2,0), 0.001);
    TS_ASSERT_DELTA(S(2,1), newS(2,1), 0.001);
    TS_ASSERT_DELTA(S(2,2), newS(2,2), 0.001);
  }    

  void test_similarity_fitting() 
  {
    static double A_data[] = {  
      0.0153, 0.9318, 1,
      0.6721, 0.6813, 1,
      0.3046, 0.6822, 1,
      0.8600, 0.7468, 1,
      0.5252, 0.8381, 1,
      0.7095, 0.1897, 1,
      0.6979, 0.8537, 1,
      0.4186, 0.2026, 1,
      0.8318, 0.4289, 1,
      0.5417, 0.3784, 1 };
    
    vw::MatrixProxy<double,10,3> A(A_data);

    // Create a fictitious similarity matrix
    double s = 2.23;
    double theta = 0.2342;
    double dx = 4.2;
    double dy = 2.9;

    Matrix3x3 S;
    S(0,0) = s*cos(theta);
    S(0,1) = s*sin(theta);
    S(1,0) = s*-sin(theta);
    S(1,1) = s*cos(theta);
    S(0,2) = dx;
    S(1,2) = dy;
    S(2,2) = 1;

    // and apply it to some points
    vw::Matrix<double> B = transpose(S*transpose(A));
    std::vector<Vector3> p1, p2;    
    for (unsigned i = 0; i < A.rows(); ++i) {
      p1.push_back(select_row(A,i));
      p2.push_back(select_row(B,i));
    }

    // compute a fit given the points only, and compare it to the true
    // similarity.
    SimilarityFittingFunctor sim_fit;
    vw::Matrix<double> newS = sim_fit(p1,p2);

    TS_ASSERT_DELTA(S(0,0), newS(0,0), 1e-6);
    TS_ASSERT_DELTA(S(0,1), newS(0,1), 1e-6);
    TS_ASSERT_DELTA(S(0,2), newS(0,2), 1e-6);
    TS_ASSERT_DELTA(S(1,0), newS(1,0), 1e-6);
    TS_ASSERT_DELTA(S(1,1), newS(1,1), 1e-6);
    TS_ASSERT_DELTA(S(1,2), newS(1,2), 1e-6);
    TS_ASSERT_DELTA(S(2,0), newS(2,0), 1e-6);
    TS_ASSERT_DELTA(S(2,1), newS(2,1), 1e-6);
    TS_ASSERT_DELTA(S(2,2), newS(2,2), 1e-6);
  }    

  

}; // class TestGeometry
