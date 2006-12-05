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
#include <vw/Math/LinearAlgebra.h>

using namespace vw;
using namespace vw::math;

class TestLinearAlgebra : public CxxTest::TestSuite
{
public:

  void test_linear_least_squares_static()
  {
    Matrix<float,4,2> A;
    A(0,0) = 4; A(0,1) = 17; 
    A(1,0) = 32; A(1,1) = 33; 
    A(2,0) = 63; A(2,1) = 3; 
    A(3,0) = 2; A(3,1) = 73; 

    Vector<float,4> b;
    b(0) = 5;
    b(1) = 95;
    b(2) = 2;
    b(3) = 1;

    Vector<float,2> x;

    // Find the least squares solution to the overconstrained problem Ax=b;
    x = least_squares(A,b);

    TS_ASSERT_DELTA(x[0], 0.52626, 0.00001);
    TS_ASSERT_DELTA(x[1], 0.37689, 0.00001);
  }    

  void test_linear_least_squares_underdetermined()
  {
    Matrix<float,2,4> A;
    A(0,0) = 4; A(1,0) = 17; 
    A(0,1) = 32; A(1,1) = 33; 
    A(0,2) = 63; A(1,2) = 3; 
    A(0,3) = 2; A(1,3) = 73; 

    Vector<float,2> b;
    b(0) = 5;
    b(1) = 95;

    Vector<float,4> x;

    // Find the least squares solution to the overconstrained problem Ax=b;
    x = least_squares(A,b);

    Vector<float,2> b_prime = A*x;
 
    TS_ASSERT_DELTA(b_prime[0], 5, 1e-4);
    TS_ASSERT_DELTA(b_prime[1], 95, 1e-4);
  }    

  void test_linear_least_squares_dynamic()
  {
    Matrix<float> A(4,2);
    A(0,0) = 4; A(0,1) = 17; 
    A(1,0) = 32; A(1,1) = 33; 
    A(2,0) = 63; A(2,1) = 3; 
    A(3,0) = 2; A(3,1) = 73; 

    Vector<float> b(4);
    b(0) = 5;
    b(1) = 95;
    b(2) = 2;
    b(3) = 1;

    Vector<float> x;

    // Find the least squares solution to the overconstrained problem Ax=b;
    x = least_squares(A,b);

    TS_ASSERT_DELTA(x[0], 0.52626, 0.00001);
    TS_ASSERT_DELTA(x[1], 0.37689, 0.00001);
  }    



  // Test the svd(A,s) routine with a float matrix
  void test_svd_float()
  {
    Matrix<float> A(4,4);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; A(0,3) = 98;
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76; A(1,3) = 66;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76; A(2,3) = 662;
    A(3,0) = 25;  A(3,1) = 62; A(3,2) = 323; A(3,3) = 23;

    Vector<float> s;

    // Find the least squares solution to the overconstrained problem Ax=b;
    svd(A,s);

    TS_ASSERT_DELTA(s[0], 744.6188, 0.0001);
    TS_ASSERT_DELTA(s[1], 336.4697, 0.0001);
    TS_ASSERT_DELTA(s[2], 262.4880, 0.0001);
    TS_ASSERT_DELTA(s[3], 5.8381, 0.0001);

    A = Matrix<float>(3,4);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; A(0,3) = 98;
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76; A(1,3) = 66;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76; A(2,3) = 662;

    // Find the least squares solution to the overconstrained problem Ax=b;
    svd(A,s);

    TS_ASSERT_DELTA(s[0], 739.8530, 0.0001);
    TS_ASSERT_DELTA(s[1], 282.1995, 0.0001);
    TS_ASSERT_DELTA(s[2], 16.0320, 0.0001);

    A = Matrix<float>(4,3);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; 
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76;
    A(3,0) = 25;  A(3,1) = 62; A(3,2) = 323; 

    // Find the least squares solution to the overconstrained problem Ax=b;
    svd(A,s);

    TS_ASSERT_DELTA(s[0], 444.7863, 0.0001);
    TS_ASSERT_DELTA(s[1], 292.8446, 0.0001);
    TS_ASSERT_DELTA(s[2], 16.6494, 0.0001);
  }    

  // Test the svd(A,s) routine with a double matrix
  void test_svd_double()
  {
    Matrix<double> A(4,4);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; A(0,3) = 98;
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76; A(1,3) = 66;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76; A(2,3) = 662;
    A(3,0) = 25;  A(3,1) = 62; A(3,2) = 323; A(3,3) = 23;

    Vector<double> s;

    // Find the least squares solution to the overconstrained problem Ax=b;
    svd(A,s);

    TS_ASSERT_DELTA(s[0], 744.6188, 0.0001);
    TS_ASSERT_DELTA(s[1], 336.4697, 0.0001);
    TS_ASSERT_DELTA(s[2], 262.4880, 0.0001);
    TS_ASSERT_DELTA(s[3], 5.8381, 0.0001);

    A = Matrix<double>(3,4);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; A(0,3) = 98;
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76; A(1,3) = 66;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76; A(2,3) = 662;

    // Find the least squares solution to the overconstrained problem Ax=b;
    svd(A,s);

    TS_ASSERT_DELTA(s[0], 739.8530, 0.0001);
    TS_ASSERT_DELTA(s[1], 282.1995, 0.0001);
    TS_ASSERT_DELTA(s[2], 16.0320, 0.0001);

    A = Matrix<double>(4,3);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; 
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76;
    A(3,0) = 25;  A(3,1) = 62; A(3,2) = 323; 

    // Find the least squares solution to the overconstrained problem Ax=b;
    svd(A,s);

    TS_ASSERT_DELTA(s[0], 444.7863, 0.0001);
    TS_ASSERT_DELTA(s[1], 292.8446, 0.0001);
    TS_ASSERT_DELTA(s[2], 16.6494, 0.0001);
  }    

  // Test the svd(A,s) routine with a mixture of float and double matrix/vectors
  void test_svd_double_float()
  {
    Matrix<double> A(4,4);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; A(0,3) = 98;
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76; A(1,3) = 66;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76; A(2,3) = 662;
    A(3,0) = 25;  A(3,1) = 62; A(3,2) = 323; A(3,3) = 23;

    Vector<float> s;

    // Find the least squares solution to the overconstrained problem Ax=b;
    svd(A,s);

    TS_ASSERT_DELTA(s[0], 744.6188, 0.0001);
    TS_ASSERT_DELTA(s[1], 336.4697, 0.0001);
    TS_ASSERT_DELTA(s[2], 262.4880, 0.0001);
    TS_ASSERT_DELTA(s[3], 5.8381, 0.0001);
  }    

  // Test the svd(A,s) routine with a mixture of float and double matrix/vectors
  void test_svd_float_double()
  {
    Matrix<float> A(4,4);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; A(0,3) = 98;
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76; A(1,3) = 66;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76; A(2,3) = 662;
    A(3,0) = 25;  A(3,1) = 62; A(3,2) = 323; A(3,3) = 23;

    Vector<double> s;

    // Find the least squares solution to the overconstrained problem Ax=b;
    svd(A,s);

    TS_ASSERT_DELTA(s[0], 744.6188, 0.0001);
    TS_ASSERT_DELTA(s[1], 336.4697, 0.0001);
    TS_ASSERT_DELTA(s[2], 262.4880, 0.0001);
    TS_ASSERT_DELTA(s[3], 5.8381, 0.0001);
  }    



  // Test the svd(A,U,s,VT) routine with a float matrix
  void test_svd_U_s_VT_float()
  {
    Matrix<float> A(4,4);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; A(0,3) = 98;
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76; A(1,3) = 66;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76; A(2,3) = 662;
    A(3,0) = 25;  A(3,1) = 62; A(3,2) = 323; A(3,3) = 23;

    Matrix<float> U, VT;
    Vector<float> s;

    // Find the least squares solution to the overconstrained problem Ax=b;
    svd(A,U,s,VT);

    TS_ASSERT_DELTA(U(0,0), -0.1357, 0.0001);
    TS_ASSERT_DELTA(U(0,1), 0.0133, 0.0001);
    TS_ASSERT_DELTA(U(1,0),  -0.2857, 0.0001);
    TS_ASSERT_DELTA(U(2,3),  0.1481, 0.0001);

    TS_ASSERT_DELTA(VT(0,0),  -0.4294, 0.0001);
    TS_ASSERT_DELTA(VT(0,1), -0.0442, 0.0001);
    TS_ASSERT_DELTA(VT(1,0), -0.3203, 0.0001);
    TS_ASSERT_DELTA(VT(2,3), 0.3173, 0.0001);

    TS_ASSERT_EQUALS(U.rows(), A.rows());
    TS_ASSERT_EQUALS(U.cols(), std::min(A.rows(), A.cols()));
    TS_ASSERT_EQUALS(VT.rows(), std::min(A.rows(), A.cols()));
    TS_ASSERT_EQUALS(VT.cols(), A.cols());

    A = Matrix<float>(3,4);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; A(0,3) = 98;
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76; A(1,3) = 66;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76; A(2,3) = 662;

    // Find the least squares solution to the overconstrained problem Ax=b;
    svd(A,U,s,VT);

    TS_ASSERT_DELTA(U(0,0),  -0.1361, 0.0001);
    TS_ASSERT_DELTA(U(0,1), 0.0682, 0.0001);
    TS_ASSERT_DELTA(U(1,0), -0.2831, 0.0001);
    TS_ASSERT_DELTA(U(1,2), 0.0271, 0.0001);

    TS_ASSERT_DELTA(VT(0,0), -0.4296, 0.0001);
    TS_ASSERT_DELTA(VT(0,1), -0.0343, 0.0001);
    TS_ASSERT_DELTA(VT(1,0), -0.8764, 0.0001);
    TS_ASSERT_DELTA(VT(1,3), 0.4472, 0.0001);

    TS_ASSERT_EQUALS(U.rows(), A.rows());
    TS_ASSERT_EQUALS(U.cols(), std::min(A.rows(), A.cols()));
    TS_ASSERT_EQUALS(VT.rows(), std::min(A.rows(), A.cols()));
    TS_ASSERT_EQUALS(VT.cols(), A.cols());

    A = Matrix<float>(4,3);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; 
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76;
    A(3,0) = 25;  A(3,1) = 62; A(3,2) = 323; 

    // Find the least squares solution to the overconstrained problem Ax=b;
    svd(A,U,s,VT);

    TS_ASSERT_DELTA(U(0,0), -0.0740, 0.0001);
    TS_ASSERT_DELTA(U(0,1), 0.0261, 0.0001);
    TS_ASSERT_DELTA(U(1,0), -0.7041 , 0.0001);
    TS_ASSERT_DELTA(U(1,2), 0.5545, 0.0001);

    TS_ASSERT_DELTA(VT(0,0), -0.8294, 0.0001);
    TS_ASSERT_DELTA(VT(0,1), -0.0988, 0.0001);
    TS_ASSERT_DELTA(VT(1,0), -0.5586, 0.0001);
    TS_ASSERT_DELTA(VT(1,2), 0.8129, 0.0001);
    
    TS_ASSERT_EQUALS(U.rows(), A.rows());
    TS_ASSERT_EQUALS(U.cols(), std::min(A.rows(), A.cols()));
    TS_ASSERT_EQUALS(VT.rows(), std::min(A.rows(), A.cols()));
    TS_ASSERT_EQUALS(VT.cols(), A.cols());

  } 


  // Test the svd(A,U,s,VT) routine with a float matrix
  void test_pseudoinverse()
  {
    Matrix<float> A(4,4);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; A(0,3) = 98;
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76; A(1,3) = 66;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76; A(2,3) = 662;
    A(3,0) = 25;  A(3,1) = 62; A(3,2) = 323; A(3,3) = 23;

    Matrix<float> pinv = pseudoinverse(A);

    TS_ASSERT_DELTA(pinv(0,0), -0.0075, 0.0001);
    TS_ASSERT_DELTA(pinv(0,1), 0.0030, 0.0001);
    TS_ASSERT_DELTA(pinv(1,0),  -0.1657, 0.0001);
    TS_ASSERT_DELTA(pinv(2,3),  0.0015, 0.0001);

    A = Matrix<float>(3,4);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; A(0,3) = 98;
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76; A(1,3) = 66;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76; A(2,3) = 662;

    pinv = pseudoinverse(A);

    TS_ASSERT_DELTA(pinv(0,0), -0.0134, 0.0001);
    TS_ASSERT_DELTA(pinv(0,1), 0.0028, 0.0001);
    TS_ASSERT_DELTA(pinv(1,0), -0.0110, 0.0001);
    TS_ASSERT_DELTA(pinv(1,2), 0.0017, 0.0001);

    TS_ASSERT_EQUALS(pinv.rows(), A.cols());
    TS_ASSERT_EQUALS(pinv.cols(), A.rows());

    A = Matrix<float>(4,3);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; 
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76;
    A(3,0) = 25;  A(3,1) = 62; A(3,2) = 323; 

    pinv = pseudoinverse(A);

    TS_ASSERT_EQUALS(pinv.rows(), A.cols());
    TS_ASSERT_EQUALS(pinv.cols(), A.rows());

  } 


  // Test the svd(A,s) routine with a double matrix
  void test_eigen()
  {
    Matrix<double> A(4,4);
    A(0,0) = 23;  A(0,1) = 1; A(0,2) = 25; A(0,3) = 98;
    A(1,0) = 327;  A(1,1) = 2; A(1,2) = 76; A(1,3) = 66;
    A(2,0) = 234;  A(2,1) = 26; A(2,2) = 76; A(2,3) = 662;
    A(3,0) = 25;  A(3,1) = 62; A(3,2) = 323; A(3,3) = 23;

    Vector<std::complex<double> > e;
    Matrix<std::complex<double> > V;
    
    // Find the least squares solution to the overconstrained problem Ax=b;
    eigen(A,e);

    TS_ASSERT_DELTA(e(0).real(), 554.90, 0.01);
    TS_ASSERT_DELTA(e(0).imag(), 0, 0.01);
    TS_ASSERT_DELTA(e(1).real(), -16.48, 0.01);
    TS_ASSERT_DELTA(e(1).imag(), 38.30, 0.01);
    TS_ASSERT_DELTA(e(2).real(), -16.48, 0.01);
    TS_ASSERT_DELTA(e(2).imag(), -38.30, 0.01);
    TS_ASSERT_DELTA(e(3).real(), -397.94, 0.01);
    TS_ASSERT_DELTA(e(3).imag(), 0, 0.01);

    // Find the least squares solution to the overconstrained problem Ax=b;
    eigen(A,e,V);

    TS_ASSERT_DELTA(V(0,0).real(), 0.134607, 0.0001);
    TS_ASSERT_DELTA(V(0,0).imag(), 0, 0.0001);
    TS_ASSERT_DELTA(V(0,1).real(), -0.0117482, 0.0001);
    TS_ASSERT_DELTA(V(0,1).imag(), 0.126045, 0.0001);
    TS_ASSERT_DELTA(V(1,0).real(), 0.252412, 0.0001);
    TS_ASSERT_DELTA(V(1,0).imag(), 0, 0.0001);
  }


}; // class TestVector
