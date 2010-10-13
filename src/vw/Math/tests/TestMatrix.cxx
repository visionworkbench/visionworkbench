// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestMatrix.h
#include <gtest/gtest.h>
#include <vw/Math/Matrix.h>

using namespace vw;

TEST(Matrix, Static) {
  // Default constructor
  Matrix<float,2,3> m1;
  ASSERT_EQ( 2u, m1.rows());
  ASSERT_EQ( 3u, m1.cols());
  EXPECT_EQ( 0, m1(0,0));
  EXPECT_EQ( 0, m1(0,1));
  EXPECT_EQ( 0, m1(0,2));
  EXPECT_EQ( 0, m1(1,0));
  EXPECT_EQ( 0, m1(1,1));
  EXPECT_EQ( 0, m1(1,2));

  // Data pointer constructor
  float data[4] = {1,2,3,4};
  Matrix2x2f m2(data);
  ASSERT_EQ( 2u, m2.rows());
  ASSERT_EQ( 2u, m2.cols());
  EXPECT_EQ( 1, m2(0, 0) );
  EXPECT_EQ( 2, m2(0, 1) );
  EXPECT_EQ( 3, m2(1, 0) );
  EXPECT_EQ( 4, m2(1, 1) );
  EXPECT_EQ( 1, m2[0][0] );
  EXPECT_EQ( 2, m2[0][1] );
  EXPECT_EQ( 3, m2[1][0] );
  EXPECT_EQ( 4, m2[1][1] );

  // Copy constructor
  Matrix2x2f m3(m2);
  ASSERT_EQ( 2u, m3.rows());
  ASSERT_EQ( 2u, m3.cols());
  EXPECT_EQ( 1, m3(0, 0) );
  EXPECT_EQ( 2, m3(0, 1) );
  EXPECT_EQ( 3, m3(1, 0) );
  EXPECT_EQ( 4, m3(1, 1) );

  // Element value constructor
  Matrix2x2f m4(1,2,3,4);
  ASSERT_EQ( 2u, m4.rows());
  ASSERT_EQ( 2u, m4.cols());
  EXPECT_EQ( 1, m4(0, 0) );
  EXPECT_EQ( 2, m4(0, 1) );
  EXPECT_EQ( 3, m4(1, 0) );
  EXPECT_EQ( 4, m4(1, 1) );

  // set_size()
  EXPECT_THROW(m1.set_size(3,3), ArgumentErr);
  EXPECT_NO_THROW(m1.set_size(2,3));

  // set_identity()
  EXPECT_NO_THROW(m2.set_identity());
  EXPECT_EQ( 1, m2(0,0) );
  EXPECT_EQ( 0, m2(0,1) );
  EXPECT_EQ( 0, m2(1,0) );
  EXPECT_EQ( 1, m2(1,1) );
  EXPECT_THROW(m2.set_identity(3), ArgumentErr);
  EXPECT_THROW(m1.set_identity(),  LogicErr);

  // Iterators
  EXPECT_EQ(&(*(m1.begin())),&(m1(0,0)));
  EXPECT_EQ(&(*(m1.begin()+1)),&(m1(0,1)));
  EXPECT_EQ(m1.end(),m1.begin()+6);
}

TEST(Matrix, Dynamic) {
  // Default constructor
  Matrix<float> m0;
  ASSERT_EQ( 0u, m0.rows() );
  ASSERT_EQ( 0u, m0.cols() );

  // Size constructor
  Matrix<float> m1(2,3);
  ASSERT_EQ( 2u , m1.rows());
  ASSERT_EQ( 3u , m1.cols());
  EXPECT_EQ( 0 , m1(0, 0) );
  EXPECT_EQ( 0 , m1(0, 1) );
  EXPECT_EQ( 0 , m1(0, 2) );
  EXPECT_EQ( 0 , m1(1, 0) );
  EXPECT_EQ( 0 , m1(1, 1) );
  EXPECT_EQ( 0 , m1(1, 2) );

  // Data pointer constructor
  float data[4] = {1,2,3,4};
  Matrix<float> m2(2,2,data);
  ASSERT_EQ( 2u, m2.rows() );
  ASSERT_EQ( 2u, m2.cols() );
  EXPECT_EQ( 1, m2(0, 0) );
  EXPECT_EQ( 2, m2(0, 1) );
  EXPECT_EQ( 3, m2(1, 0) );
  EXPECT_EQ( 4, m2(1, 1) );
  EXPECT_EQ( 1, m2[0][0] );
  EXPECT_EQ( 2, m2[0][1] );
  EXPECT_EQ( 3, m2[1][0] );
  EXPECT_EQ( 4, m2[1][1] );

  // Copy constructor
  Matrix<float> m3(m2);
  ASSERT_EQ( 2u, m3.rows());
  ASSERT_EQ( 2u, m3.cols());
  EXPECT_EQ( 1, m3(0, 0) );
  EXPECT_EQ( 2, m3(0, 1) );
  EXPECT_EQ( 3, m3(1, 0) );
  EXPECT_EQ( 4, m3(1, 1) );

  // set_size()
  EXPECT_NO_THROW(m2.set_size(2,1,true));
  ASSERT_EQ( 2u, m2.rows());
  ASSERT_EQ( 1u, m2.cols());
  EXPECT_EQ( 1, m2(0, 0) );
  EXPECT_EQ( 3, m2(1, 0) );
  EXPECT_NO_THROW(m2.set_size(2,2));
  ASSERT_EQ( 2u, m2.rows() );
  ASSERT_EQ( 2u, m2.cols() );

  // set_identity()
  EXPECT_NO_THROW(m2.set_identity());
  ASSERT_EQ( 2u, m2.rows() );
  ASSERT_EQ( 2u, m2.cols() );
  EXPECT_EQ( 1, m2(0, 0) );
  EXPECT_EQ( 0, m2(0, 1) );
  EXPECT_EQ( 0, m2(1, 0) );
  EXPECT_EQ( 1, m2(1, 1) );
  EXPECT_NO_THROW(m2.set_identity(3));
  ASSERT_EQ( 3u, m2.rows() );
  ASSERT_EQ( 3u, m2.cols() );
  EXPECT_EQ( 1, m2(0, 0) );
  EXPECT_EQ( 0, m2(0, 1) );
  EXPECT_EQ( 0, m2(0, 2) );
  EXPECT_EQ( 0, m2(1, 0) );
  EXPECT_EQ( 1, m2(1, 1) );
  EXPECT_EQ( 0, m2(1, 2) );
  EXPECT_EQ( 0, m2(2, 0) );
  EXPECT_EQ( 0, m2(2, 1) );
  EXPECT_EQ( 1, m2(2, 2) );

  // Iterators
  EXPECT_EQ(&(*(m1.begin())),&(m1(0,0)));
  EXPECT_EQ(&(*(m1.begin()+1)),&(m1(0,1)));
  EXPECT_EQ(m1.end(),m1.begin()+6);
}

TEST(Matrix, VectorEQ) {
  Matrix<float,1,2> m1(1,2), m2(1.1f,1.9f), m3(1,2);
  EXPECT_FALSE( m1==m2 );
  EXPECT_TRUE ( m1==m3 );
  EXPECT_FALSE( equal(m1,m2) );
  EXPECT_TRUE ( equal(m1,m3) );
  EXPECT_FALSE( equal(m1,m2,.05) );
  EXPECT_TRUE ( equal(m1,m3,.05) );
  EXPECT_TRUE ( equal(m1,m2,.5) );
  EXPECT_TRUE ( equal(m1,m3,.5) );
  EXPECT_TRUE ( m1!=m2 );
  EXPECT_FALSE( m1!=m3 );
  EXPECT_TRUE ( not_equal(m1,m2) );
  EXPECT_FALSE( not_equal(m1,m3) );
  EXPECT_TRUE ( not_equal(m1,m2,.05) );
  EXPECT_FALSE( not_equal(m1,m3,.05) );
  EXPECT_FALSE( not_equal(m1,m2,.5) );
  EXPECT_FALSE( not_equal(m1,m3,.5) );
}

TEST(Matrix, BasicMath) {
  Matrix2x2f m1(1,2,3,4), m2(2,1,3,-1);

  EXPECT_EQ(Matrix2x2f(-1,-2,-3,-4), -m1);
  EXPECT_EQ(Matrix2x2f(2,1,3,1), abs(m2) );

  EXPECT_EQ( 2, elem_sum(m1,1)(0,0) );
  EXPECT_EQ( 3, elem_sum(m1,1)(0,1) );
  EXPECT_EQ( 4, elem_sum(m1,1)(1,0) );
  EXPECT_EQ( 5, elem_sum(m1,1)(1,1) );

  EXPECT_EQ( Matrix2x2f(3,3,6,3)      , elem_sum(m1,m2) );
  EXPECT_EQ( Matrix2x2f(3,3,6,3)      , m1+m2 );
  EXPECT_EQ( Matrix2x2f(2,3,4,5)      , elem_sum(m1,1) );
  EXPECT_EQ( Matrix2x2f(3,4,5,6)      , elem_sum(2,m1) );

  EXPECT_EQ( Matrix2x2f(-1,1,0,5)     , elem_diff(m1,m2) );
  EXPECT_EQ( Matrix2x2f(-1,1,0,5)     , m1-m2 );
  EXPECT_EQ( Matrix2x2f(0,1,2,3)      , elem_diff(m1,1) );
  EXPECT_EQ( Matrix2x2f(1,0,-1,-2)    , elem_diff(2,m1) );

  EXPECT_EQ( Matrix2x2f(2,2,9,-4)     , elem_prod(m1,m2) );
  EXPECT_EQ( Matrix2x2f(2,4,6,8)      , elem_prod(m1,2) );
  EXPECT_EQ( Matrix2x2f(2,4,6,8)      , m1*2 );
  EXPECT_EQ( Matrix2x2f(3,6,9,12)     , elem_prod(3,m1) );
  EXPECT_EQ( Matrix2x2f(3,6,9,12)     , 3*m1 );

  EXPECT_EQ( Matrix2x2f(0.5,2,1,-4)   , elem_quot(m1,m2) );
  EXPECT_EQ( Matrix2x2f(0.5,1,1.5,2)  , elem_quot(m1,2) );
  EXPECT_EQ( Matrix2x2f(0.5,1,1.5,2)  , m1/2 );
  EXPECT_EQ( Matrix2x2f(3,1.5,1,0.75) , elem_quot(3,m1) );
}

TEST(Matrix, BasicComparison) {
  Matrix2x2f v1(1,2,3,4), v2(2,2,2,5);
  typedef Matrix<bool,2,2> Matrix2x2b;

  EXPECT_EQ( Matrix2x2b(false,true,false,false) , elem_eq(v1,v2) );
  EXPECT_EQ( Matrix2x2b(false,true,false,false) , elem_eq(v1,2) );
  EXPECT_EQ( Matrix2x2b(false,true,false,false) , elem_eq(2,v1) );

  EXPECT_EQ( Matrix2x2b(true,false,true,true)   , elem_neq(v1,v2) );
  EXPECT_EQ( Matrix2x2b(true,false,true,true)   , elem_neq(v1,2) );
  EXPECT_EQ( Matrix2x2b(true,false,true,true)   , elem_neq(2,v1) );

  EXPECT_EQ( Matrix2x2b(true,false,false,true)  , elem_lt(v1,v2) );
  EXPECT_EQ( Matrix2x2b(true,false,false,false) , elem_lt(v1,2) );
  EXPECT_EQ( Matrix2x2b(false,false,true,true)  , elem_lt(2,v1) );

  EXPECT_EQ( Matrix2x2b(false,false,true,false) , elem_gt(v1,v2) );
  EXPECT_EQ( Matrix2x2b(false,false,true,true)  , elem_gt(v1,2) );
  EXPECT_EQ( Matrix2x2b(true,false,false,false) , elem_gt(2,v1) );

  EXPECT_EQ( Matrix2x2b(true,true,false,true)   , elem_lte(v1,v2) );
  EXPECT_EQ( Matrix2x2b(true,true,false,false)  , elem_lte(v1,2) );
  EXPECT_EQ( Matrix2x2b(false,true,true,true)   , elem_lte(2,v1) );

  EXPECT_EQ( Matrix2x2b(false,true,true,false)  , elem_gte(v1,v2) );
  EXPECT_EQ( Matrix2x2b(false,true,true,true)   , elem_gte(v1,2) );
  EXPECT_EQ( Matrix2x2b(true,true,false,false)  , elem_gte(2,v1) );
}

TEST(Matrix, Norms) {
  Matrix2x2f m(1,2,3,4);

  //EXPECT_DOUBLE_EQ(5.46499, norm_2(m)); // This norm is not yet supported
  EXPECT_EQ( 6, norm_1(m) );
  EXPECT_EQ( 7, norm_inf(m) );
  EXPECT_EQ( 30, norm_frobenius_sqr(m) );
  EXPECT_FLOAT_EQ( 5.4772255750f, norm_frobenius(m) );

  EXPECT_EQ( 10, sum(m) );
  EXPECT_EQ( 24, prod(m) );
  EXPECT_EQ( 5, trace(m) );
}

TEST(Matrix, SubMatrix) {
  Matrix2x2f m(1,2,3,4);

  Matrix<float> sm = submatrix(m,0,1,2,1);
  ASSERT_EQ( 2u, sm.rows() );
  ASSERT_EQ( 1u, sm.cols() );
  EXPECT_EQ( 2, sm(0,0) );
  EXPECT_EQ( 4, sm(1,0) );

  submatrix(m,1,0,1,2) = Matrix<float,1,2>(4,5);
  EXPECT_EQ( 1, m(0,0) );
  EXPECT_EQ( 2, m(0,1) );
  EXPECT_EQ( 4, m(1,0) );
  EXPECT_EQ( 5, m(1,1) );

  Vector<float> cv = select_col(m,1);
  ASSERT_EQ( 2u, cv.size() );
  EXPECT_EQ( 2, cv(0) );
  EXPECT_EQ( 5, cv(1) );

  select_col(m,0) = Vector<float,2>(2,3);
  EXPECT_EQ( 2, m(0,0) );
  EXPECT_EQ( 2, m(0,1) );
  EXPECT_EQ( 3, m(1,0) );
  EXPECT_EQ( 5, m(1,1) );

  Vector<float> rv = select_row(m,1);
  ASSERT_EQ( 2u, rv.size() );
  EXPECT_EQ( 3, rv(0) );
  EXPECT_EQ( 5, rv(1) );

  select_row(m,0) = Vector<float,2>(2,3);
  EXPECT_EQ( 2, m(0,0) );
  EXPECT_EQ( 3, m(0,1) );
  EXPECT_EQ( 3, m(1,0) );
  EXPECT_EQ( 5, m(1,1) );

  m = Matrix2x2(1,2,3,4);
  Matrix2x2f dest;
  select_row(dest,0) = select_row(m,0);
  EXPECT_EQ( 1, dest(0,0) );
  EXPECT_EQ( 2, dest(0,1) );
  EXPECT_EQ( 0, dest(1,0) );
  EXPECT_EQ( 0, dest(1,1) );

  dest = Matrix2x2f();
  select_col(dest,0) = select_col(m,0);
  EXPECT_EQ(1, dest(0,0) );
  EXPECT_EQ(0, dest(0,1) );
  EXPECT_EQ(3, dest(1,0) );
  EXPECT_EQ(0, dest(1,1) );

  dest = Matrix2x2();
  select_col(dest,0) = select_row(m,1);
  EXPECT_EQ(3, dest(0,0) );
  EXPECT_EQ(0, dest(0,1) );
  EXPECT_EQ(4, dest(1,0) );
  EXPECT_EQ(0, dest(1,1) );

  // Subvector of select vector
  /*
  dest = Matrix3x3();
  subvector(select_row(dest,0),0,2)[0] = 3;
  subvector(select_row(dest,0),0,2) = Vector2(4,5);
  subvector(select_col(dest,2),1,2) = Vector2(5,6);
  EXPECT_EQ(4, dest(0,0));
  EXPECT_EQ(5, dest(0,1));
  EXPECT_EQ(0, dest(0,2));
  EXPECT_EQ(0, dest(1,0));
  EXPECT_EQ(0, dest(1,1));
  EXPECT_EQ(5, dest(1,2));
  EXPECT_EQ(0, dest(2,0));
  EXPECT_EQ(0, dest(2,1));
  EXPECT_EQ(6, dest(2,2));
  */
}

TEST(Matrix, Products) {
  Matrix2x2f m(1,2,3,4);
  Vector2f v(1,2);

  // Matrix*Vector
  Vector<float> r1 = m*v;
  ASSERT_EQ( 2u, r1.size() );
  EXPECT_EQ( 5, r1(0) );
  EXPECT_EQ( 11, r1(1) );

  // Matrix*Vector
  Vector<float> r2 = transpose(transpose(v)*m);
  ASSERT_EQ( 2u, r2.size() );
  EXPECT_EQ( 7, r2(0) );
  EXPECT_EQ( 10, r2(1) );

  // Matrix*Matrix
  Matrix<float> r3 = m*m;
  ASSERT_EQ( 2u, r3.rows() );
  ASSERT_EQ( 2u, r3.cols() );
  EXPECT_EQ( 7, r3(0,0) );
  EXPECT_EQ( 10, r3(0,1) );
  EXPECT_EQ( 15, r3(1,0) );
  EXPECT_EQ( 22, r3(1,1) );

  // Vector*VectorTranspose (i.e. outer product)
  Matrix<float> r4 = Vector2f(1,2)*transpose(Vector2f(2,3));
  ASSERT_EQ( 2u, r4.rows() );
  ASSERT_EQ( 2u, r4.cols() );
  EXPECT_EQ( 2, r4(0,0) );
  EXPECT_EQ( 3, r4(0,1) );
  EXPECT_EQ( 4, r4(1,0) );
  EXPECT_EQ( 6, r4(1,1) );

  // Matrix*Matrix self-assignment
  Matrix<float> r5 = m;
  r5 = r5*r5;
  ASSERT_EQ( 2u, r5.rows() );
  ASSERT_EQ( 2u, r5.cols() );
  EXPECT_EQ( 7, r5(0,0) );
  EXPECT_EQ( 10, r5(0,1) );
  EXPECT_EQ( 15, r5(1,0) );
  EXPECT_EQ( 22, r5(1,1) );

  // Matrix*Matrix self-assignment (no temporary)
  Matrix<float> r6 = m;
  r6 = no_tmp( r6*r6 );
  ASSERT_EQ( 2u, r6.rows() );
  ASSERT_EQ( 2u, r6.cols() );
  EXPECT_EQ( 7, r6(0,0) );
  EXPECT_EQ( 22, r6(0,1) );
  EXPECT_EQ( 33, r6(1,0) );
  EXPECT_EQ( 742, r6(1,1) );

  // Matrix*Matrix self-assignment
  Matrix2x2f r7 = m;
  r7 = r7*r7;
  ASSERT_EQ( 2u, r7.rows() );
  ASSERT_EQ( 2u, r7.cols() );
  EXPECT_EQ( 7, r7(0,0) );
  EXPECT_EQ( 10, r7(0,1) );
  EXPECT_EQ( 15, r7(1,0) );
  EXPECT_EQ( 22, r7(1,1) );

  // Matrix*Matrix self-assignment (no temporary)
  Matrix2x2f r8 = m;
  r8 = no_tmp( r8*r8 );
  ASSERT_EQ( 2u, r8.rows() );
  ASSERT_EQ( 2u, r8.cols() );
  EXPECT_EQ( 7, r8(0,0) );
  EXPECT_EQ( 22, r8(0,1) );
  EXPECT_EQ( 33, r8(1,0) );
  EXPECT_EQ( 742, r8(1,1) );
}

TEST(Matrix, Transpose) {
  Matrix2x2f m(1,2,3,4);
  Matrix<float> r = transpose(m);

  ASSERT_EQ( 2u, r.rows() );
  ASSERT_EQ( 2u, r.cols() );
  EXPECT_EQ( 1, r(0,0) );
  EXPECT_EQ( 3, r(0,1) );
  EXPECT_EQ( 2, r(1,0) );
  EXPECT_EQ( 4, r(1,1) );

  // Invoking the const iterator (in MatrixCol)
  math::MatrixTranspose<Matrix2x2f > trans(m);
  Vector<float> cv = select_col(trans,0);
  ASSERT_EQ( 2u, cv.size() );
  EXPECT_EQ( 1, cv(0) );
  EXPECT_EQ( 2, cv(1) );

  // Invoking the const iterator (in MatrixRow)
  Vector<float> rv = select_row(trans,0);
  ASSERT_EQ( 2u, rv.size() );
  EXPECT_EQ( 1, rv(0) );
  EXPECT_EQ( 3, rv(1) );

  // Mashing more unique operations
  Matrix2x2 dest;
  select_row(dest,0) = transpose(Vector2(4,5));
  EXPECT_EQ(4, dest(0,0));
  EXPECT_EQ(5, dest(0,1));
  EXPECT_EQ(0, dest(1,0));
  EXPECT_EQ(0, dest(1,1));
}

TEST(Matrix, Inverse) {
  Matrix2x2f m1(1,2,3,4);
  Matrix<float> i1=inverse(m1);
  EXPECT_FLOAT_EQ( -2   , i1(0,0) );
  EXPECT_FLOAT_EQ( 1    , i1(0,1) );
  EXPECT_FLOAT_EQ( 1.5  , i1(1,0) );
  EXPECT_FLOAT_EQ( -0.5 , i1(1,1) );

  Matrix3x3f m2;
  m2(0,0)=5; m2(0,1)=4; m2(0,2)=3;
  m2(1,0)=2; m2(1,1)=1; m2(1,2)=0;
  m2(2,0)=9; m2(2,1)=6; m2(2,2)=8;

  Matrix<float> i2=inverse(m2);
  EXPECT_FLOAT_EQ(   8.0f / -15.0f , i2(0,0) );
  EXPECT_FLOAT_EQ( -14.0f / -15.0f , i2(0,1) );
  EXPECT_FLOAT_EQ( -3.0f  / -15.0f , i2(0,2) );
  EXPECT_FLOAT_EQ( -16.0f / -15.0f , i2(1,0) );
  EXPECT_FLOAT_EQ(  13.0f / -15.0f , i2(1,1) );
  EXPECT_FLOAT_EQ(   6.0f / -15.0f , i2(1,2) );
  EXPECT_FLOAT_EQ(   3.0f / -15.0f , i2(2,0) );
  EXPECT_FLOAT_EQ(   6.0f / -15.0f , i2(2,1) );
  EXPECT_FLOAT_EQ(  -3.0f / -15.0f , i2(2,2) );
}
