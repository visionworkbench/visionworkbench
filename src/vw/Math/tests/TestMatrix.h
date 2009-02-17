// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestMatrix.h
#include <cxxtest/TestSuite.h>
#include <vw/Math/Matrix.h>

using namespace vw;

class TestMatrix : public CxxTest::TestSuite
{
public:

  void test_static_matrix()
  {
    // Default constructor
    Matrix<float,2,3> m1;
    TS_ASSERT_EQUALS(m1.rows(),2);
    TS_ASSERT_EQUALS(m1.cols(),3);
    TS_ASSERT_EQUALS(m1(0,0),0);
    TS_ASSERT_EQUALS(m1(0,1),0);
    TS_ASSERT_EQUALS(m1(0,2),0);
    TS_ASSERT_EQUALS(m1(1,0),0);
    TS_ASSERT_EQUALS(m1(1,1),0);
    TS_ASSERT_EQUALS(m1(1,2),0);

    // Data pointer constructor
    float data[4] = {1,2,3,4};
    Matrix<float,2,2> m2(data);
    TS_ASSERT_EQUALS(m2.rows(),2);
    TS_ASSERT_EQUALS(m2.cols(),2);
    TS_ASSERT_EQUALS(m2(0,0),1);
    TS_ASSERT_EQUALS(m2(0,1),2);
    TS_ASSERT_EQUALS(m2(1,0),3);
    TS_ASSERT_EQUALS(m2(1,1),4);
    TS_ASSERT_EQUALS(m2[0][0],1);
    TS_ASSERT_EQUALS(m2[0][1],2);
    TS_ASSERT_EQUALS(m2[1][0],3);
    TS_ASSERT_EQUALS(m2[1][1],4);

    // Copy constructor
    Matrix<float,2,2> m3(m2);
    TS_ASSERT_EQUALS(m3.rows(),2);
    TS_ASSERT_EQUALS(m3.cols(),2);
    TS_ASSERT_EQUALS(m3(0,0),1);
    TS_ASSERT_EQUALS(m3(0,1),2);
    TS_ASSERT_EQUALS(m3(1,0),3);
    TS_ASSERT_EQUALS(m3(1,1),4);

    // Element value constructor
    Matrix<float,2,2> m4(1,2,3,4);
    TS_ASSERT_EQUALS(m4.rows(),2);
    TS_ASSERT_EQUALS(m4.cols(),2);
    TS_ASSERT_EQUALS(m4(0,0),1);
    TS_ASSERT_EQUALS(m4(0,1),2);
    TS_ASSERT_EQUALS(m4(1,0),3);
    TS_ASSERT_EQUALS(m4(1,1),4);

    // set_size()
    TS_ASSERT_THROWS(m1.set_size(3,3), ArgumentErr);
    TS_ASSERT_THROWS_NOTHING(m1.set_size(2,3));

    // set_identity()
    TS_ASSERT_THROWS_NOTHING(m2.set_identity());
    TS_ASSERT_EQUALS(m2(0,0),1);
    TS_ASSERT_EQUALS(m2(0,1),0);
    TS_ASSERT_EQUALS(m2(1,0),0);
    TS_ASSERT_EQUALS(m2(1,1),1);
    TS_ASSERT_THROWS(m2.set_identity(3), ArgumentErr);
    TS_ASSERT_THROWS(m1.set_identity(), LogicErr);

    // Iterators
    TS_ASSERT_EQUALS(&(*(m1.begin())),&(m1(0,0)));
    TS_ASSERT_EQUALS(&(*(m1.begin()+1)),&(m1(0,1)));
    TS_ASSERT_EQUALS(m1.end(),m1.begin()+6);
  }

  void test_dynamic_matrix()
  {
    // Default constructor
    Matrix<float> m0;
    TS_ASSERT_EQUALS(m0.rows(),0);
    TS_ASSERT_EQUALS(m0.cols(),0);

    // Size constructor
    Matrix<float> m1(2,3);
    TS_ASSERT_EQUALS(m1.rows(),2);
    TS_ASSERT_EQUALS(m1.cols(),3);
    TS_ASSERT_EQUALS(m1(0,0),0);
    TS_ASSERT_EQUALS(m1(0,1),0);
    TS_ASSERT_EQUALS(m1(0,2),0);
    TS_ASSERT_EQUALS(m1(1,0),0);
    TS_ASSERT_EQUALS(m1(1,1),0);
    TS_ASSERT_EQUALS(m1(1,2),0);

    // Data pointer constructor
    float data[4] = {1,2,3,4};
    Matrix<float> m2(2,2,data);
    TS_ASSERT_EQUALS(m2.rows(),2);
    TS_ASSERT_EQUALS(m2.cols(),2);
    TS_ASSERT_EQUALS(m2(0,0),1);
    TS_ASSERT_EQUALS(m2(0,1),2);
    TS_ASSERT_EQUALS(m2(1,0),3);
    TS_ASSERT_EQUALS(m2(1,1),4);
    TS_ASSERT_EQUALS(m2[0][0],1);
    TS_ASSERT_EQUALS(m2[0][1],2);
    TS_ASSERT_EQUALS(m2[1][0],3);
    TS_ASSERT_EQUALS(m2[1][1],4);

    // Copy constructor
    Matrix<float> m3(m2);
    TS_ASSERT_EQUALS(m3.rows(),2);
    TS_ASSERT_EQUALS(m3.cols(),2);
    TS_ASSERT_EQUALS(m3(0,0),1);
    TS_ASSERT_EQUALS(m3(0,1),2);
    TS_ASSERT_EQUALS(m3(1,0),3);
    TS_ASSERT_EQUALS(m3(1,1),4);

    // set_size()
    TS_ASSERT_THROWS_NOTHING(m2.set_size(2,1,true));
    TS_ASSERT_EQUALS(m2.rows(),2);
    TS_ASSERT_EQUALS(m2.cols(),1);
    TS_ASSERT_EQUALS(m2(0,0),1);
    TS_ASSERT_EQUALS(m2(1,0),3);
    TS_ASSERT_THROWS_NOTHING(m2.set_size(2,2));
    TS_ASSERT_EQUALS(m2.rows(),2);
    TS_ASSERT_EQUALS(m2.cols(),2);

     // set_identity()
    TS_ASSERT_THROWS_NOTHING(m2.set_identity());
    TS_ASSERT_EQUALS(m2.rows(),2);
    TS_ASSERT_EQUALS(m2.cols(),2);
    TS_ASSERT_EQUALS(m2(0,0),1);
    TS_ASSERT_EQUALS(m2(0,1),0);
    TS_ASSERT_EQUALS(m2(1,0),0);
    TS_ASSERT_EQUALS(m2(1,1),1);
    TS_ASSERT_THROWS_NOTHING(m2.set_identity(3));
    TS_ASSERT_EQUALS(m2.rows(),3);
    TS_ASSERT_EQUALS(m2.cols(),3);
    TS_ASSERT_EQUALS(m2(0,0),1);
    TS_ASSERT_EQUALS(m2(0,1),0);
    TS_ASSERT_EQUALS(m2(0,2),0);
    TS_ASSERT_EQUALS(m2(1,0),0);
    TS_ASSERT_EQUALS(m2(1,1),1);
    TS_ASSERT_EQUALS(m2(1,2),0);
    TS_ASSERT_EQUALS(m2(2,0),0);
    TS_ASSERT_EQUALS(m2(2,1),0);
    TS_ASSERT_EQUALS(m2(2,2),1);

   // Iterators
    TS_ASSERT_EQUALS(&(*(m1.begin())),&(m1(0,0)));
    TS_ASSERT_EQUALS(&(*(m1.begin()+1)),&(m1(0,1)));
    TS_ASSERT_EQUALS(m1.end(),m1.begin()+6);
  }

  void test_vector_equality()
  {
    Matrix<float,1,2> m1(1,2), m2(1.1,1.9), m3(1,2);
    TS_ASSERT_EQUALS( m1==m2, false );
    TS_ASSERT_EQUALS( m1==m3, true );
    TS_ASSERT_EQUALS( equal(m1,m2), false );
    TS_ASSERT_EQUALS( equal(m1,m3), true );
    TS_ASSERT_EQUALS( equal(m1,m2,.05), false );
    TS_ASSERT_EQUALS( equal(m1,m3,.05), true );
    TS_ASSERT_EQUALS( equal(m1,m2,.5), true );
    TS_ASSERT_EQUALS( equal(m1,m3,.5), true );
    TS_ASSERT_EQUALS( m1!=m2, true );
    TS_ASSERT_EQUALS( m1!=m3, false );
    TS_ASSERT_EQUALS( not_equal(m1,m2), true );
    TS_ASSERT_EQUALS( not_equal(m1,m3), false );
    TS_ASSERT_EQUALS( not_equal(m1,m2,.05), true );
    TS_ASSERT_EQUALS( not_equal(m1,m3,.05), false );
    TS_ASSERT_EQUALS( not_equal(m1,m2,.5), false );
    TS_ASSERT_EQUALS( not_equal(m1,m3,.5), false );
  }

  void test_basic_matrix_math()
  {
    Matrix<float,2,2> m1(1,2,3,4), m2(2,1,3,-1);

    TS_ASSERT_EQUALS(-m1,(Matrix<float,2,2>(-1,-2,-3,-4)));

    TS_ASSERT_EQUALS(elem_sum(m1,m2),(Matrix<float,2,2>(3,3,6,3)));
    TS_ASSERT_EQUALS(m1+m2,(Matrix<float,2,2>(3,3,6,3)));
    TS_ASSERT_EQUALS(elem_sum(m1,1)(0,0),2);
    TS_ASSERT_EQUALS(elem_sum(m1,1)(0,1),3);
    TS_ASSERT_EQUALS(elem_sum(m1,1)(1,0),4);
    TS_ASSERT_EQUALS(elem_sum(m1,1)(1,1),5);
    TS_ASSERT_EQUALS(elem_sum(m1,1),(Matrix<float,2,2>(2,3,4,5)));
    TS_ASSERT_EQUALS(elem_sum(2,m1),(Matrix<float,2,2>(3,4,5,6)));

    TS_ASSERT_EQUALS(elem_diff(m1,m2),(Matrix<float,2,2>(-1,1,0,5)));
    TS_ASSERT_EQUALS(m1-m2,(Matrix<float,2,2>(-1,1,0,5)));
    TS_ASSERT_EQUALS(elem_diff(m1,1),(Matrix<float,2,2>(0,1,2,3)));
    TS_ASSERT_EQUALS(elem_diff(2,m1),(Matrix<float,2,2>(1,0,-1,-2)));

    TS_ASSERT_EQUALS(elem_prod(m1,m2),(Matrix<float,2,2>(2,2,9,-4)));
    TS_ASSERT_EQUALS(elem_prod(m1,2),(Matrix<float,2,2>(2,4,6,8)));
    TS_ASSERT_EQUALS(m1*2,(Matrix<float,2,2>(2,4,6,8)));
    TS_ASSERT_EQUALS(elem_prod(3,m1),(Matrix<float,2,2>(3,6,9,12)));
    TS_ASSERT_EQUALS(3*m1,(Matrix<float,2,2>(3,6,9,12)));

    TS_ASSERT_EQUALS(elem_quot(m1,m2),(Matrix<float,2,2>(0.5,2,1,-4)));
    TS_ASSERT_EQUALS(elem_quot(m1,2),(Matrix<float,2,2>(0.5,1,1.5,2)));
    TS_ASSERT_EQUALS(m1/2,(Matrix<float,2,2>(0.5,1,1.5,2)));
    TS_ASSERT_EQUALS(elem_quot(3,m1),(Matrix<float,2,2>(3,1.5,1,0.75)));
  }

  void test_basic_matrix_comparison()
  {
    Matrix<float,2,2> v1(1,2,3,4), v2(2,2,2,5);

    TS_ASSERT_EQUALS(elem_eq(v1,v2),(Matrix<bool,2,2>(false,true,false,false)));
    TS_ASSERT_EQUALS(elem_eq(v1,2),(Matrix<bool,2,2>(false,true,false,false)));
    TS_ASSERT_EQUALS(elem_eq(2,v1),(Matrix<bool,2,2>(false,true,false,false)));

    TS_ASSERT_EQUALS(elem_neq(v1,v2),(Matrix<bool,2,2>(true,false,true,true)));
    TS_ASSERT_EQUALS(elem_neq(v1,2),(Matrix<bool,2,2>(true,false,true,true)));
    TS_ASSERT_EQUALS(elem_neq(2,v1),(Matrix<bool,2,2>(true,false,true,true)));

    TS_ASSERT_EQUALS(elem_lt(v1,v2),(Matrix<bool,2,2>(true,false,false,true)));
    TS_ASSERT_EQUALS(elem_lt(v1,2),(Matrix<bool,2,2>(true,false,false,false)));
    TS_ASSERT_EQUALS(elem_lt(2,v1),(Matrix<bool,2,2>(false,false,true,true)));

    TS_ASSERT_EQUALS(elem_gt(v1,v2),(Matrix<bool,2,2>(false,false,true,false)));
    TS_ASSERT_EQUALS(elem_gt(v1,2),(Matrix<bool,2,2>(false,false,true,true)));
    TS_ASSERT_EQUALS(elem_gt(2,v1),(Matrix<bool,2,2>(true,false,false,false)));

    TS_ASSERT_EQUALS(elem_lte(v1,v2),(Matrix<bool,2,2>(true,true,false,true)));
    TS_ASSERT_EQUALS(elem_lte(v1,2),(Matrix<bool,2,2>(true,true,false,false)));
    TS_ASSERT_EQUALS(elem_lte(2,v1),(Matrix<bool,2,2>(false,true,true,true)));

    TS_ASSERT_EQUALS(elem_gte(v1,v2),(Matrix<bool,2,2>(false,true,true,false)));
    TS_ASSERT_EQUALS(elem_gte(v1,2),(Matrix<bool,2,2>(false,true,true,true)));
    TS_ASSERT_EQUALS(elem_gte(2,v1),(Matrix<bool,2,2>(true,true,false,false)));
  }

  void test_matrix_norms()
  {
    Matrix<float,2,2> m(1,2,3,4);

    TS_ASSERT_EQUALS(norm_1(m), 6);
    //TS_ASSERT_DELTA(norm_2(m), 5.46499, .0001); This norm is not yet supported
    TS_ASSERT_EQUALS(norm_inf(m), 7);
    TS_ASSERT_DELTA(norm_frobenius(m), 5.47723, .0001);
    TS_ASSERT_EQUALS(norm_frobenius_sqr(m), 30);

    TS_ASSERT_EQUALS(sum(m),10);
    TS_ASSERT_EQUALS(prod(m),24);
    TS_ASSERT_EQUALS(trace(m),5);
  }

  void test_submatrix()
  {
    Matrix<float,2,2> m(1,2,3,4);

    Matrix<float> sm = submatrix(m,0,1,2,1);
    TS_ASSERT_EQUALS( sm.rows(), 2 );
    TS_ASSERT_EQUALS( sm.cols(), 1 );
    TS_ASSERT_EQUALS( sm(0,0), 2 );
    TS_ASSERT_EQUALS( sm(1,0), 4 );

    submatrix(m,1,0,1,2) = Matrix<float,1,2>(4,5);
    TS_ASSERT_EQUALS( m(0,0), 1 );
    TS_ASSERT_EQUALS( m(0,1), 2 );
    TS_ASSERT_EQUALS( m(1,0), 4 );
    TS_ASSERT_EQUALS( m(1,1), 5 );

    Vector<float> cv = select_col(m,1);
    TS_ASSERT_EQUALS( int(cv.size()), 2 );
    TS_ASSERT_EQUALS( cv(0), 2 );
    TS_ASSERT_EQUALS( cv(1), 5 );

    select_col(m,0) = Vector<float,2>(2,3);
    TS_ASSERT_EQUALS( m(0,0), 2 );
    TS_ASSERT_EQUALS( m(0,1), 2 );
    TS_ASSERT_EQUALS( m(1,0), 3 );
    TS_ASSERT_EQUALS( m(1,1), 5 );

    Vector<float> rv = select_row(m,1);
    TS_ASSERT_EQUALS( int(rv.size()), 2 );
    TS_ASSERT_EQUALS( rv(0), 3 );
    TS_ASSERT_EQUALS( rv(1), 5 );

    select_row(m,0) = Vector<float,2>(2,3);
    TS_ASSERT_EQUALS( m(0,0), 2 );
    TS_ASSERT_EQUALS( m(0,1), 3 );
    TS_ASSERT_EQUALS( m(1,0), 3 );
    TS_ASSERT_EQUALS( m(1,1), 5 );
  }

  void test_matrix_products()
  {
    Matrix<float,2,2> m(1,2,3,4);
    Vector<float,2> v(1,2);

    // Matrix*Vector
    Vector<float> r1 = m*v;
    TS_ASSERT_EQUALS( int(r1.size()), 2 );
    TS_ASSERT_EQUALS( r1(0), 5 );
    TS_ASSERT_EQUALS( r1(1), 11 );

    // Matrix*Vector
    Vector<float> r2 = transpose(transpose(v)*m);
    TS_ASSERT_EQUALS( int(r2.size()), 2 );
    TS_ASSERT_EQUALS( r2(0), 7 );
    TS_ASSERT_EQUALS( r2(1), 10 );

    // Matrix*Matrix
    Matrix<float> r3 = m*m;
    TS_ASSERT_EQUALS( r3.rows(), 2 );
    TS_ASSERT_EQUALS( r3.cols(), 2 );
    TS_ASSERT_EQUALS( r3(0,0), 7 );
    TS_ASSERT_EQUALS( r3(0,1), 10 );
    TS_ASSERT_EQUALS( r3(1,0), 15 );
    TS_ASSERT_EQUALS( r3(1,1), 22 );

    // Vector*VectorTranspose (i.e. outer product)
    Matrix<float> r4 = Vector<float,2>(1,2)*transpose(Vector<float,2>(2,3));
    TS_ASSERT_EQUALS( r4.rows(), 2 );
    TS_ASSERT_EQUALS( r4.cols(), 2 );
    TS_ASSERT_EQUALS( r4(0,0), 2 );
    TS_ASSERT_EQUALS( r4(0,1), 3 );
    TS_ASSERT_EQUALS( r4(1,0), 4 );
    TS_ASSERT_EQUALS( r4(1,1), 6 );

    // Matrix*Matrix self-assignment
    Matrix<float> r5 = m;
    r5 = r5*r5;
    TS_ASSERT_EQUALS( r5.rows(), 2 );
    TS_ASSERT_EQUALS( r5.cols(), 2 );
    TS_ASSERT_EQUALS( r5(0,0), 7 );
    TS_ASSERT_EQUALS( r5(0,1), 10 );
    TS_ASSERT_EQUALS( r5(1,0), 15 );
    TS_ASSERT_EQUALS( r5(1,1), 22 );
    /*
    // Matrix*Matrix self-assignment (no temporary)
    Matrix<float> r6 = m;
    r6 = no_tmp( r6*r6 );
    TS_ASSERT_EQUALS( r6.rows(), 2 );
    TS_ASSERT_EQUALS( r6.cols(), 2 );
    TS_ASSERT_EQUALS( r6(0,0), 7 );
    TS_ASSERT_EQUALS( r6(0,1), 22 );
    TS_ASSERT_EQUALS( r6(1,0), 33 );
    TS_ASSERT_EQUALS( r6(1,1), 742 );
    */
    // Matrix*Matrix self-assignment
    Matrix<float,2,2> r7 = m;
    r7 = r7*r7;
    TS_ASSERT_EQUALS( r7.rows(), 2 );
    TS_ASSERT_EQUALS( r7.cols(), 2 );
    TS_ASSERT_EQUALS( r7(0,0), 7 );
    TS_ASSERT_EQUALS( r7(0,1), 10 );
    TS_ASSERT_EQUALS( r7(1,0), 15 );
    TS_ASSERT_EQUALS( r7(1,1), 22 );
    /*
    // Matrix*Matrix self-assignment (no temporary)
    Matrix<float,2,2> r8 = m;
    r8 = no_tmp( r8*r8 );
    TS_ASSERT_EQUALS( r8.rows(), 2 );
    TS_ASSERT_EQUALS( r8.cols(), 2 );
    TS_ASSERT_EQUALS( r8(0,0), 7 );
    TS_ASSERT_EQUALS( r8(0,1), 22 );
    TS_ASSERT_EQUALS( r8(1,0), 33 );
    TS_ASSERT_EQUALS( r8(1,1), 742 );
    */
  }

  void test_matrix_transpose()
  {
    Matrix<float,2,2> m(1,2,3,4);
    Matrix<float> r = transpose(m);

    TS_ASSERT_EQUALS( r.rows(), 2 );
    TS_ASSERT_EQUALS( r.cols(), 2 );
    TS_ASSERT_EQUALS( r(0,0), 1 );
    TS_ASSERT_EQUALS( r(0,1), 3 );
    TS_ASSERT_EQUALS( r(1,0), 2 );
    TS_ASSERT_EQUALS( r(1,1), 4 );
  }

  void test_matrix_inverse()
  {
    Matrix<float,2,2> m1(1,2,3,4);
    Matrix<float> i1=inverse(m1);
    TS_ASSERT_DELTA( i1(0,0), -2,   .0001 );
    TS_ASSERT_DELTA( i1(0,1),  1,   .0001 );
    TS_ASSERT_DELTA( i1(1,0),  1.5, .0001 );
    TS_ASSERT_DELTA( i1(1,1), -0.5, .0001 );

    Matrix<float,3,3> m2;
    m2(0,0)=5; m2(0,1)=4; m2(0,2)=3;
    m2(1,0)=2; m2(1,1)=1; m2(1,2)=0;
    m2(2,0)=9; m2(2,1)=6; m2(2,2)=8;
    Matrix<float> i2=inverse(m2);
    TS_ASSERT_DELTA( i2(0,0), -0.53333, .0001 );
    TS_ASSERT_DELTA( i2(0,1),  0.93333, .0001 );
    TS_ASSERT_DELTA( i2(0,2),  0.2,     .0001 );
    TS_ASSERT_DELTA( i2(1,0),  1.06666, .0001 );
    TS_ASSERT_DELTA( i2(1,1), -0.86666, .0001 );
    TS_ASSERT_DELTA( i2(1,2), -0.4,     .0001 );
    TS_ASSERT_DELTA( i2(2,0), -0.2,     .0001 );
    TS_ASSERT_DELTA( i2(2,1), -0.4,     .0001 );
    TS_ASSERT_DELTA( i2(2,2),  0.2,     .0001 );
  }

}; // class TestMatrix
