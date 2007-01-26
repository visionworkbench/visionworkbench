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

// TestVector.h
#include <cxxtest/TestSuite.h>
#include <vw/Math/Vector.h>

using namespace vw;

class TestVector : public CxxTest::TestSuite
{
public:

  void test_static_vector()
  {
    // Default constructor
    Vector<float,4> v;
    TS_ASSERT_EQUALS(v.size(),4);
    TS_ASSERT_EQUALS(v(0),0);
    TS_ASSERT_EQUALS(v(1),0);
    TS_ASSERT_EQUALS(v(2),0);
    TS_ASSERT_EQUALS(v(3),0);

    // Values constructor
    Vector<float,4> v2(6,7,8,9);
    TS_ASSERT_EQUALS(v2.size(),4);
    TS_ASSERT_EQUALS(v2[0],6);
    TS_ASSERT_EQUALS(v2[1],7);
    TS_ASSERT_EQUALS(v2[2],8);
    TS_ASSERT_EQUALS(v2[3],9);
    TS_ASSERT_EQUALS(v2(0),6);
    TS_ASSERT_EQUALS(v2(1),7);
    TS_ASSERT_EQUALS(v2(2),8);
    TS_ASSERT_EQUALS(v2(3),9);
    TS_ASSERT_EQUALS(v2.x(),6);
    TS_ASSERT_EQUALS(v2.y(),7);
    TS_ASSERT_EQUALS(v2.z(),8);

    // Copy constructor
    Vector<float,4> v3(v2);
    TS_ASSERT_EQUALS(v3.size(),4);
    TS_ASSERT_EQUALS(v3(0),6);
    TS_ASSERT_EQUALS(v3(1),7);
    TS_ASSERT_EQUALS(v3(2),8);
    TS_ASSERT_EQUALS(v3(3),9);

    // set_size()
    TS_ASSERT_THROWS(v.set_size(3), ArgumentErr);
    TS_ASSERT_THROWS_NOTHING(v.set_size(4));

    // Iterators
    TS_ASSERT_EQUALS(&(*(v.begin())),&(v(0)));
    TS_ASSERT_EQUALS(&(*(v.begin()+1)),&(v(1)));
    TS_ASSERT_EQUALS(v.end(),v.begin()+4);
  }    

  void test_dynamic_vector()
  {
    Vector<float> v(4);
    TS_ASSERT_EQUALS(v.size(),4);
    TS_ASSERT_EQUALS(v(0),0);
    TS_ASSERT_EQUALS(v(1),0);
    TS_ASSERT_EQUALS(v(2),0);
    TS_ASSERT_EQUALS(v(3),0);

    v = Vector<float,4>(6,7,8,9);
    TS_ASSERT_EQUALS(v[0],6);
    TS_ASSERT_EQUALS(v[1],7);
    TS_ASSERT_EQUALS(v[2],8);
    TS_ASSERT_EQUALS(v[3],9);
    TS_ASSERT_EQUALS(v(0),6);
    TS_ASSERT_EQUALS(v(1),7);
    TS_ASSERT_EQUALS(v(2),8);
    TS_ASSERT_EQUALS(v(3),9);

    v.set_size(3,true);
    TS_ASSERT_EQUALS(v.size(),3);
    TS_ASSERT_EQUALS(v(0),6);
    TS_ASSERT_EQUALS(v(1),7);
    TS_ASSERT_EQUALS(v(2),8);

    TS_ASSERT_EQUALS(&(*(v.begin())),&(v(0)));
    TS_ASSERT_EQUALS(&(*(v.begin()+1)),&(v(1)));
    TS_ASSERT_EQUALS(&(*(v.end()-1)),&(v(2)));
  }    

  void test_vector_proxy()
  {
    float data[] = {1,2,3,4};

    VectorProxy<float,4> vp1(data);
    TS_ASSERT_EQUALS(vp1.size(),4);
    TS_ASSERT_EQUALS(vp1(0),1);
    TS_ASSERT_EQUALS(vp1(1),2);
    TS_ASSERT_EQUALS(vp1(2),3);
    TS_ASSERT_EQUALS(vp1(3),4);

    TS_ASSERT_THROWS_NOTHING( (vp1=Vector<float,4>(5,6,7,8)) );
    TS_ASSERT_EQUALS(vp1.size(),4);
    TS_ASSERT_EQUALS(vp1(0),5);
    TS_ASSERT_EQUALS(vp1(1),6);
    TS_ASSERT_EQUALS(vp1(2),7);
    TS_ASSERT_EQUALS(vp1(3),8);

    VectorProxy<float> vp2(4,data);
    TS_ASSERT_EQUALS(vp2.size(),4);
    TS_ASSERT_EQUALS(vp2(0),5);
    TS_ASSERT_EQUALS(vp2(1),6);
    TS_ASSERT_EQUALS(vp2(2),7);
    TS_ASSERT_EQUALS(vp2(3),8);

    TS_ASSERT_THROWS_NOTHING( (vp2=Vector<float,4>(1,2,3,4)) );
    TS_ASSERT_EQUALS(vp1.size(),4);
    TS_ASSERT_EQUALS(vp2(0),1);
    TS_ASSERT_EQUALS(vp2(1),2);
    TS_ASSERT_EQUALS(vp2(2),3);
    TS_ASSERT_EQUALS(vp2(3),4);
  }

  void test_subvector()
  {
    Vector<float,4> v(1,2,3,4);
    Vector<float> sv = subvector(v,2,2);
    TS_ASSERT_EQUALS( sv.size(), 2 );
    TS_ASSERT_EQUALS( sv(0), 3 );
    TS_ASSERT_EQUALS( sv(1), 4 );
    subvector(v,1,2) = Vector<float,2>(4,5);
    TS_ASSERT_EQUALS( v(0), 1 );
    TS_ASSERT_EQUALS( v(1), 4 );
    TS_ASSERT_EQUALS( v(2), 5 );
    TS_ASSERT_EQUALS( v(3), 4 );
  }

  void test_vector_iostream()
  {
    std::ostringstream oss;
    Vector<float,3> v1(1,2,3);
    oss << v1;
    TS_ASSERT_EQUALS( oss.str(), "Vector3(1,2,3)" );
    oss.str("");
    Vector<float> v2=v1;
    oss << v2;
    TS_ASSERT_EQUALS( oss.str(), "Vector3(1,2,3)" );
  }

  void test_vector_equality()
  {
    Vector<float,3> v1(1,2,3), v2(1.1,1.9,3), v3(1,2,3);
    TS_ASSERT_EQUALS( v1==v2, false );
    TS_ASSERT_EQUALS( v1==v3, true );
    TS_ASSERT_EQUALS( equal(v1,v2), false );
    TS_ASSERT_EQUALS( equal(v1,v3), true );
    TS_ASSERT_EQUALS( equal(v1,v2,.05), false );
    TS_ASSERT_EQUALS( equal(v1,v3,.05), true );
    TS_ASSERT_EQUALS( equal(v1,v2,.5), true );
    TS_ASSERT_EQUALS( equal(v1,v3,.5), true );
    TS_ASSERT_EQUALS( v1!=v2, true );
    TS_ASSERT_EQUALS( v1!=v3, false );
    TS_ASSERT_EQUALS( not_equal(v1,v2), true );
    TS_ASSERT_EQUALS( not_equal(v1,v3), false );
    TS_ASSERT_EQUALS( not_equal(v1,v2,.05), true );
    TS_ASSERT_EQUALS( not_equal(v1,v3,.05), false );
    TS_ASSERT_EQUALS( not_equal(v1,v2,.5), false );
    TS_ASSERT_EQUALS( not_equal(v1,v3,.5), false );
  }

  void test_basic_vector_math()
  {
    Vector<float,3> v1(1,2,3), v2(2,4,4);

    TS_ASSERT_EQUALS(-v1,(Vector<float,3>(-1,-2,-3)));

    TS_ASSERT_EQUALS(elem_sum(v1,v2),(Vector<float,3>(3,6,7)));
    TS_ASSERT_EQUALS(v1+v2,(Vector<float,3>(3,6,7)));
    TS_ASSERT_EQUALS(elem_sum(v1,1),(Vector<float,3>(2,3,4)));
    TS_ASSERT_EQUALS(elem_sum(2,v1),(Vector<float,3>(3,4,5)));

    TS_ASSERT_EQUALS(elem_diff(v1,v2),(Vector<float,3>(-1,-2,-1)));
    TS_ASSERT_EQUALS(v1-v2,(Vector<float,3>(-1,-2,-1)));
    TS_ASSERT_EQUALS(elem_diff(v1,1),(Vector<float,3>(0,1,2)));
    TS_ASSERT_EQUALS(elem_diff(2,v1),(Vector<float,3>(1,0,-1)));

    TS_ASSERT_EQUALS(elem_prod(v1,v2),(Vector<float,3>(2,8,12)));
    TS_ASSERT_EQUALS(elem_prod(v1,2),(Vector<float,3>(2,4,6)));
    TS_ASSERT_EQUALS(v1*2,(Vector<float,3>(2,4,6)));
    TS_ASSERT_EQUALS(elem_prod(3,v1),(Vector<float,3>(3,6,9)));
    TS_ASSERT_EQUALS(3*v1,(Vector<float,3>(3,6,9)));

    TS_ASSERT_EQUALS(elem_quot(v1,v2),(Vector<float,3>(0.5,0.5,0.75)));
    TS_ASSERT_EQUALS(elem_quot(v1,2),(Vector<float,3>(0.5,1,1.5)));
    TS_ASSERT_EQUALS(v1/2,(Vector<float,3>(0.5,1,1.5)));
    TS_ASSERT_EQUALS(elem_quot(3,v1),(Vector<float,3>(3,1.5,1)));
  }

  void test_vector_elem_comparison()
  {
    Vector<float,3> v1(1,2,3), v2(2,2,2);

    TS_ASSERT_EQUALS(elem_eq(v1,v2),(Vector<bool,3>(false,true,false)));
    TS_ASSERT_EQUALS(elem_eq(v1,2),(Vector<bool,3>(false,true,false)));
    TS_ASSERT_EQUALS(elem_eq(2,v1),(Vector<bool,3>(false,true,false)));

    TS_ASSERT_EQUALS(elem_neq(v1,v2),(Vector<bool,3>(true,false,true)));
    TS_ASSERT_EQUALS(elem_neq(v1,2),(Vector<bool,3>(true,false,true)));
    TS_ASSERT_EQUALS(elem_neq(2,v1),(Vector<bool,3>(true,false,true)));

    TS_ASSERT_EQUALS(elem_lt(v1,v2),(Vector<bool,3>(true,false,false)));
    TS_ASSERT_EQUALS(elem_lt(v1,2),(Vector<bool,3>(true,false,false)));
    TS_ASSERT_EQUALS(elem_lt(2,v1),(Vector<bool,3>(false,false,true)));

    TS_ASSERT_EQUALS(elem_gt(v1,v2),(Vector<bool,3>(false,false,true)));
    TS_ASSERT_EQUALS(elem_gt(v1,2),(Vector<bool,3>(false,false,true)));
    TS_ASSERT_EQUALS(elem_gt(2,v1),(Vector<bool,3>(true,false,false)));

    TS_ASSERT_EQUALS(elem_lte(v1,v2),(Vector<bool,3>(true,true,false)));
    TS_ASSERT_EQUALS(elem_lte(v1,2),(Vector<bool,3>(true,true,false)));
    TS_ASSERT_EQUALS(elem_lte(2,v1),(Vector<bool,3>(false,true,true)));

    TS_ASSERT_EQUALS(elem_gte(v1,v2),(Vector<bool,3>(false,true,true)));
    TS_ASSERT_EQUALS(elem_gte(v1,2),(Vector<bool,3>(false,true,true)));
    TS_ASSERT_EQUALS(elem_gte(2,v1),(Vector<bool,3>(true,true,false)));
  }

  void test_vector_norms()
  {
    Vector<float,3> v(1,2,-1);

    TS_ASSERT_EQUALS(norm_1(v), 4);
    TS_ASSERT_DELTA(norm_2(v), 2.44948, .0001);
    TS_ASSERT_EQUALS(norm_2_sqr(v), 6);
    TS_ASSERT_EQUALS(norm_inf(v), 2);
    TS_ASSERT_EQUALS(index_norm_inf(v), 1);

    TS_ASSERT_EQUALS(sum(v),2);
    TS_ASSERT_EQUALS(prod(v),-2);
  }

  void test_vector_funcs()
  {
    Vector<float,3> v1(1,2,3), v2(2,4,4);

    TS_ASSERT_EQUALS( dot_prod(v1,v2), 22 );
    TS_ASSERT_EQUALS( cross_prod(v1,v2), (Vector<float,3>(-4,2,0)) );

    Vector<float,3> vn = normalize(v1);
    TS_ASSERT_DELTA( vn(0), 0.26726, 0.0001 );
    TS_ASSERT_DELTA( vn(1), 0.53452, 0.0001 );
    TS_ASSERT_DELTA( vn(2), 0.80178, 0.0001 );
  }

  void test_vector_transpose()
  {
    Vector<float,3> v1(1,2,3), v2(2,4,4);

    TS_ASSERT_EQUALS( transpose(v1)*v2, 22 );
    TS_ASSERT_EQUALS( (-transpose(v1))*v2, -22 );
    TS_ASSERT_EQUALS( (transpose(v1)+transpose(v2))*v2, 58 );
    TS_ASSERT_EQUALS( (transpose(v1)-transpose(v2))*v2, -14 );
    TS_ASSERT_EQUALS( (2*transpose(v1))*v2, 44 );
    TS_ASSERT_EQUALS( (transpose(v1)*3)*v2, 66 );
    TS_ASSERT_EQUALS( (transpose(v1)/2)*v2, 11 );
  }

}; // class TestVector
