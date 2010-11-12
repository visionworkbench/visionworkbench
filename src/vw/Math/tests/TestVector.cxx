// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestVector.h
#include <gtest/gtest.h>
#include <vw/Math/Vector.h>

using namespace vw;

TEST(Vector, Static) {
  // Default constructor
  Vector4f v;
  ASSERT_EQ( 4u, v.size() );
  EXPECT_EQ( 0, v(0) );
  EXPECT_EQ( 0, v(1) );
  EXPECT_EQ( 0, v(2) );
  EXPECT_EQ( 0, v(3) );

  // Values constructor
  Vector4f v2(6,7,8,9);
  ASSERT_EQ( 4u, v2.size() );
  EXPECT_EQ( 6, v2[0] );
  EXPECT_EQ( 7, v2[1] );
  EXPECT_EQ( 8, v2[2] );
  EXPECT_EQ( 9, v2[3] );
  EXPECT_EQ( 6, v2(0) );
  EXPECT_EQ( 7, v2(1) );
  EXPECT_EQ( 8, v2(2) );
  EXPECT_EQ( 9, v2(3) );
  EXPECT_EQ( 6, v2.x() );
  EXPECT_EQ( 7, v2.y() );
  EXPECT_EQ( 8, v2.z() );

  // Copy constructor
  Vector4f v3(v2);
  ASSERT_EQ( 4u, v3.size() );
  EXPECT_EQ( 6, v3(0) );
  EXPECT_EQ( 7, v3(1) );
  EXPECT_EQ( 8, v3(2) );
  EXPECT_EQ( 9, v3(3) );

  // set_size()
  EXPECT_THROW(v.set_size(3), ArgumentErr);
  ASSERT_NO_THROW(v.set_size(4));

  // Iterators
  EXPECT_EQ(&(*(v.begin())),&(v(0)));
  EXPECT_EQ(&(*(v.begin()+1)),&(v(1)));
  EXPECT_EQ(v.end(),v.begin()+4);
}

TEST(Vector, Dynamic) {
  Vector<float> v(4);
  ASSERT_EQ( 4u, v.size() );
  EXPECT_EQ( 0, v(0) );
  EXPECT_EQ( 0, v(1) );
  EXPECT_EQ( 0, v(2) );
  EXPECT_EQ( 0, v(3) );

  v = Vector4f(6,7,8,9);
  EXPECT_EQ( 6, v[0] );
  EXPECT_EQ( 7, v[1] );
  EXPECT_EQ( 8, v[2] );
  EXPECT_EQ( 9, v[3] );
  EXPECT_EQ( 6, v(0) );
  EXPECT_EQ( 7, v(1) );
  EXPECT_EQ( 8, v(2) );
  EXPECT_EQ( 9, v(3) );

  v.set_size(3,true);
  ASSERT_EQ( 3u, v.size() );
  EXPECT_EQ( 6, v(0) );
  EXPECT_EQ( 7, v(1) );
  EXPECT_EQ( 8, v(2) );

  EXPECT_EQ(&(*(v.begin())),&(v(0)));
  EXPECT_EQ(&(*(v.begin()+1)),&(v(1)));
  EXPECT_EQ(&(*(v.end()-1)),&(v(2)));
}

TEST(Vector, Proxy) {
  float data[] = {1,2,3,4};

  VectorProxy<float,4> vp1(data);
  ASSERT_EQ( 4u, vp1.size() );
  EXPECT_EQ( 1, vp1(0) );
  EXPECT_EQ( 2, vp1(1) );
  EXPECT_EQ( 3, vp1(2) );
  EXPECT_EQ( 4, vp1(3) );

  ASSERT_NO_THROW( (vp1=Vector4f(5,6,7,8)) );
  ASSERT_EQ( 4u, vp1.size() );
  EXPECT_EQ( 5, vp1(0) );
  EXPECT_EQ( 6, vp1(1) );
  EXPECT_EQ( 7, vp1(2) );
  EXPECT_EQ( 8, vp1(3) );

  VectorProxy<float> vp2(4,data);
  ASSERT_EQ( 4u, vp2.size() );
  EXPECT_EQ( 5, vp2(0) );
  EXPECT_EQ( 6, vp2(1) );
  EXPECT_EQ( 7, vp2(2) );
  EXPECT_EQ( 8, vp2(3) );

  ASSERT_NO_THROW( (vp2=Vector4f(1,2,3,4)) );
  ASSERT_EQ( 4u, vp1.size() );
  EXPECT_EQ( 1, vp2(0) );
  EXPECT_EQ( 2, vp2(1) );
  EXPECT_EQ( 3, vp2(2) );
  EXPECT_EQ( 4, vp2(3) );
}

TEST(Vector, SubVector) {
  Vector4f v(1,2,3,4);
  Vector<float> sv = subvector(v,2,2);
  ASSERT_EQ( 2u, sv.size() );
  EXPECT_EQ( 3, sv(0) );
  EXPECT_EQ( 4, sv(1) );
  subvector(v,1,2) = Vector2f(4,5);
  EXPECT_EQ( 1, v(0) );
  EXPECT_EQ( 4, v(1) );
  EXPECT_EQ( 5, v(2) );
  EXPECT_EQ( 4, v(3) );
}

TEST(Vector, IOStream) {
  std::ostringstream oss;
  Vector3f v1(1,2,3);
  oss << v1;
  EXPECT_EQ( oss.str(), "Vector3(1,2,3)" );
  oss.str("");
  Vector<float> v2=v1;
  oss << v2;
  EXPECT_EQ( oss.str(), "Vector3(1,2,3)" );
}

TEST(Vector, Equality) {
  Vector3f v1(1,2,3), v2(1.1f,1.9f,3), v3(1,2,3);
  EXPECT_FALSE ( v1 == v2 );
  EXPECT_TRUE  ( v1 == v3 );
  EXPECT_FALSE ( equal(v1,v2) );
  EXPECT_TRUE  ( equal(v1,v3) );
  EXPECT_FALSE ( equal(v1,v2,.05) );
  EXPECT_TRUE  ( equal(v1,v3,.05) );
  EXPECT_TRUE  ( equal(v1,v2,.5)  );
  EXPECT_TRUE  ( equal(v1,v3,.5)  );
  EXPECT_TRUE  ( v1!=v2 );
  EXPECT_FALSE ( v1!=v3 );
  EXPECT_TRUE  ( not_equal(v1,v2) );
  EXPECT_FALSE ( not_equal(v1,v3) );
  EXPECT_TRUE  ( not_equal(v1,v2,.05) );
  EXPECT_FALSE ( not_equal(v1,v3,.05) );
  EXPECT_FALSE ( not_equal(v1,v2,.5) );
  EXPECT_FALSE ( not_equal(v1,v3,.5) );
}

TEST(Vector, BasicMath) {
  Vector3f v1(1,2,3), v2(2,4,4), v3(4.2f,-1.1f,-3.3f);

  EXPECT_EQ( Vector3f(-1, -2, -3)     , -v1 );

  EXPECT_EQ( Vector3f(3, 6, 7)        , elem_sum(v1, v2) );
  EXPECT_EQ( Vector3f(3, 6, 7)        , v1+v2 );
  EXPECT_EQ( Vector3f(2, 3, 4)        , elem_sum(v1, 1) );
  EXPECT_EQ( Vector3f(3, 4, 5)        , elem_sum(2, v1) );

  EXPECT_EQ( Vector3f(-1, -2, -1)     , elem_diff(v1, v2) );
  EXPECT_EQ( Vector3f(-1, -2, -1)     , v1-v2 );
  EXPECT_EQ( Vector3f(0, 1, 2)        , elem_diff(v1, 1) );
  EXPECT_EQ( Vector3f(1, 0, -1)       , elem_diff(2, v1) );

  EXPECT_EQ( Vector3f(2, 8, 12)       , elem_prod(v1, v2) );
  EXPECT_EQ( Vector3f(2, 4, 6)        , elem_prod(v1, 2) );
  EXPECT_EQ( Vector3f(2, 4, 6)        , v1*2 );
  EXPECT_EQ( Vector3f(3, 6, 9)        , elem_prod(3, v1) );
  EXPECT_EQ( Vector3f(3, 6, 9)        , 3*v1 );

  EXPECT_EQ( Vector3f(0.5, 0.5, 0.75) , elem_quot(v1, v2) );
  EXPECT_EQ( Vector3f(0.5, 1, 1.5)    , elem_quot(v1, 2) );
  EXPECT_EQ( Vector3f(0.5, 1, 1.5)    , v1/2 );
  EXPECT_EQ( Vector3f(3, 1.5, 1)      , elem_quot(3, v1) );

  EXPECT_EQ( Vector3f(4.2f,1.1f,3.3f) , abs(v3) );
}

TEST(Vector, ElemCompare) {
  typedef Vector<bool,3> Vector3b;

  Vector3f v1(1,2,3), v2(2,2,2);

  EXPECT_EQ( Vector3b(false,true,false) , elem_eq(v1,v2) );
  EXPECT_EQ( Vector3b(false,true,false) , elem_eq(v1,2) );
  EXPECT_EQ( Vector3b(false,true,false) , elem_eq(2,v1) );

  EXPECT_EQ( Vector3b(true,false,true)  , elem_neq(v1,v2) );
  EXPECT_EQ( Vector3b(true,false,true)  , elem_neq(v1,2) );
  EXPECT_EQ( Vector3b(true,false,true)  , elem_neq(2,v1) );

  EXPECT_EQ( Vector3b(true,false,false) , elem_lt(v1,v2) );
  EXPECT_EQ( Vector3b(true,false,false) , elem_lt(v1,2) );
  EXPECT_EQ( Vector3b(false,false,true) , elem_lt(2,v1) );

  EXPECT_EQ( Vector3b(false,false,true) , elem_gt(v1,v2) );
  EXPECT_EQ( Vector3b(false,false,true) , elem_gt(v1,2) );
  EXPECT_EQ( Vector3b(true,false,false) , elem_gt(2,v1) );

  EXPECT_EQ( Vector3b(true,true,false)  , elem_lte(v1,v2) );
  EXPECT_EQ( Vector3b(true,true,false)  , elem_lte(v1,2) );
  EXPECT_EQ( Vector3b(false,true,true)  , elem_lte(2,v1) );

  EXPECT_EQ( Vector3b(false,true,true)  , elem_gte(v1,v2) );
  EXPECT_EQ( Vector3b(false,true,true)  , elem_gte(v1,2) );
  EXPECT_EQ( Vector3b(true,true,false)  , elem_gte(2,v1) );
}

TEST(Vector, Norms) {
  Vector3f v(1,2,-1);

  EXPECT_EQ( 4, norm_1(v) );
  EXPECT_DOUBLE_EQ(2.44948974278317809819, norm_2(v));
  EXPECT_EQ( 6, norm_2_sqr(v) );
  EXPECT_EQ( 2, norm_inf(v) );
  EXPECT_EQ( 1u, index_norm_inf(v) );

  EXPECT_EQ( 2, sum(v) );
  EXPECT_EQ( -2, prod(v) );
}

TEST(Vector, Funcs) {
  Vector3f v1(1,2,3), v2(2,4,4);

  EXPECT_EQ(               22, dot_prod(v1,v2) );
  EXPECT_EQ( Vector3f(-4,2,0), cross_prod(v1,v2) );

  Vector3f vn = normalize(v1);
  EXPECT_FLOAT_EQ( 0.2672612419f, vn(0) );
  EXPECT_FLOAT_EQ( 0.5345224838f, vn(1) );
  EXPECT_FLOAT_EQ( 0.8017837257f, vn(2) );
}

TEST(Vector, Transpose) {
  Vector3f v1(1,2,3), v2(2,4,4), v3(4.2f,-1.1f,-3.3f);

  EXPECT_EQ( 22  , transpose(v1)*v2 );
  EXPECT_EQ( -22 , (-transpose(v1))*v2 );
  EXPECT_EQ( 58  , (transpose(v1)+transpose(v2))*v2 );
  EXPECT_EQ( -14 , (transpose(v1)-transpose(v2))*v2 );
  EXPECT_EQ( 44  , (2*transpose(v1))*v2 );
  EXPECT_EQ( 66  , (transpose(v1)*3)*v2 );
  EXPECT_EQ( 11  , (transpose(v1)/2)*v2 );
  EXPECT_NEAR( 16.3, abs(transpose(v3))*v1, 1e-6 );
}

TEST(Vector, Real) {
  Vector<std::complex<float>,3> v(std::complex<float>(1,2),
      std::complex<float>(2,3),
      std::complex<float>(3,4));

  EXPECT_EQ( 1, real(v)(0) );
  EXPECT_EQ( 2, real(v)(1) );
  EXPECT_EQ( 3, real(v)(2) );
}

TEST(Vector, Imag) {
  Vector<std::complex<float>,3> v(std::complex<float>(1,2),
      std::complex<float>(2,3),
      std::complex<float>(3,4));

  EXPECT_EQ( 2, imag(v)(0) );
  EXPECT_EQ( 3, imag(v)(1) );
  EXPECT_EQ( 4, imag(v)(2) );
}

TEST(Vector, IndexingIterator) {
  typedef Vector2 Vec;
  typedef math::IndexingVectorIterator<Vec> Iter;

  Vec v1(3,4);

  // Construct
  Iter i(v1, 0);
  EXPECT_EQ(3, *i);

  // Copy
  Iter j(i);
  EXPECT_EQ(3, *j);
  EXPECT_NE(&i, &j);
  EXPECT_TRUE(i == j);

  // Assign
  Iter k(v1, 0);
  k = j;
  EXPECT_EQ(3, *k);
  EXPECT_NE(&j, &k);
  EXPECT_TRUE(j == k);
}
