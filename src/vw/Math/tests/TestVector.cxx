// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


// TestVector.h
#include <test/Helpers.h>
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

  // set_all
  v3.set_all(6);
  EXPECT_EQ( 6, v3(0) );
  EXPECT_EQ( 6, v3(1) );
  EXPECT_EQ( 6, v3(2) );
  EXPECT_EQ( 6, v3(3) );

  // set_size()
  EXPECT_THROW(v.set_size(3), ArgumentErr);
  ASSERT_NO_THROW(v.set_size(4));

  // Iterators
  EXPECT_EQ(&(*(v.begin())),&(v(0)));
  EXPECT_EQ(&(*(v.begin()+1)),&(v(1)));
  EXPECT_EQ(v.end(),v.begin()+4);
}

TEST(Vector, Static_Construction) {

  size_t value = 0;
  Vector3 v(value, value, value); // Should not compile
  Vector2 v2(0, 0);
  //Vector3 v3(v2); // Should not compile
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

  // set_all
  v.set_all(6);
  EXPECT_EQ( 6, v(0) );
  EXPECT_EQ( 6, v(1) );
  EXPECT_EQ( 6, v(2) );

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
  // Check writing a vector to a stream
  std::stringstream ss;
  Vector3f v1(1,2,3);
  ss << v1;
  EXPECT_EQ( ss.str(), "Vector3(1,2,3)" );
  ss.str("");
  Vector<float> v2=v1;
  ss << v2;
  EXPECT_EQ( ss.str(), "Vector3(1,2,3)" );
  
  // Check read from a stream
  Vector<float32,0> vread;
  ss >> vread;
  EXPECT_EQ(v1, vread);
  
  ss.str("");
  vread = Vector3f(.1336185e-8,-0.5226175e-12, 0);
  const size_t ACCURATE_DIGITS = 17; // = std::numeric_limits<double>::max_digits10
  ss << std::setprecision(ACCURATE_DIGITS) << vread;
  Vector3f vread3;
  ss >> vread3;
  EXPECT_EQ(vread, vread3);
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

TEST(Vector, Polyval) {
  Vector3 v(3,2,1);
  EXPECT_EQ( 86,  polyval(v,5) );
  EXPECT_EQ( 162, polyval(v,7) );
  EXPECT_EQ( 262, polyval(v,9) );
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

  EXPECT_VECTOR_EQ( real(v), Vector3f(1,2,3) );
}

TEST(Vector, Imag) {
  Vector<std::complex<float>,3> v(std::complex<float>(1,2),
      std::complex<float>(2,3),
      std::complex<float>(3,4));

  EXPECT_VECTOR_EQ( imag(v), Vector3f(2,3,4) );
}

TEST(Vector, Floor) {
  Vector3 v(1.4,25.3,-13.7);
  Vector3 o = floor(v);

  EXPECT_VECTOR_EQ( o, Vector3(1,25,-14) );
}

TEST(Vector, Ceil) {
  Vector3 v(1.4,25.3,-13.7);
  Vector3 o = ceil(v);

  EXPECT_VECTOR_EQ( o, Vector3(2,26,-13) );
}

// Make sure that the set_all function works on basic types
TEST(Vector, set_all) {

  Vector2 v2(1, 2);
  set_all(v2, 3);
  EXPECT_VECTOR_EQ( v2, Vector2(3,3) );
  
  Vector<float> v4(4);
  set_all(v4, 3);
  EXPECT_VECTOR_EQ( v4, Vector4(3,3,3,3) );
  
  double d=4;
  set_all(d, 3.0);
  EXPECT_EQ(d, 3);
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

TEST(Vector, SmallIntArithmetic) {
  {
    Vector3 small(1e-5,1e-5,2e-5);
    int64 scalar = 10000;
    EXPECT_VECTOR_NEAR( Vector3(.1,.1,.2),
                        (small*scalar), 1e-6 );
    EXPECT_VECTOR_NEAR( Vector3(.1,.1,.2),
                        (scalar*small), 1e-6 );
    EXPECT_VECTOR_NEAR( Vector3(.1,.1,.2),
                        double(scalar)*small, 1e-6 );
    EXPECT_VECTOR_NEAR( Vector3(.1,.1,.2),
                        small*double(scalar), 1e-6 );
    // These evaluate differently somehow
    const Vector3 iresult = scalar*small;
    Vector3 fresult = double(scalar)*small;
    EXPECT_VECTOR_NEAR( Vector3(.1,.1,.2), iresult, 1e-6 );
    EXPECT_VECTOR_NEAR( Vector3(.1,.1,.2), fresult, 1e-6 );
  }
  {
    const Vector3 small(0.1,0.1,0.2);
    int64 scalar = 10;
    EXPECT_VECTOR_NEAR( Vector3(.01,.01,.02),
                        small/scalar, 1e-6 );
    EXPECT_VECTOR_NEAR( Vector3(.01,.01,.02),
                        small/double(scalar),1e-6 );
    Vector3 iresult = small/scalar;
    Vector3 fresult = small/double(scalar);
    EXPECT_VECTOR_NEAR( Vector3(.01,.01,.02), iresult, 1e-6 );
    EXPECT_VECTOR_NEAR( Vector3(.01,.01,.02), fresult, 1e-6 );
  }
}
