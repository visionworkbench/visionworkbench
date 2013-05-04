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


#include <gtest/gtest_VW.h>
#include <vw/Core/Functors.h>

using namespace vw;

template <class T1, class T2>
static bool is_of_type( T2 ) {
  return boost::is_same<T1,T2>::value;
}

TEST(Functors, Negation) {
  ArgNegationFunctor f;
  EXPECT_EQ( -1, f(1) );
  EXPECT_TRUE( is_of_type<int>( f(int(1)) ) );
  EXPECT_TRUE( is_of_type<float>( f(float(1)) ) );
  EXPECT_TRUE( is_of_type<double>( f(double(1)) ) );
}

// We use the value 1 instead of the default-constructed value here so
// that division by zero doesn't bite us when optimization is disabled.
#define TEST_BINARY_MATH_FUNCTOR(name,arg1,arg2,result)                                      \
  do {                                                                                       \
    ArgArg##name##Functor f1;                                                                \
    EXPECT_EQ( (result), f1((arg1),(arg2)) );                                                \
    EXPECT_TRUE( is_of_type<int>( f1(int(1),int(1)) ) );                                     \
    EXPECT_TRUE( is_of_type<float>( f1(int(1),float(1)) ) );                                 \
    EXPECT_TRUE( is_of_type<double>( f1(double(1),float(1)) ) );                             \
    EXPECT_EQ( (result), ValArg##name##Functor<double>(arg1)(arg2) );                        \
    EXPECT_TRUE( is_of_type<int>( ValArg##name##Functor<int>(int(1))(int(1)) ) );            \
    EXPECT_TRUE( is_of_type<float>( ValArg##name##Functor<float>(float(1))(int(1)) ) );      \
    EXPECT_TRUE( is_of_type<float>( ValArg##name##Functor<int>(int(1))(float(1)) ) );        \
    EXPECT_TRUE( is_of_type<double>( ValArg##name##Functor<float>(float(1))(double(1)) ) );  \
    EXPECT_TRUE( is_of_type<double>( ValArg##name##Functor<double>(double(1))(float(1)) ) ); \
    EXPECT_EQ( (result), ArgVal##name##Functor<double>(arg2)(arg1) );                        \
    EXPECT_TRUE( is_of_type<int>( ArgVal##name##Functor<int>(int(1))(int(1)) ) );            \
    EXPECT_TRUE( is_of_type<float>( ArgVal##name##Functor<float>(float(1))(int(1)) ) );      \
    EXPECT_TRUE( is_of_type<float>( ArgVal##name##Functor<int>(int(1))(float(1)) ) );        \
    EXPECT_TRUE( is_of_type<double>( ArgVal##name##Functor<float>(float(1))(double(1)) ) );  \
    EXPECT_TRUE( is_of_type<double>( ArgVal##name##Functor<double>(double(1))(float(1)) ) ); \
  } while(false)

TEST(Functors, Sum)        { TEST_BINARY_MATH_FUNCTOR(Sum,1,2,3); }
TEST(Functors, Difference) { TEST_BINARY_MATH_FUNCTOR(Difference,1,2,-1); }
TEST(Functors, Product)    { TEST_BINARY_MATH_FUNCTOR(Product,2,3,6); }
TEST(Functors, Quotient)   { TEST_BINARY_MATH_FUNCTOR(Quotient,6,3,2); }

#define TEST_BINARY_BOOL_FUNCTOR(name,arg1true,arg2true,arg1false,arg2false)             \
  do {                                                                                   \
    ArgArg##name##Functor f1;                                                            \
    EXPECT_TRUE(  f1((arg1true),(arg2true)) );                                           \
    EXPECT_FALSE( f1((arg1false),(arg2false)) );                                         \
    EXPECT_TRUE( is_of_type<bool>( f1(int(),int()) ) );                                  \
    EXPECT_TRUE( is_of_type<bool>( f1(int(),float()) ) );                                \
    EXPECT_TRUE( is_of_type<bool>( f1(double(),float()) ) );                             \
    EXPECT_TRUE(  ValArg##name##Functor<double>(arg1true)(arg2true) );                   \
    EXPECT_FALSE( ValArg##name##Functor<double>(arg1false)(arg2false) );                 \
    EXPECT_TRUE( is_of_type<bool>( ValArg##name##Functor<int>(int())(int()) ) );         \
    EXPECT_TRUE( is_of_type<bool>( ValArg##name##Functor<float>(float())(int()) ) );     \
    EXPECT_TRUE( is_of_type<bool>( ValArg##name##Functor<int>(int())(float()) ) );       \
    EXPECT_TRUE( is_of_type<bool>( ValArg##name##Functor<float>(float())(double()) ) );  \
    EXPECT_TRUE( is_of_type<bool>( ValArg##name##Functor<double>(double())(float()) ) ); \
    EXPECT_TRUE(  ArgVal##name##Functor<double>(arg2true)(arg1true) );                   \
    EXPECT_FALSE( ArgVal##name##Functor<double>(arg2false)(arg1false) );                 \
    EXPECT_TRUE( is_of_type<bool>( ArgVal##name##Functor<int>(int())(int()) ) );         \
    EXPECT_TRUE( is_of_type<bool>( ArgVal##name##Functor<float>(float())(int()) ) );     \
    EXPECT_TRUE( is_of_type<bool>( ArgVal##name##Functor<int>(int())(float()) ) );       \
    EXPECT_TRUE( is_of_type<bool>( ArgVal##name##Functor<float>(float())(double()) ) );  \
    EXPECT_TRUE( is_of_type<bool>( ArgVal##name##Functor<double>(double())(float()) ) ); \
  } while(false)

TEST(Functors, Equality)           { TEST_BINARY_BOOL_FUNCTOR(Equality,3,3,2,4); }
TEST(Functors, Inequality)         { TEST_BINARY_BOOL_FUNCTOR(Inequality,3,4,2,2); }
TEST(Functors, LessThan)           { TEST_BINARY_BOOL_FUNCTOR(LessThan,3,4,2,1); }
TEST(Functors, LessThanOrEqual)    { TEST_BINARY_BOOL_FUNCTOR(LessThanOrEqual,3,3,2,1); }
TEST(Functors, GreaterThan)        { TEST_BINARY_BOOL_FUNCTOR(GreaterThan,4,3,2,4); }
TEST(Functors, GreaterThanOrEqual) { TEST_BINARY_BOOL_FUNCTOR(GreaterThanOrEqual,3,3,2,4); }
