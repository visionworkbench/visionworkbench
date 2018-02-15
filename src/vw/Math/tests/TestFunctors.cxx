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
#include <boost/random.hpp>
#include <vw/config.h>
#include <vw/Math/Functors.h>

using namespace vw;
using namespace vw::math;

static const double DELTA = 1e-5;

template <class T1, class T2>
static bool is_of_type( T2 ) {
  return boost::is_same<T1,T2>::value;
}

TEST(Functors, Real) {
  ArgRealFunctor f;
  EXPECT_EQ( f(2.0), 2.0 );
  EXPECT_EQ( f(std::complex<double>(2.0,3.0)), 2.0 );
  EXPECT_TRUE( is_of_type<float>( f(float()) ) );
  EXPECT_TRUE( is_of_type<double>( f(double()) ) );
  EXPECT_TRUE( is_of_type<long double>( f((long double)(0)) ) );
  EXPECT_TRUE( is_of_type<int>( f(int()) ) );
  EXPECT_TRUE( is_of_type<float>( f(std::complex<float> ()) ) );
  EXPECT_TRUE( is_of_type<double>( f(std::complex<double> ()) ) );
  EXPECT_TRUE( is_of_type<long double>( f(std::complex<long double> ()) ) );
  EXPECT_TRUE( is_of_type<int>( f(std::complex<int> ()) ) );
}

TEST(Functors, Imag) {
  ArgImagFunctor f;
  EXPECT_EQ( f(2.0), 0.0 );
  EXPECT_EQ( f(std::complex<double>(2.0,3.0)), 3.0 );
  EXPECT_TRUE( is_of_type<float>( f(float()) ) );
  EXPECT_TRUE( is_of_type<double>( f(double()) ) );
  EXPECT_TRUE( is_of_type<long double>( f((long double)(0)) ) );
  EXPECT_TRUE( is_of_type<int>( f(int()) ) );
  EXPECT_TRUE( is_of_type<float>( f(std::complex<float>()) ) );
  EXPECT_TRUE( is_of_type<double>( f(std::complex<double>()) ) );
  EXPECT_TRUE( is_of_type<long double>( f(std::complex<long double>()) ) );
  EXPECT_TRUE( is_of_type<int>( f(std::complex<int>()) ) );
}

TEST(Functors, Abs) {
  ArgAbsFunctor f;
  EXPECT_EQ( f(2.0), 2.0 );
  EXPECT_EQ( f(-2.0), 2.0 );
  EXPECT_DOUBLE_EQ( f(std::complex<double>(3.0,4.0)), 5.0 );
  EXPECT_TRUE( is_of_type<float>( f(float()) ) );
  EXPECT_TRUE( is_of_type<double>( f(double()) ) );
  EXPECT_TRUE( is_of_type<long double>( f((long double)(0)) ) );
  EXPECT_TRUE( is_of_type<int>( f(int()) ) );
  EXPECT_TRUE( is_of_type<float>( f(std::complex<float>()) ) );
  EXPECT_TRUE( is_of_type<double>( f(std::complex<double>()) ) );
  EXPECT_TRUE( is_of_type<long double>( f(std::complex<long double>()) ) );
}

TEST(Functors, Conj) {
  ArgConjFunctor f;
  EXPECT_EQ( f(2.0), 2.0 );
  EXPECT_EQ( f(std::complex<double>(1.0,2.0)), std::complex<double>(1.0,-2.0) );
  EXPECT_TRUE( is_of_type<float>( f(float()) ) );
  EXPECT_TRUE( is_of_type<double>( f(double()) ) );
  EXPECT_TRUE( is_of_type<long double>( f((long double)(0)) ) );
  EXPECT_TRUE( is_of_type<int>( f(int()) ) );
  EXPECT_TRUE( is_of_type<std::complex<float> >( f(std::complex<float>()) ) );
  EXPECT_TRUE( is_of_type<std::complex<double> >( f(std::complex<double>()) ) );
  EXPECT_TRUE( is_of_type<std::complex<long double> >( f(std::complex<long double>()) ) );
  EXPECT_TRUE( is_of_type<std::complex<int> >( f(std::complex<int>()) ) );
}

TEST(Functors, Median){
  MedianAccumulator<double> V;
  V(8);
  V(9);
  V(3);
  V(5);

  // The median better be 6.5
  EXPECT_TRUE( V.value() == 6.5 );
}

TEST(Functors, StdDev){
  StdDevAccumulator<double> V;
  V(8);
  V(9);
  V(3);
  V(5);

  EXPECT_TRUE( V.value() == 2.3848480035423640366 );
}

TEST(Functors, DestructiveMedian){

  std::vector<double> V;
  V.push_back(8);
  V.push_back(9);
  V.push_back(3);
  V.push_back(5);

  // The median better be 6.5
  EXPECT_TRUE( destructive_median(V) == 6.5 );
}

TEST(Functors, DestructiveNmad){

  std::vector<double> V;
  V.push_back(8);
  V.push_back(9);
  V.push_back(3);
  V.push_back(5);

  EXPECT_TRUE( destructive_nmad(V) == 2.9652 );
}

TEST(Functors, DestructivePercentile){

  std::vector<double> V, W;
  V.push_back(8);
  V.push_back(9);
  V.push_back(3);
  V.push_back(5);

  // Start with new W each time, as it will be messed up

  W = V; EXPECT_TRUE( destructive_percentile(V,   0) == 3 );
  W = V; EXPECT_TRUE( destructive_percentile(V,  20) == 3 );
  W = V; EXPECT_TRUE( destructive_percentile(V,  25) == 3 );
  W = V; EXPECT_TRUE( destructive_percentile(V,  26) == 5 );
  W = V; EXPECT_TRUE( destructive_percentile(V,  45) == 5 );
  W = V; EXPECT_TRUE( destructive_percentile(V,  50) == 5 );
  W = V; EXPECT_TRUE( destructive_percentile(V,  51) == 8 );
  W = V; EXPECT_TRUE( destructive_percentile(V,  70) == 8 );
  W = V; EXPECT_TRUE( destructive_percentile(V,  75) == 8 );
  W = V; EXPECT_TRUE( destructive_percentile(V,  76) == 9 );
  W = V; EXPECT_TRUE( destructive_percentile(V, 100) == 9 );
  
}

#define TEST_UNARY_MATH_FUNCTOR(func,arg,result)                        \
  do {                                                                                                    \
    Arg##func##Functor f;                                                                                 \
    EXPECT_NEAR( (result), f((float)(arg)),       DELTA );                                                        \
    EXPECT_NEAR( (result), f((double)(arg)),      DELTA );                                                       \
    EXPECT_NEAR( (result), f((long double)(arg)), DELTA );                                                  \
    EXPECT_TRUE( is_of_type<float>(f((float)(arg))) );                                                      \
    EXPECT_TRUE( is_of_type<double>(f((double)(arg))) );                                                    \
    EXPECT_TRUE( is_of_type<long double>(f((long double)(arg))) );                                          \
    EXPECT_TRUE( is_of_type<double>(f((int)(arg))) );                                                       \
  } while(false)

#define TEST_BINARY_MATH_FUNCTOR(func,arg1,arg2,result)                                                   \
  do {                                                                                                    \
    ArgArg##func##Functor f;                                                                              \
    EXPECT_NEAR( (result), f((float)(arg1),(float)(arg2)),             DELTA );                                         \
    EXPECT_NEAR( (result), f((double)(arg1),(double)(arg2)),           DELTA );                                       \
    EXPECT_NEAR( (result), f((long double)(arg1),(long double)(arg2)), DELTA );                             \
    EXPECT_TRUE( is_of_type<float>(f((float)(arg1),(float)(arg2))) );                                       \
    EXPECT_TRUE( is_of_type<double>(f((double)(arg1),(double)(arg2))) );                                    \
    EXPECT_TRUE( is_of_type<long double>(f((long double)(arg1),(long double)(arg2))) );                     \
    EXPECT_TRUE( is_of_type<double>( f((float)(arg1),(double)(arg2))) );                                    \
    EXPECT_TRUE( is_of_type<long double>( f((float)(arg1),(long double)(arg2))) );                          \
    EXPECT_TRUE( is_of_type<double>( f((int)(arg1),(int)(arg2))) );                                         \
    EXPECT_TRUE( is_of_type<float>( f((int)(arg1),(float)(arg2))) );                                        \
    EXPECT_NEAR( (result), ArgVal##func##Functor<float>(arg2)((float)(arg1)),             DELTA );                      \
    EXPECT_NEAR( (result), ArgVal##func##Functor<double>(arg2)((double)(arg1)),           DELTA );                    \
    EXPECT_NEAR( (result), ArgVal##func##Functor<long double>(arg2)((long double)(arg1)), DELTA );          \
    EXPECT_TRUE( is_of_type<float>(ArgVal##func##Functor<float>(arg2)((float)(arg1))) );                    \
    EXPECT_TRUE( is_of_type<double>(ArgVal##func##Functor<double>(arg2)((double)(arg1))) );                 \
    EXPECT_TRUE( is_of_type<long double>(ArgVal##func##Functor<long double>(arg2)((long double)(arg1))) );  \
    EXPECT_TRUE( is_of_type<double>( ArgVal##func##Functor<double>(arg2)((float)(arg1))) );                 \
    EXPECT_TRUE( is_of_type<long double>( ArgVal##func##Functor<long double>(arg2)((float)(arg1))) );       \
    EXPECT_TRUE( is_of_type<double>( ArgVal##func##Functor<int>((int)(arg2))((int)(arg1))) );               \
    EXPECT_TRUE( is_of_type<float>( ArgVal##func##Functor<float>(arg2)((int)(arg1))) );                     \
    EXPECT_NEAR( (result), ValArg##func##Functor<float>(arg1)((float)(arg2)),             DELTA );                      \
    EXPECT_NEAR( (result), ValArg##func##Functor<double>(arg1)((double)(arg2)),           DELTA );                    \
    EXPECT_NEAR( (result), ValArg##func##Functor<long double>(arg1)((long double)(arg2)), DELTA );          \
    EXPECT_TRUE( is_of_type<float>(ValArg##func##Functor<float>(arg1)((float)(arg2))) );                    \
    EXPECT_TRUE( is_of_type<double>(ValArg##func##Functor<double>(arg1)((double)(arg2))) );                 \
    EXPECT_TRUE( is_of_type<long double>(ValArg##func##Functor<long double>(arg1)((long double)(arg2))) );  \
    EXPECT_TRUE( is_of_type<double>( ValArg##func##Functor<float>(arg1)((double)(arg2))) );                 \
    EXPECT_TRUE( is_of_type<long double>( ValArg##func##Functor<float>(arg1)((long double)(arg2))) );       \
    EXPECT_TRUE( is_of_type<double>( ValArg##func##Functor<int>((int)(arg1))((int)(arg2))) );               \
    EXPECT_TRUE( is_of_type<float>( ValArg##func##Functor<int>((int)arg1)((float)(arg2))) );                \
  } while(false)

TEST(Functors, Acos)  { TEST_UNARY_MATH_FUNCTOR(Acos,  0.5, 1.04719755);   }
TEST(Functors, Asin)  { TEST_UNARY_MATH_FUNCTOR(Asin,  0.5, 0.523598776); }
TEST(Functors, Atan)  { TEST_UNARY_MATH_FUNCTOR(Atan,  1.0, 0.785398163); }
TEST(Functors, Cos)   { TEST_UNARY_MATH_FUNCTOR(Cos,   1.0, 0.540302306);  }
TEST(Functors, Sin)   { TEST_UNARY_MATH_FUNCTOR(Sin,   1.0, 0.841471);  }
TEST(Functors, Tan)   { TEST_UNARY_MATH_FUNCTOR(Tan,   1.0, 1.55741);   }
TEST(Functors, Cosh)  { TEST_UNARY_MATH_FUNCTOR(Cosh,  1.0, 1.54308);  }
TEST(Functors, Sinh)  { TEST_UNARY_MATH_FUNCTOR(Sinh,  1.0, 1.1752);   }
TEST(Functors, Tanh)  { TEST_UNARY_MATH_FUNCTOR(Tanh,  1.0, 0.761594); }
TEST(Functors, Exp)   { TEST_UNARY_MATH_FUNCTOR(Exp,   1.0, 2.718281);  }
TEST(Functors, Log)   { TEST_UNARY_MATH_FUNCTOR(Log,   2.0, 0.693147);  }
TEST(Functors, Log10) { TEST_UNARY_MATH_FUNCTOR(Log10, 2.0, 0.30103); }
TEST(Functors, Sqrt)  { TEST_UNARY_MATH_FUNCTOR(Sqrt,  2.0, 1.41421);  }
TEST(Functors, Ceil)  { TEST_UNARY_MATH_FUNCTOR(Ceil,  1.5, 2.0);
                        TEST_UNARY_MATH_FUNCTOR(Ceil, -1.5, -1.0);    }
TEST(Functors, Floor) { TEST_UNARY_MATH_FUNCTOR(Floor, 1.5, 1.0);
                        TEST_UNARY_MATH_FUNCTOR(Floor,-1.5, -2.0);   }

TEST(Functors, Atan2) { TEST_BINARY_MATH_FUNCTOR(Atan2,2.0,1.0,1.10715); }
TEST(Functors, Pow)   { TEST_BINARY_MATH_FUNCTOR(Pow,3.0,2.0,9.0);  }

#ifndef WIN32
TEST(Functors, Acosh) { TEST_UNARY_MATH_FUNCTOR(Acosh,1.5,0.962424); }
TEST(Functors, Asinh) { TEST_UNARY_MATH_FUNCTOR(Asinh,1.0,0.881374); }
TEST(Functors, Atanh) { TEST_UNARY_MATH_FUNCTOR(Atanh,0.5,0.549306); }

#ifdef VW_HAVE_EXP2
TEST(Functors, Exp2)  { TEST_UNARY_MATH_FUNCTOR(Exp2,1.0,2.0); }
#else
TEST(Functors, DISABLED_Exp2)  { }
#endif

#ifdef VW_HAVE_LOG2
TEST(Functors, Log2) { TEST_UNARY_MATH_FUNCTOR(Log2,2.0,1.0); }
#else
TEST(Functors, DISABLED_Log2) { }
#endif

#ifdef VW_HAVE_TGAMMA
TEST(Functors, Tgamma) { TEST_UNARY_MATH_FUNCTOR(Tgamma,1.5,0.886227); }
#else
TEST(Functors, DISABLED_Tgamma) { }
#endif

TEST(Functors, Expm1)    { TEST_UNARY_MATH_FUNCTOR(Expm1,1.0,1.718281);  }
TEST(Functors, Log1p)    { TEST_UNARY_MATH_FUNCTOR(Log1p,1.0,0.693147);  }
TEST(Functors, Cbrt)     { TEST_UNARY_MATH_FUNCTOR(Cbrt,2.0,1.25992);    }
TEST(Functors, Erf)      { TEST_UNARY_MATH_FUNCTOR(Erf,1.0,0.842701);    }
TEST(Functors, Erfc)     { TEST_UNARY_MATH_FUNCTOR(Erfc,1.0,0.157299);   }
TEST(Functors, Lgamma)   { TEST_UNARY_MATH_FUNCTOR(Lgamma,2.5,0.284683); }
TEST(Functors, Round)    { TEST_UNARY_MATH_FUNCTOR(Round,1.4,1.0);
                           TEST_UNARY_MATH_FUNCTOR(Round,1.5,2.0);       }
TEST(Functors, Trunc)    { TEST_UNARY_MATH_FUNCTOR(Trunc,1.5,1.0);
                           TEST_UNARY_MATH_FUNCTOR(Trunc,-1.5,-1.0);     }

TEST(Functors, Hypot)    { TEST_BINARY_MATH_FUNCTOR(Hypot,2.0,1.0,2.23607); }
TEST(Functors, Copysign) { TEST_BINARY_MATH_FUNCTOR(Copysign,3.0,-2.0,-3.0);
                           TEST_BINARY_MATH_FUNCTOR(Copysign,3.0,2.0,3.0); }
TEST(Functors, Fdim)     { TEST_BINARY_MATH_FUNCTOR(Fdim,3.0,2.0,1.0);
                           TEST_BINARY_MATH_FUNCTOR(Fdim,2.0,3.0,0.0); }

#endif


TEST(Functiors, Median) {
  MedianAccumulator<double> median;

  boost::mt19937 random_gen(42);
  boost::cauchy_distribution<double> cauchy(35,80);
  boost::variate_generator<boost::mt19937&,
    boost::cauchy_distribution<double> > generator(random_gen, cauchy);

  for ( uint16 i = 0; i < 50000; i++ )
    median( generator() );

  EXPECT_NEAR( median.value(), 35.0, 1.5 );
}

