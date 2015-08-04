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
#include <vw/Math/Functors.h>
#include <boost/random.hpp>

using namespace vw;
using namespace vw::math;

static const double DELTA = 1e-5;

TEST(Accumulators, CDF_cauchy) {
  boost::mt19937 random_gen(42);
  boost::cauchy_distribution<double> cauchy(35,80);
  boost::variate_generator<boost::mt19937&,
    boost::cauchy_distribution<double> > generator(random_gen, cauchy);

  { // Default settings
    CDFAccumulator<double> cdf;
    for ( uint16 i = 0; i < 50000; i++ )
      cdf( generator() );

    EXPECT_NEAR( cdf.median(), 35.0, 2.0 );
    EXPECT_NEAR( cdf.first_quartile(), -45, 4.0 );
    EXPECT_NEAR( cdf.third_quartile(), 115, 4.0 );
  }
  { // More quantiles == more precision
    CDFAccumulator<double> cdf(2000,500);
    for ( uint16 i = 0; i < 50000; i++ )
      cdf( generator() );

    EXPECT_NEAR( cdf.median(), 35.0, 1.0 );
    EXPECT_NEAR( cdf.first_quartile(), -45, 2.0 );
    EXPECT_NEAR( cdf.third_quartile(), 115, 2.0 );
  }
}

TEST(Accumulators, CDF_triangular) {
  boost::mt19937 random_gen(42);
  boost::triangle_distribution<double> triangular(10,60,80);
  boost::variate_generator<boost::mt19937&,
    boost::triangle_distribution<double> > generator(random_gen, triangular);

  { // Default settings
    CDFAccumulator<double> cdf;
    for ( uint16 i = 0; i < 50000; i++ )
      cdf( generator() );

    EXPECT_NEAR( cdf.median(), 51.833, 1.0 );
    EXPECT_NEAR( cdf.first_quartile(), 39.6, 2.0 );
    EXPECT_NEAR( cdf.third_quartile(), 61.3, 2.0 );
    EXPECT_NEAR( cdf.approximate_mean(), 50, 1.0 );
    EXPECT_NEAR( cdf.approximate_mean(0.05), 50, 0.5 );
    EXPECT_NEAR( cdf.approximate_stddev(), 14.7196, 2.0 );
    EXPECT_NEAR( cdf.approximate_stddev(0.05), 14.7196, 1.0 );
  }
  { // More quantiles == more precision
    CDFAccumulator<double> cdf(2000,500);
    for ( uint16 i = 0; i < 50000; i++ )
      cdf( generator() );

    EXPECT_NEAR( cdf.median(), 51.833, 0.5 );
    EXPECT_NEAR( cdf.first_quartile(), 39.6, 1.0 );
    EXPECT_NEAR( cdf.third_quartile(), 61.3, 1.0 );
    EXPECT_NEAR( cdf.approximate_mean(), 50, 0.5 );
    EXPECT_NEAR( cdf.approximate_mean(0.05), 50, 0.25 );
    EXPECT_NEAR( cdf.approximate_stddev(), 14.7196, 3.0 );
    EXPECT_NEAR( cdf.approximate_stddev(0.05), 14.7196, 0.5 );
  }
}

TEST(Accumulators, Median) {
  MedianAccumulator<double> median;

  boost::mt19937 random_gen(42);
  boost::cauchy_distribution<double> cauchy(35,80);
  boost::variate_generator<boost::mt19937&,
    boost::cauchy_distribution<double> > generator(random_gen, cauchy);

  for ( uint16 i = 0; i < 50000; i++ )
    median( generator() );

  EXPECT_NEAR( median.value(), 35.0, 1.5 );
}

TEST(Accumulators, CDF_Merge ) {
  boost::mt19937 random_gen(42);
  boost::normal_distribution<double> norm1(0, 3);
  boost::normal_distribution<double> norm2(5, 3);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > generator1( random_gen, norm1 );
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > generator2( random_gen, norm2 );

  CDFAccumulator<double> cdf1, cdf2, cdf3;
  // 1 will only see gen1 .. but will later be merged with 2.
  // 2 will only see gen2
  // 3 is our control and sees both gen1 and gen2;
  for ( size_t i = 0; i < 50000; i++ ) {
    double sample1 = generator1(), sample2 = generator2();
    cdf1( sample1 );
    cdf2( sample2 );
    cdf3( sample1 );
    cdf3( sample2 );
  }

  cdf1.update();
  cdf2.update();
  cdf3.update();
  cdf1(cdf2);

  EXPECT_NEAR( cdf1.median(), cdf3.median(), 0.01 );
  EXPECT_NEAR( cdf1.first_quartile(), cdf3.first_quartile(), 0.01 );
  EXPECT_NEAR( cdf1.third_quartile(), cdf3.third_quartile(), 0.01 );
  EXPECT_NEAR( cdf1.approximate_mean(),
               cdf3.approximate_mean(), 0.01 );
}
