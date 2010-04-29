// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <vw/Math/Functors.h>
#include <boost/random.hpp>

using namespace vw;
using namespace vw::math;

static const double DELTA = 1e-5;

TEST(Accumulators, CDF_cauchy) {


  boost::minstd_rand random_gen(42u);
  random_gen.seed(0);
  boost::cauchy_distribution<double> cauchy(35,80);
  boost::variate_generator<boost::minstd_rand, boost::cauchy_distribution<double> > generator(random_gen, cauchy);

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
  boost::minstd_rand random_gen(42u);
  random_gen.seed(0);
  boost::triangle_distribution<double> triangular(10,60,80);
  boost::variate_generator<boost::minstd_rand, boost::triangle_distribution<double> > generator(random_gen, triangular);

  { // Default settings
    CDFAccumulator<double> cdf;
    for ( uint16 i = 0; i < 50000; i++ )
      cdf( generator() );

    EXPECT_NEAR( cdf.median(), 51.833, 1.0 );
    EXPECT_NEAR( cdf.first_quartile(), 39.6, 2.0 );
    EXPECT_NEAR( cdf.third_quartile(), 61.3, 2.0 );
    EXPECT_NEAR( cdf.approximate_mean(), 50, 1.0 );
  }
  { // More quantiles == more precision
    CDFAccumulator<double> cdf(2000,500);
    for ( uint16 i = 0; i < 50000; i++ )
      cdf( generator() );

    EXPECT_NEAR( cdf.median(), 51.833, 0.5 );
    EXPECT_NEAR( cdf.first_quartile(), 39.6, 1.0 );
    EXPECT_NEAR( cdf.third_quartile(), 61.3, 1.0 );
    EXPECT_NEAR( cdf.approximate_mean(), 50, 0.5 );
  }
}

TEST(Accumulators, Median) {
  MedianAccumulator<double> median;

  boost::minstd_rand random_gen(42u);
  random_gen.seed(0);
  boost::cauchy_distribution<double> cauchy(35,80);
  boost::variate_generator<boost::minstd_rand, boost::cauchy_distribution<double> > generator(random_gen, cauchy);

  for ( uint16 i = 0; i < 50000; i++ )
    median( generator() );

  EXPECT_NEAR( median.value(), 35.0, 1.5 );
}
