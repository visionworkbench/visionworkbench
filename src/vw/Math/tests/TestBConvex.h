// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestBConvex.h

#include <sstream>
#include <limits>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm> // std::sort

#include <boost/algorithm/string.hpp>

#include <cxxtest/TestSuite.h>
#include <vw/Math/Vector.h>
#include <vw/Math/BConvex.h>

using namespace vw;

bool check_bconvex_print(std::string a, std::string b)
{
  std::vector<std::string> split_a, split_b;
  std::vector<std::string>::iterator i;
  boost::split(split_a, a, boost::is_any_of(","));
  for (i = split_a.begin(); i != split_a.end(); i++)
    boost::trim(*i);
  std::sort(split_a.begin(), split_a.end());
  boost::split(split_b, b, boost::is_any_of(","));
  for (i = split_b.begin(); i != split_b.end(); i++)
    boost::trim(*i);
  std::sort(split_b.begin(), split_b.end());
  return (split_a == split_b);
}

class TestBConvex : public CxxTest::TestSuite
{
public:

  void test_bconvex()
  {
    std::vector<Vector3> v;

    BConvex c0;
    TS_ASSERT( c0.empty() );
    c0.grow(Vector3(0, 0, 0));
    TS_ASSERT( !c0.empty() );

    v.clear();
    // box with min (0,0,0) and max (1,1,1)
    v.push_back(Vector3(0, 0, 0));
    v.push_back(Vector3(0, 0, 1));
    v.push_back(Vector3(0, 1, 0));
    v.push_back(Vector3(0, 1, 1));
    v.push_back(Vector3(1, 0, 0));
    v.push_back(Vector3(1, 0, 1));
    v.push_back(Vector3(1, 1, 0));
    v.push_back(Vector3(1, 1, 1));
    BConvex c1(v);
    TS_ASSERT( c1.contains(Vector3(0.5, 0.5, 0.5)) );
    TS_ASSERT( c1.contains(Vector3(0.9, 0.9, 0.5)) );
    TS_ASSERT( !c1.contains(Vector3(0.9, 0.9, 1.5)) );

    c1.grow(Vector3(0.5, 0.5, 0.5));
    std::ostringstream os1;
    os1 << c1;
    TS_ASSERT( check_bconvex_print(os1.str(), std::string("-C >= -1, -B >= -1, -A >= -1, A >= 0, C >= 0, B >= 0")) || check_bconvex_print(os1.str(), std::string("-1*X + 0*Y + 0*Z + 1 >= 0, 0*X + -1*Y + 0*Z + 1 >= 0, 0*X + 0*Y + -1*Z + 1 >= 0, 0*X + 0*Y + 1*Z + 0 >= 0, 0*X + 1*Y + 0*Z + 0 >= 0, 1*X + 0*Y + 0*Z + 0 >= 0")) );
    std::ostringstream os2;
    os2 << (c1 * 4);
    TS_ASSERT( check_bconvex_print(os2.str(), std::string("-C >= -4, -B >= -4, -A >= -4, A >= 0, C >= 0, B >= 0")) || check_bconvex_print(os2.str(), std::string("-1*X + 0*Y + 0*Z + 4 >= 0, 0*X + -1*Y + 0*Z + 4 >= 0, 0*X + 0*Y + -1*Z + 4 >= 0, 0*X + 0*Y + 1*Z + 0 >= 0, 0*X + 1*Y + 0*Z + 0 >= 0, 1*X + 0*Y + 0*Z + 0 >= 0")) );
    std::ostringstream os3;
    os3 << (4 * c1);
    TS_ASSERT( check_bconvex_print(os3.str(), std::string("-C >= -4, -B >= -4, -A >= -4, A >= 0, C >= 0, B >= 0")) || check_bconvex_print(os3.str(), std::string("-1*X + 0*Y + 0*Z + 4 >= 0, 0*X + -1*Y + 0*Z + 4 >= 0, 0*X + 0*Y + -1*Z + 4 >= 0, 0*X + 0*Y + 1*Z + 0 >= 0, 0*X + 1*Y + 0*Z + 0 >= 0, 1*X + 0*Y + 0*Z + 0 >= 0")) );
    std::ostringstream os4;
    os4 << ((c1 * 4) / 4);
    TS_ASSERT( check_bconvex_print(os4.str(), std::string("-C >= -1, -B >= -1, -A >= -1, A >= 0, C >= 0, B >= 0")) || check_bconvex_print(os4.str(), std::string("-1*X + 0*Y + 0*Z + 1 >= 0, 0*X + -1*Y + 0*Z + 1 >= 0, 0*X + 0*Y + -1*Z + 1 >= 0, 0*X + 0*Y + 1*Z + 0 >= 0, 0*X + 1*Y + 0*Z + 0 >= 0, 1*X + 0*Y + 0*Z + 0 >= 0")) );
    std::ostringstream os5;
    os5 << (c1 + Vector3(1, -2, 3));
    TS_ASSERT( check_bconvex_print(os5.str(), std::string("-C >= -4, -B >= 1, -A >= -2, A >= 1, C >= 3, B >= -2")) || check_bconvex_print(os5.str(), std::string("-1*X + 0*Y + 0*Z + 2 >= 0, 0*X + -1*Y + 0*Z + -1 >= 0, 0*X + 0*Y + -1*Z + 4 >= 0, 0*X + 0*Y + 1*Z + -3 >= 0, 0*X + 1*Y + 0*Z + 2 >= 0, 1*X + 0*Y + 0*Z + -1 >= 0")) );
    std::ostringstream os6;
    os6 << ((c1 + Vector3(1, -2, 3)) - Vector3(1, -2, 3));
    TS_ASSERT( check_bconvex_print(os6.str(), std::string("-C >= -1, -B >= -1, -A >= -1, A >= 0, C >= 0, B >= 0")) || check_bconvex_print(os6.str(), std::string("-1*X + 0*Y + 0*Z + 1 >= 0, 0*X + -1*Y + 0*Z + 1 >= 0, 0*X + 0*Y + -1*Z + 1 >= 0, 0*X + 0*Y + 1*Z + 0 >= 0, 0*X + 1*Y + 0*Z + 0 >= 0, 1*X + 0*Y + 0*Z + 0 >= 0")) );

    v.clear();
    // box with min (-2,-2,-2) and max (2,2,2)
    v.push_back(Vector3(-2, -2, -2));
    v.push_back(Vector3(-2, -2, 2));
    v.push_back(Vector3(-2, 2, -2));
    v.push_back(Vector3(-2, 2, 2));
    v.push_back(Vector3(2, -2, -2));
    v.push_back(Vector3(2, -2, 2));
    v.push_back(Vector3(2, 2, -2));
    v.push_back(Vector3(2, 2, 2));
    BConvex c2(v);
    TS_ASSERT( c2.contains(c1) );
    TS_ASSERT( !c1.contains(c2) );
    TS_ASSERT( c2.intersects(c1) );
    TS_ASSERT( c1.intersects(c2) );

    v.clear();
    // pyramid
    v.push_back(Vector3(-1.5, -1.5, -1.5));
    v.push_back(Vector3(-1.5, -1.5, -1.0));
    v.push_back(Vector3(-1.5, -1.0, -1.5));
    v.push_back(Vector3(-1.5, -1.0, -1.0));
    v.push_back(Vector3(-3, -1.25, -1.25));
    BConvex c3(v);
    TS_ASSERT( c3.contains(Vector3(-2, -1.1, -1.1)) );
    TS_ASSERT( !c3.contains(Vector3(-2, -1.01, -1.01)) );
    TS_ASSERT( !c1.contains(c3) );
    TS_ASSERT( !c3.contains(c1) );
    TS_ASSERT( !c2.contains(c3) );
    TS_ASSERT( !c3.contains(c2) );
    TS_ASSERT( !c1.intersects(c3) );
    TS_ASSERT( !c3.intersects(c1) );
    TS_ASSERT( c2.intersects(c3) );
    TS_ASSERT( c3.intersects(c2) );

    BConvex c4(v);
    TS_ASSERT( c3 == c4 );
    TS_ASSERT( !(c3 != c4) );

    BConvex c5(c4);
    TS_ASSERT( c3 == c5 );
    TS_ASSERT( !(c3 != c5) );

    c5 = c2;
    TS_ASSERT( c2 == c5 );
    TS_ASSERT( !(c2 != c5) );
    c5.grow(Vector3(-5, -5, -5));
    TS_ASSERT( !(c2 == c5) );
    TS_ASSERT( c2 != c5 );

    c5.crop(c2);
    TS_ASSERT( c2 == c5 );
    TS_ASSERT( !(c2 != c5) );

    c5 = c3;
    TS_ASSERT( !c5.contains(Vector3(0.5, 0.5, 0.5)) );
    TS_ASSERT( !c5.contains(Vector3(0.9, 0.9, 0.5)) );
    TS_ASSERT( !c5.contains(Vector3(0.9, 0.9, 1.5)) );
    c5.grow(c1);
    TS_ASSERT( c5.contains(Vector3(0.5, 0.5, 0.5)) );
    TS_ASSERT( c5.contains(Vector3(0.9, 0.9, 0.5)) );
    TS_ASSERT( !c5.contains(Vector3(0.9, 0.9, 1.5)) );

    BConvex c6;
    c6.grow(Vector3i(-2, 0, 0));
    c6.grow(Vector3i(-2, 0, 2));
    c6.grow(Vector3i(-2, 2, 0));
    c6.grow(Vector3i(-2, 2, 2));
    c6.grow(Vector3i(2, 0, 0));
    c6.grow(Vector3i(2, 0, 2));
    c6.grow(Vector3i(2, 2, 0));
    c6.grow(Vector3i(2, 2, 2));
    TS_ASSERT_EQUALS( c6.center(), Vector3i(0, 1, 1) );
    TS_ASSERT_EQUALS( c6.size(), 2*std::sqrt(6) );

    BConvex c7;
    c7.grow(c6);
    TS_ASSERT( c7 == c6 );

    BConvex c8;
  }

}; // class TestBConvex
