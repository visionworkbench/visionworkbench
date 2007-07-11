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

  void test_rational()
  {
    using math::bconvex_rational::operator<<;

    math::bconvex_rational::Rational r;
    std::ostringstream os, os1;

    math::bconvex_rational::RationalFuncs<int32> f1;
    f1.set(0, 1, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("0/1") );
    TS_ASSERT_EQUALS( f1.get(r), 0 );
    f1.set(1, 1, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("1/1") );
    TS_ASSERT_EQUALS( f1.get(r), 1 );
    f1.set(10, 10, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("10/1") );
    TS_ASSERT_EQUALS( f1.get(r), 10 );
    f1.set(-11, 11, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("-11/1") );
    TS_ASSERT_EQUALS( f1.get(r), -11 );
    f1.set(std::numeric_limits<int32>::max(), std::numeric_limits<int32>::max(), r);
    os.str("");
    os << r;
    os1.str("");
    os1 << std::numeric_limits<int32>::max() << std::string("/1");
    TS_ASSERT_EQUALS( os.str(), os1.str() );
    TS_ASSERT_EQUALS( f1.get(r), std::numeric_limits<int32>::max() );
    f1.set(std::numeric_limits<int32>::min(), -std::numeric_limits<int32>::min(), r);
    os.str("");
    os << r;
    os1.str("");
    os1 << std::numeric_limits<int32>::min() << std::string("/1");
    TS_ASSERT_EQUALS( os.str(), os1.str() );
    TS_ASSERT_EQUALS( f1.get(r), std::numeric_limits<int32>::min() );

    math::bconvex_rational::RationalFuncs<uint32> f2;
    f2.set(0, 1, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("0/1") );
    TS_ASSERT_EQUALS( f2.get(r), 0 );
    f2.set(1, 1, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("1/1") );
    TS_ASSERT_EQUALS( f2.get(r), 1 );
    f2.set(10, 10, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("10/1") );
    TS_ASSERT_EQUALS( f2.get(r), 10 );
    f2.set(std::numeric_limits<uint32>::max(), std::numeric_limits<uint32>::max(), r);
    os.str("");
    os << r;
    os1.str("");
    os1 << std::numeric_limits<uint32>::max() << std::string("/1");
    TS_ASSERT_EQUALS( os.str(), os1.str() );
    TS_ASSERT_EQUALS( f2.get(r), std::numeric_limits<uint32>::max() );
    f2.set(std::numeric_limits<uint32>::min(), -std::numeric_limits<uint32>::min(), r);
    os.str("");
    os << r;
    os1.str("");
    os1 << std::numeric_limits<uint32>::min() << std::string("/1");
    TS_ASSERT_EQUALS( os.str(), os1.str() );
    TS_ASSERT_EQUALS( f2.get(r), std::numeric_limits<uint32>::min() );

    math::bconvex_rational::RationalFuncs<int64> f3;
    f3.set(0, 1, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("0/1") );
    TS_ASSERT_EQUALS( f3.get(r), 0 );
    f3.set(1, 1, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("1/1") );
    TS_ASSERT_EQUALS( f3.get(r), 1 );
    f3.set(10, 10, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("10/1") );
    TS_ASSERT_EQUALS( f3.get(r), 10 );
    f3.set(-11, 11, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("-11/1") );
    TS_ASSERT_EQUALS( f3.get(r), -11 );
    f3.set(std::numeric_limits<int32>::max(), std::numeric_limits<int32>::max(), r);
    os.str("");
    os << r;
    os1.str("");
    os1 << std::numeric_limits<int32>::max() << std::string("/1");
    TS_ASSERT_EQUALS( os.str(), os1.str() );
    TS_ASSERT_EQUALS( f3.get(r), std::numeric_limits<int32>::max() );
    //std::cout << std::numeric_limits<int32>::min() << -std::numeric_limits<int32>::min() << std::endl;
    f3.set((uint64)std::numeric_limits<int32>::min(), -(uint64)std::numeric_limits<int32>::min(), r);
    os.str("");
    os << r;
    os1.str("");
    os1 << std::numeric_limits<int32>::min() << std::string("/1");
    TS_ASSERT_EQUALS( os.str(), os1.str() );
    TS_ASSERT_EQUALS( f3.get(r), std::numeric_limits<int32>::min() );
    //f3.set(std::numeric_limits<int64>::max(), std::numeric_limits<int64>::max(), r);
    //os.str("");
    //os << r;
    //os1.str("");
    //os1 << std::numeric_limits<int64>::max() << std::string("/1");
    //TS_ASSERT_EQUALS( os.str(), os1.str() );
    //TS_ASSERT_EQUALS( f3.get(r), std::numeric_limits<int64>::max() );
    //f3.set(std::numeric_limits<int64>::min(), -std::numeric_limits<int64>::min(), r);
    //os.str("");
    //os << r;
    //os1.str("");
    //os1 << std::numeric_limits<int64>::min() << std::string("/1");
    //TS_ASSERT_EQUALS( os.str(), os1.str() );
    //TS_ASSERT_EQUALS( f3.get(r), std::numeric_limits<int64>::min() );

    math::bconvex_rational::RationalFuncs<uint64> f4;
    f4.set(0, 1, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("0/1") );
    TS_ASSERT_EQUALS( f4.get(r), 0 );
    f4.set(1, 1, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("1/1") );
    TS_ASSERT_EQUALS( f4.get(r), 1 );
    f4.set(10, 10, r);
    os.str("");
    os << r;
    TS_ASSERT_EQUALS( os.str(), std::string("10/1") );
    TS_ASSERT_EQUALS( f4.get(r), 10 );
    f4.set(std::numeric_limits<uint32>::max(), std::numeric_limits<uint32>::max(), r);
    os.str("");
    os << r;
    os1.str("");
    os1 << std::numeric_limits<uint32>::max() << std::string("/1");
    TS_ASSERT_EQUALS( os.str(), os1.str() );
    TS_ASSERT_EQUALS( f4.get(r), std::numeric_limits<uint32>::max() );
    f4.set(std::numeric_limits<uint32>::min(), -std::numeric_limits<uint32>::min(), r);
    os.str("");
    os << r;
    os1.str("");
    os1 << std::numeric_limits<uint32>::min() << std::string("/1");
    TS_ASSERT_EQUALS( os.str(), os1.str() );
    TS_ASSERT_EQUALS( f4.get(r), std::numeric_limits<uint32>::min() );
    //f4.set(std::numeric_limits<uint64>::max(), std::numeric_limits<uint64>::max(), r);
    //os.str("");
    //os << r;
    //os1.str("");
    //os1 << std::numeric_limits<uint64>::max() << std::string("/1");
    //TS_ASSERT_EQUALS( os.str(), os1.str() );
    //TS_ASSERT_EQUALS( f4.get(r), std::numeric_limits<uint64>::max() );
    //f4.set(std::numeric_limits<uint64>::min(), -std::numeric_limits<uint64>::min(), r);
    //os.str("");
    //os << r;
    //os1.str("");
    //os1 << std::numeric_limits<uint64>::min() << std::string("/1");
    //TS_ASSERT_EQUALS( os.str(), os1.str() );
    //TS_ASSERT_EQUALS( f4.get(r), std::numeric_limits<uint64>::min() );

    math::bconvex_rational::RationalFuncs<float64> f5;
    float64 v5;
    f5.set(0, 1, r);
    TS_ASSERT_DELTA( f5.get(r), 0, 0.001 );
    f5.set(1, 1, r);
    TS_ASSERT_DELTA( f5.get(r), 1, 0.001 );
    f5.set(10, 10, r);
    TS_ASSERT_DELTA( f5.get(r), 10, 0.01 );
    f5.set(-11, 11, r);
    TS_ASSERT_DELTA( f5.get(r), -11, 0.01 );
    f5.set(9.187, 10, r);
    TS_ASSERT_DELTA( f5.get(r), 9.187, 0.01 );
    f5.set(-10.63412, 11, r);
    TS_ASSERT_DELTA( f5.get(r), -10.63412, 0.01 );
    f5.set(std::numeric_limits<int32>::max() / 2, std::numeric_limits<int32>::max() / 2, r);
    v5 = std::numeric_limits<int32>::max() / 4;
    TS_ASSERT( f5.get(r) > v5 );
    f5.set(std::numeric_limits<int32>::min() / 2, -(std::numeric_limits<int32>::min() / 2), r);
    v5 = std::numeric_limits<int32>::min() / 4;
    TS_ASSERT( f5.get(r) < v5 );
    //f5.set(std::numeric_limits<int64>::max() / 2, std::numeric_limits<int64>::max() / 2, r);
    //v5 = std::numeric_limits<int64>::max() / 4;
    //TS_ASSERT( f5.get(r) > v5 );
    //f5.set(std::numeric_limits<int64>::min() / 2, -(std::numeric_limits<int64>::min() / 2), r);
    //v5 = std::numeric_limits<int64>::min() / 4;
    //TS_ASSERT( f5.get(r) < v5 );
  }

  void test_bconvex()
  {
    std::vector<Vector3> v;
   
    BConvex c0(3); 
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
    TS_ASSERT( check_bconvex_print(os1.str(), std::string("-C >= -1, -B >= -1, -A >= -1, A >= 0, C >= 0, B >= 0")) );
    std::ostringstream os2;
    os2 << (c1 * 4);
    TS_ASSERT( check_bconvex_print(os2.str(), std::string("-C >= -4, -B >= -4, -A >= -4, A >= 0, C >= 0, B >= 0")) );
    std::ostringstream os3;
    os3 << (4 * c1);
    TS_ASSERT( check_bconvex_print(os3.str(), std::string("-C >= -4, -B >= -4, -A >= -4, A >= 0, C >= 0, B >= 0")) );
    std::ostringstream os4;
    os4 << ((c1 * 4) / 4);
    TS_ASSERT( check_bconvex_print(os4.str(), std::string("-C >= -1, -B >= -1, -A >= -1, A >= 0, C >= 0, B >= 0")) );
    std::ostringstream os5;
    os5 << (c1 + Vector3(1, -2, 3));
    TS_ASSERT( check_bconvex_print(os5.str(), std::string("-C >= -4, -B >= 1, -A >= -2, A >= 1, C >= 3, B >= -2")) );
    std::ostringstream os6;
    os6 << ((c1 + Vector3(1, -2, 3)) - Vector3(1, -2, 3));
    TS_ASSERT( check_bconvex_print(os6.str(), std::string("-C >= -1, -B >= -1, -A >= -1, A >= 0, C >= 0, B >= 0")) );

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
    
    BConvex c6(3);
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
  }

}; // class TestBConvex
