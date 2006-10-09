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

// Math/TestFunctions.h
#include <cxxtest/TestSuite.h>

#include <vw/Math/Functions.h>

using namespace vw;

class TestFunctions : public CxxTest::TestSuite
{
public:

  void test_erf() {
    TS_ASSERT_DELTA( vw::math::erf(-3.0), -0.9999779095030014, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erf(-1e0), -0.8427007929497149, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erf(-1e-1), -0.1124629160182849, 1e-13 );
    TS_ASSERT_DELTA( vw::math::erf(-1e-2), -0.01128341555584962, 1e-14 );
    TS_ASSERT_DELTA( vw::math::erf(-1e-4), -0.0001128379163334249, 1e-16 );
    TS_ASSERT_EQUALS( vw::math::erf(0), 0 );
    TS_ASSERT_DELTA( vw::math::erf(1e-4), 0.0001128379163334249, 1e-16 );
    TS_ASSERT_DELTA( vw::math::erf(1e-2), 0.01128341555584962, 1e-14 );
    TS_ASSERT_DELTA( vw::math::erf(1e-1), 0.1124629160182849, 1e-13 );
    TS_ASSERT_DELTA( vw::math::erf(1e0), 0.8427007929497149, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erf(3.0), 0.9999779095030014, 1e-12 );
  }

  void test_erfc() {
    TS_ASSERT_DELTA( vw::math::erfc(-3.0), 1.999977909503001, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erfc(-1e0), 1.842700792949715, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erfc(-1e-1), 1.112462916018285, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erfc(-1e-2), 1.011283415555850, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erfc(-1e-4), 1.000112837916333, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erfc(0), 1, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erfc(1e-4), 0.9998871620836666, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erfc(1e-2), 0.9887165844441504, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erfc(1e-1), 0.8875370839817151, 1e-12 );
    TS_ASSERT_DELTA( vw::math::erfc(1e0), 0.1572992070502851, 1e-13 );
    TS_ASSERT_DELTA( vw::math::erfc(3.0), 0.00002209049699858544, 1e-17 );
  }

}; // class TestFunctions
