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

// TestFunctors.h
#include <cxxtest/TestSuite.h>

#include <vw/Core/Exception.h>

using namespace std;
using namespace vw;

VW_DEFINE_EXCEPTION(Level1Err, vw::Exception);
VW_DEFINE_EXCEPTION(Level2Err, Level1Err);

class TestExceptions : public CxxTest::TestSuite
{
public:
  void test_exception_hierarchy()
  {
#if defined(VW_ENABLE_EXCEPTIONS) && (VW_ENABLE_EXCEPTIONS==1)
    TS_ASSERT_THROWS(throw Level1Err(), vw::Exception);
    TS_ASSERT_THROWS(throw Level1Err(), Level1Err);
    TS_ASSERT_THROWS(throw Level2Err(), vw::Exception);
    TS_ASSERT_THROWS(throw Level2Err(), Level1Err);
    TS_ASSERT_THROWS(throw Level2Err(), Level2Err);
#endif
  }

};
