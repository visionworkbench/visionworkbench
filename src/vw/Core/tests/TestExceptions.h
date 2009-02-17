// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
