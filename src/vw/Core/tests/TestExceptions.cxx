// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <vw/Core/Exception.h>
#include <test/Helpers.h>

using namespace vw;

VW_DEFINE_EXCEPTION(Level1Err, vw::Exception);
VW_DEFINE_EXCEPTION(Level2Err, Level1Err);

TEST(Exceptions, HAS_EXCEPTIONS(Hierarchy)) {

  // Test exception names
  Level1Err l1;
  Level2Err l2;
  l1 << "Test message1.";
  l2 << "Test message2.";
  EXPECT_EQ(l1.name() , "Level1Err" );
  EXPECT_EQ(l2.name() , "Level2Err" );

  EXPECT_THROW(throw Level1Err(), vw::Exception);
  EXPECT_THROW(throw Level1Err(), Level1Err);
  EXPECT_THROW(throw Level2Err(), vw::Exception);
  EXPECT_THROW(throw Level2Err(), Level1Err);
  EXPECT_THROW(throw Level2Err(), Level2Err);
}
