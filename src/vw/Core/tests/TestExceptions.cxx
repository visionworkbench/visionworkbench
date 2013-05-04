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


#include <vw/Core/Exception.h>
#include <test/Helpers.h>

using namespace vw;

VW_DEFINE_EXCEPTION(Level1Err, Exception);
VW_DEFINE_EXCEPTION(Level2Err, Level1Err);

VW_DEFINE_EXCEPTION_EXT(Code0, Exception) {
  VW_EXCEPTION_API(Code0);
  virtual int32 code() const {return -1;}
};

VW_DEFINE_EXCEPTION(Code1, Code0);

VW_DEFINE_EXCEPTION_EXT(Code2, Code1) {
  VW_EXCEPTION_API(Code2);
  int32 code() const {return 4;}
};

TEST(Exceptions, HAS_EXCEPTIONS(Hierarchy)) {

  // Test exception names
  Level1Err l1;
  Level2Err l2;
  l1 << "Test message1.";
  l2 << "Test message2.";
  EXPECT_EQ(l1.name() , "Level1Err" );
  EXPECT_EQ(l2.name() , "Level2Err" );

  EXPECT_THROW(throw Level1Err(), Exception);
  EXPECT_THROW(throw Level1Err(), Level1Err);
  EXPECT_THROW(throw Level2Err(), Exception);
  EXPECT_THROW(throw Level2Err(), Level1Err);
  EXPECT_THROW(throw Level2Err(), Level2Err);
}

TEST(Exceptions, HAS_EXCEPTIONS(Ext)) {
  Code2 c;
  c << "rawr";
  EXPECT_EQ("Code2", c.name());
  EXPECT_THROW(throw Code2(), Code2);
  EXPECT_THROW(throw Code2(), Code1);
  EXPECT_THROW(throw Code2(), Code0);
  EXPECT_THROW(throw Code2(), Exception);

  EXPECT_THROW(vw_throw( Code2() << "Rawr" ), Code2);
  EXPECT_THROW(vw_throw( Code2() << "Rawr" ), Code1);
  EXPECT_THROW(vw_throw( Code2() << "Rawr" ), Code0);
  EXPECT_THROW(vw_throw( Code2() << "Rawr" ), Exception);

  try {
  } catch (const Code0& c) {
    EXPECT_EQ(2, c.code());
    EXPECT_EQ("Code2", c.name());
  }
}
