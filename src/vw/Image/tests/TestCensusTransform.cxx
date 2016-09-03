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


// TestCensusTransform.h
#include <gtest/gtest_VW.h>
#include <vw/Image/CensusTransform.h>

using namespace vw;

TEST( CensusTransform, PointTests ) {
  ImageView<double> src(8,8); 
  src(0,0)=1; src(1,0)=2; src(2,0)=7; src(3,0)=2; src(4,0)=2; src(5,0)=8; src(6,0)=5; src(7,0)=2;
  src(0,1)=1; src(1,1)=4; src(2,1)=2; src(3,1)=9; src(4,1)=8; src(5,1)=8; src(6,1)=2; src(7,1)=6;
  src(0,2)=5; src(1,2)=2; src(2,2)=7; src(3,2)=2; src(4,2)=2; src(5,2)=2; src(6,2)=4; src(7,2)=6;
  src(0,3)=1; src(1,3)=2; src(2,3)=2; src(3,3)=2; src(4,3)=1; src(5,3)=4; src(6,3)=5; src(7,3)=2;
  src(0,4)=6; src(1,4)=6; src(2,4)=3; src(3,4)=7; src(4,4)=2; src(5,4)=2; src(6,4)=5; src(7,4)=5;
  src(0,5)=1; src(1,5)=2; src(2,5)=9; src(3,5)=2; src(4,5)=2; src(5,5)=2; src(6,5)=2; src(7,5)=2;
  src(0,6)=7; src(1,6)=9; src(2,6)=2; src(3,6)=8; src(4,6)=5; src(5,6)=2; src(6,6)=3; src(7,6)=2;
  src(0,7)=1; src(1,7)=2; src(2,7)=2; src(3,7)=2; src(4,7)=2; src(5,7)=2; src(6,7)=2; src(7,7)=1;

  
  EXPECT_EQ(get_census_value_3x3(src, 2,2), 0x20); // 00100000
  EXPECT_EQ(get_census_value_3x3(src, 4,5), 0x86); // 10000110
  EXPECT_EQ(get_census_value_3x3(src, 6,1), 0xDB); // 11011011
  
  EXPECT_EQ(get_census_value_5x5(src, 4,4), 0x0088F60D); // 00000000 100010001111011000001101
  EXPECT_EQ(get_census_value_5x5(src, 2,3), 0x005D03C4); // 00000000 010111010000001111000100
  
  EXPECT_EQ(get_census_value_7x7(src, 3,4), 0x00001C0000041400); // 0000000000000000 000111000000000000000000000001000001010000000000

}

TEST( HammingDist, Tests) {

  EXPECT_EQ(hamming_distance(uint8(0x01), uint8(0x00)), 1);
  EXPECT_EQ(hamming_distance(uint8(0xF0), uint8(0x00)), 4);
  EXPECT_EQ(hamming_distance(uint8(0xF0), uint8(0xB1)), 2);
  
  EXPECT_EQ(hamming_distance(uint32(0x00B06FFF), uint32(0x00B06F11)), 6);
  EXPECT_EQ(hamming_distance(uint32(0x0033C8BA), uint32(0x0023C0BA)), 2);
  
  EXPECT_EQ(hamming_distance(uint64(0x00002820A0F038), uint64(0x00002820000030)), 7);
}



