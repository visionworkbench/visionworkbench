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


#include <test/Helpers.h>
#include <vw/Image/BlockRasterize.h>

using namespace vw;
using namespace std;

TEST(BlockRasterize, Basic) {
  typedef ImageView<uint32> Image;
  typedef BlockRasterizeView<Image> Block;
  Vector2i block(1,1);

  Image img1(2,2), img2;
  img1(0,0) = 1; img1(0,1) = 2; img1(1,0) = 3; img1(1,1) = 4;

  Block b1(img1, block, 4);
  Block b2 = block_rasterize(img1, block, 4);
  Block b3 = block_cache(img1, block, 4);
  Block b4 = block_cache(img1, block, 4, vw_system_cache());

  img2 = b1;
  EXPECT_RANGE_EQ(img1.begin(), img1.end(), img2.begin(), img2.end());
  img2 = b2;
  EXPECT_RANGE_EQ(img1.begin(), img1.end(), img2.begin(), img2.end());
  img2 = b3;
  EXPECT_RANGE_EQ(img1.begin(), img1.end(), img2.begin(), img2.end());
  img2 = b4;
  EXPECT_RANGE_EQ(img1.begin(), img1.end(), img2.begin(), img2.end());
}
