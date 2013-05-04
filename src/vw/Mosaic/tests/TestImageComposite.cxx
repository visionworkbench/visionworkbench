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


#include <gtest/gtest_VW.h>
#include <vw/Mosaic/ImageComposite.h>

using namespace std;
using namespace vw;
using namespace vw::mosaic;

ImageView<uint32> make(uint32 x) {
  ImageView<uint32> img(8,8);
  fill(img, x);
  return img;
}

TEST(TestImageComposite, Basic) {
  ImageComposite<uint32> c;
  c.set_draft_mode(true);

  c.insert(make(2), 3, 0);
  EXPECT_EQ(8, c.rows());
  EXPECT_EQ(8, c.cols());
  EXPECT_EQ(1, c.planes());

  c.insert(make(3), 0, 0);
  EXPECT_EQ(8, c.rows());
  EXPECT_EQ(11, c.cols());
  EXPECT_EQ(1, c.planes());

  for (int32 row = 0; row < 8; ++row) {
    for (int32 col = 0; col < 8; ++col)
      EXPECT_EQ(3, c(col, row)) << "at (" << col << "," << row << ")";
    for (int32 col = 8; col < 11; ++col)
      EXPECT_EQ(2, c(col, row)) << "at (" << col << "," << row << ")";
  }
}
