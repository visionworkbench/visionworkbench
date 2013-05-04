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
#include <test/Helpers.h>
#include <vw/Plate/TileManipulation.h>
#include <boost/foreach.hpp>

using namespace std;
using namespace vw;
using namespace vw::platefile;
using namespace vw::test;

class BlobManagerTest : public ::testing::Test {
  protected:

#define bbox_does_tile_area(bbox, width, height, total_tiles) do {\
  SCOPED_TRACE("");\
  bbox_does_tile_area_(bbox, width, height, total_tiles);\
} while(0)

  void bbox_does_tile_area_(const BBox2i& orig_bbox, unsigned width, unsigned height, unsigned total_tiles) {
    list<BBox2i> tiles = bbox_tiles(orig_bbox, width, height);
    EXPECT_EQ(total_tiles, tiles.size());

    BBox2i total_bbox;

    // It should completely tile the space. No tile should overlap.
    BOOST_FOREACH( const BBox2i& bbox, tiles ) {
      SCOPED_TRACE(Message() << "Outer box[" << bbox << "]");

      // Tiles are allowed to be smaller than the requested size. Just not larger.
      EXPECT_LE(bbox.width(),  width);
      EXPECT_LE(bbox.height(), height);

      // One of the boxes should be the same.
      bool intersected = false;

      // Make sure it doesn't intersect any other box. (Since we're comparing a
      // list to itself, that means there should be just one intersection)
      BOOST_FOREACH( const BBox2i& new_box, tiles ) {
        SCOPED_TRACE(Message() << "Inner box[" << new_box << "]");
        if (bbox.intersects(new_box)) {
          EXPECT_FALSE(intersected);
          intersected = true;
        }
      }

      total_bbox.grow(bbox);
    }

    // Make sure that the union of all tiles is the original bbox
    EXPECT_EQ(orig_bbox, total_bbox);
  }
};

TEST_F(BlobManagerTest, Tiles) {
  bbox_does_tile_area(BBox2i(0,0,8,12), 4, 4, 6);
}

TEST_F(BlobManagerTest, Weird) {
  bbox_does_tile_area(BBox2i(7,17,31,13), 7, 5, 15);
}
