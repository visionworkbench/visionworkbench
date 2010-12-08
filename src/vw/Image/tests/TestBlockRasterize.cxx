// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
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
