#include <gtest/gtest.h>
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
