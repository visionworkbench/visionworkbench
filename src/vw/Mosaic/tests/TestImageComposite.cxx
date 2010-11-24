#include <gtest/gtest.h>
#include <vw/Mosaic/ImageComposite.h>

using namespace std;
using namespace vw;
using namespace vw::mosaic;

ImageView<uint32> make(uint32 x) {
  ImageView<uint32> img(1,1);
  img(0,0) = x;
  return img;
}

TEST(TestImageComposite, Basic) {
  ImageComposite<uint32> c;
  c.set_draft_mode(true);

  c.insert(make(2), 3, 0);
  EXPECT_EQ(1, c.rows());
  EXPECT_EQ(1, c.cols());
  EXPECT_EQ(1, c.planes());

  c.insert(make(3), 0, 0);
  EXPECT_EQ(1, c.rows());
  EXPECT_EQ(4, c.cols());
  EXPECT_EQ(1, c.planes());

}
