// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>
#include <vw/FileIO/EncodedImageResource.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ViewImageResource.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageIO.h>
#include <fstream>

#include <vw/FileIO/EncodedImageResourceJPEG.h>
using namespace vw;
using namespace std;

struct EncodedImageResourceTest : public ::testing::Test {
  void slurp(const string& filename, VarArray<uint8>& data) {
    ifstream f(filename.c_str(), ios::binary);
    ASSERT_TRUE(f.is_open());

    f.seekg(0, ios_base::end);
    ASSERT_FALSE(f.fail());

    size_t size = f.tellg();
    ASSERT_FALSE(f.fail());

    f.seekg(0, ios_base::beg);
    ASSERT_FALSE(f.fail());

    data.resize(size, false);
    f.read(reinterpret_cast<char*>(data.begin()), size);
    ASSERT_FALSE(f.fail());
  }
};

TEST_F(EncodedImageResourceTest, BasicRead) {
  const std::string fn("rgb2x2.jpg");

  VarArray<uint8> raw;
  slurp(fn, raw);

  SrcEncodedImageResourceJPEG r(raw.begin(), raw.size());

  ImageView<PixelRGB<uint8> > img1, img2;

  read_image(img1, fn);
  read_image(img2, r);

  ASSERT_EQ(img1.rows(), img2.rows());
  ASSERT_EQ(img1.cols(), img2.cols());
  ASSERT_EQ(img1.planes(), img2.planes());

  EXPECT_RANGE_EQ(img1.begin(), img1.end(), img2.begin(), img2.end());
}

TEST_F(EncodedImageResourceTest, BasicWrite) {
  const std::string src("rgb2x2.jpg");

  VarArray<uint8> data;
  ImageView<PixelRGB<uint8> > img1;

  DstEncodedImageResourceJPEG r(data);

  read_image(img1, src);
  write_image(r, img1);

  VarArray<uint8> raw;
  slurp(src, raw);

  ASSERT_EQ(raw.size(), data.size());
  EXPECT_RANGE_EQ(raw.begin(), raw.end(), data.begin(), data.end());
}
