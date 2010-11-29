// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>
#include <vw/FileIO/MemoryImageResource.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ViewImageResource.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageMath.h>
#include <boost/filesystem/convenience.hpp>
#include <fstream>

namespace fs = boost::filesystem;
using namespace vw;
using namespace std;

struct MemoryImageResourceTest : public ::testing::TestWithParam<std::string> {
  void slurp(const string& filename, vector<uint8>& data) {
    ifstream f(filename.c_str(), ios::binary);
    ASSERT_TRUE(f.is_open());

    f.seekg(0, ios_base::end);
    ASSERT_FALSE(f.fail());

    size_t size = f.tellg();
    ASSERT_FALSE(f.fail());

    f.seekg(0, ios_base::beg);
    ASSERT_FALSE(f.fail());

    data.resize(size, false);
    f.read(reinterpret_cast<char*>(&*data.begin()), size);
    ASSERT_FALSE(f.fail());
  }
};

TEST_P(MemoryImageResourceTest, BasicRead) {
  vector<uint8> raw;
  slurp(GetParam(), raw);

  boost::scoped_ptr<SrcMemoryImageResource> r(SrcMemoryImageResource::open(fs::extension(GetParam()), &*raw.begin(), raw.size()));

  ImageView<PixelRGB<uint8> > img1, img2;
  ImageView<PixelRGB<float> > img3, img4;

  read_image(img1, GetParam());
  // simple read
  read_image(img2, *r);
  // read w/ convert (and also reread the same resource)
  read_image(img4, *r);

  EXPECT_SEQ_EQ(img1, img2);

  img3 = channel_cast_rescale<float>(img1);
  EXPECT_SEQ_NEAR(img3, img4, 1e-7);
}

TEST_P(MemoryImageResourceTest, BasicWriteRead) {
  typedef PixelRGB<uint8> Px;
  ImageView<Px> src(64,64);
  std::string type(fs::extension(GetParam()));

  Px gray(0x7f, 0x7f, 0x7f), red(0x7f, 0, 0), green(0, 0x7f, 0), blue(0, 0, 0x7f);

  vw::fill(crop(src,  0,  0, 32, 32), gray);
  vw::fill(crop(src, 32,  0, 32, 32), red);
  vw::fill(crop(src,  0, 32, 32, 32), green);
  vw::fill(crop(src, 32, 32, 32, 32), blue);


  vector<uint8> data;
  {
    // Set up dst to bytes go into data
    boost::scoped_ptr<DstMemoryImageResource> r(DstMemoryImageResource::create(type, &data, src.format()));
    write_image(*r, src);
    // data now contains the jpeg-encoded bytes
  }

  ImageView<Px> img1;
  {
    boost::scoped_ptr<SrcMemoryImageResource> r(SrcMemoryImageResource::open(type, &data[0], data.size()));
    read_image(img1, *r);
  }

  ASSERT_EQ(src.rows(),   img1.rows());
  ASSERT_EQ(src.cols(),   img1.cols());
  ASSERT_EQ(src.planes(), img1.planes());

  EXPECT_RANGE_NEAR(src.begin(), src.end(), img1.begin(), img1.end(), 2);
}

vector<string> test_paths() {
  vector<string> v;
#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
  v.push_back("rgb2x2.jpg");
#endif
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
  v.push_back("rgb2x2.png");
#endif
    return v;
}

INSTANTIATE_TEST_CASE_P(FileNames, MemoryImageResourceTest, ::testing::ValuesIn(test_paths()));

