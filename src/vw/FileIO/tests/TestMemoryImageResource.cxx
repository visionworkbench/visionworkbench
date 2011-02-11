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
#include <vw/Image/Filter.h>
#include <vw/Image/UtilityViews.h>
#include <boost/filesystem/convenience.hpp>
#include <boost/random/linear_congruential.hpp>
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
  typedef PixelRGBA<uint8> Px;
  ImageView<Px> src;

  {
    typedef PixelRGBA<float> Py;
    const size_t SIZE = 64;
    ImageView<Py> src_(SIZE,SIZE);
    for (size_t row = 0; row < SIZE; ++row) {
      for (size_t col = 0; col < SIZE; ++col) {
        src_(col, row) =
          Py(float(row)/SIZE, float(col)/SIZE, 1 - ((float(row) + col) / 2 / SIZE), 1);
      }
    }
    // jpeg is lossy, and has trouble with noise-free images. Add some noise and blur to help it out.
    boost::rand48 gen(uint64(test::get_random_seed()));
    src_ += gaussian_noise_view(gen, 0.008, 0.004, src_);
    src = gaussian_filter(pixel_cast<Px>(normalize(src_, 0, 255)), 2, 2, 4, 4);
    vw::fill(vw::select_channel(src, 3), 255);
  }

  std::string type(fs::extension(GetParam()));

  boost::scoped_ptr<SrcImageResource> src2;
  boost::scoped_ptr<DstMemoryImageResource> dst;

  ASSERT_NO_THROW(dst.reset(DstMemoryImageResource::create(type, src.format())));
  EXPECT_NO_THROW(write_image(*dst, src));
  ASSERT_NO_THROW(src2.reset(SrcMemoryImageResource::open(type, dst->data(), dst->size())));

  ImageView<Px> img1;
  read_image(img1, *src2);

  EXPECT_SEQ_NEAR(src, img1, 6);
}

vector<string> test_paths() {
  vector<string> v;
#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
  v.push_back("rgb2x2.jpg");
#endif
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
  v.push_back("rgb2x2.png");
#endif
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  v.push_back("rgb2x2.tif");
#endif
    return v;
}

INSTANTIATE_TEST_CASE_P(FileNames, MemoryImageResourceTest, ::testing::ValuesIn(test_paths()));

