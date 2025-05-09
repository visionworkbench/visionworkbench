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
#include <vw/Image/Algorithms.h>
#include <vw/Image/Filter.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageResource.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/UtilityViews.h>
#include <vw/FileIO/MemoryImageResource.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/config.h>

#include <fstream>

#include <boost/filesystem/path.hpp>
#include <boost/random/linear_congruential.hpp>


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

TEST_P(MemoryImageResourceTest, Zero) {
  // wrong data, zero len
  EXPECT_THROW(SrcMemoryImageResource::open(fs::path(GetParam()).extension().string(), (uint8*)42,   0),  ArgumentErr);
  // zero data, wrong len
  EXPECT_THROW(SrcMemoryImageResource::open(fs::path(GetParam()).extension().string(), NULL, 42),         ArgumentErr);
}

TEST_P(MemoryImageResourceTest, BasicRead) {
  vector<uint8> raw;
  slurp(GetParam(), raw);

  boost::scoped_ptr<SrcMemoryImageResource> r(SrcMemoryImageResource::open(fs::path(GetParam()).extension().string(), &*raw.begin(), raw.size()));

  ImageView<PixelRGBA<uint8> > disk, mem;
  ImageView<PixelRGBA<float> > diskf, memf;

  read_image(disk, GetParam());
  // simple read
  read_image(mem, *r);
  // read w/ convert (and also reread the same resource)
  read_image(memf, *r);

  ASSERT_SEQ_EQ(disk, mem);

  diskf = channel_cast_rescale<float>(disk);
  EXPECT_SEQ_NEAR(diskf, memf, 2e-3);
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
    boost::rand48 gen((uint64(t::get_random_seed())));
    src_ += gaussian_noise_view(gen, 0.008, 0.004, src_);
    src = gaussian_filter(pixel_cast<Px>(normalize(src_, 0, 255)), 2, 2, 4, 4);
    vw::fill(vw::select_channel(src, 3), 255);
  }

  std::string type(fs::path(GetParam()).extension().string());

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
  v.push_back("rgb4x4_halfalpha.png");
  v.push_back("rgb4x4_alpha.png");
#endif
#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  v.push_back("rgb2x2.tif");
  v.push_back("rgb4x4_halfalpha.tif");
  v.push_back("rgb4x4_alpha.tif");
  v.push_back("rgb4x4f_alpha.tif");
  v.push_back("rgb4x4f_band.tif");
#endif
  return v;
}

INSTANTIATE_TEST_CASE_P(FileNames, MemoryImageResourceTest, ::testing::ValuesIn(test_paths()));

#if defined(VW_HAVE_PKG_OPENEXR) && VW_HAVE_PKG_OPENEXR==1

// OpenEXR doesn't truly support uint8 only floats. So we'll write our own test.
template <typename PixelT>
class MemoryOpenEXR : public ::testing::Test {
protected:

  MemoryOpenEXR() {}

  virtual void SetUp() {
    ImageView<PixelT> src;

    {
      const size_t SIZE = 64;
      ImageView<PixelT> src_(SIZE,SIZE);
      for (size_t row = 0; row < SIZE; ++row) {
        for (size_t col = 0; col < SIZE; ++col) {
          for ( size_t ch = 0; ch < CompoundNumChannels<PixelT>::value; ch++ ) {
            switch (ch) {
            case 0:
              src_(col,row)[ch] = float(row)/SIZE; break;
            case 1:
              src_(col,row)[ch] = float(col)/SIZE; break;
            case 2:
              src_(col,row)[ch] = 1 - ((float(row) + col) / 2 / SIZE); break;
            default:
              src_(col,row)[ch] = float( ch + row ) / SIZE;
            }
          }
        }
      }

      boost::rand48 gen((uint64(t::get_random_seed())));
      src_ += gaussian_noise_view(gen, 0.008, 0.004, src_);
      src = gaussian_filter(src_, 2, 2, 4, 4);
      vw::fill(vw::select_channel(src,CompoundNumChannels<PixelT>::value - 1), 0.788);
    }

    std::string type("exr");

    boost::scoped_ptr<SrcImageResource> src2;
    boost::scoped_ptr<DstMemoryImageResource> dst;

    ASSERT_NO_THROW(dst.reset(DstMemoryImageResource::create(type, src.format())));
    EXPECT_NO_THROW(write_image(*dst, src));
    ASSERT_NO_THROW(src2.reset(SrcMemoryImageResource::open(type, dst->data(), dst->size())));

    ImageView<PixelT> img1;
    read_image(img1, *src2);

    EXPECT_SEQ_NEAR(src, img1, 1);
  }
};

// There's an internal switch in the code that treats RGB and Grays differently.
typedef MemoryOpenEXR<PixelRGB<float> > MemoryOpenEXRRGBF32;
TEST_F( MemoryOpenEXRRGBF32, RGB_F32 ) {}
typedef MemoryOpenEXR<PixelRGBA<float> > MemoryOpenEXRRGBAF32;
TEST_F( MemoryOpenEXRRGBAF32, RGBA_F32 ) {}
typedef MemoryOpenEXR<PixelGray<float> > MemoryOpenEXRGRAYF32;
TEST_F( MemoryOpenEXRGRAYF32, GRAY_F32 ) {}
typedef MemoryOpenEXR<PixelGrayA<float> > MemoryOpenEXRGRAYAF32;
TEST_F( MemoryOpenEXRGRAYAF32, GRAYA_F32 ) {}

#endif
