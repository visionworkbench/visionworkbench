// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <vw/Image/ImageResource.h>
#include <vw/Image/ImageResourceImpl.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelTypes.h>

#include <test/Helpers.h>

#if defined(VW_HAVE_PKG_OPENCV) && (VW_HAVE_PKG_OPENCV==1)
# include <opencv/cxcore.h>
#endif

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
namespace io = boost::iostreams;

#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

using namespace vw;
using namespace vw::test;

// This tests whether premultiplication preserves integer data.
TEST( ImageResource, PreMultiply ) {
  typedef PixelGrayA<uint16> Px;

  ImageFormat fmt;
  fmt.rows = fmt.cols = 2;
  fmt.planes = 1;
  fmt.pixel_format = PixelFormatID<Px>::value;
  fmt.channel_type = ChannelTypeID<PixelChannelType<Px>::type>::value;

  // Some arbitrary data (48195 triggered a rounding bug)
  // Small alpha AND value will basically always truncate.
  Px buf1_data[4] = {
    Px(10023,255), Px(48195,48195), Px(65535, 3), Px(1,65535)
  };
  Px buf2_data[4], buf3_data[4];

  ImageBuffer buf1(fmt, buf1_data, true);
  ImageBuffer buf2(fmt, buf2_data, false);
  ImageBuffer buf3(fmt, buf3_data, true);

  convert(buf2, buf1, true);
  convert(buf3, buf2, true);

  for (size_t i = 0; i < 4; ++i)
    EXPECT_PIXEL_EQ( buf3_data[i], buf1_data[i] );
}

class SrcNoopResource : public SrcImageResource {
  private:
    const ImageFormat& m_fmt;
    const uint8* m_data;
    bool m_copy;

  public:
    SrcNoopResource(const ImageFormat& fmt, const uint8* data, bool copy)
      : m_fmt(fmt), m_data(data), m_copy(copy) {}

    virtual ImageFormat format() const {return m_fmt;}
    virtual void read( ImageBuffer const& buf, BBox2i const& /*bbox*/ ) const { vw::convert(buf, ImageBuffer(m_fmt, const_cast<uint8*>(m_data))); }
    virtual bool has_block_read() const  {return false;}
    virtual bool has_nodata_read() const {return false;}
    virtual boost::shared_array<const uint8> native_ptr() const {
      return m_copy ? SrcImageResource::native_ptr() : boost::shared_array<const uint8>(m_data, NOP());
    }
};

TEST( ImageResource, NativePtr ) {
  ImageFormat fmt;
  fmt.cols = fmt.rows = 2;
  fmt.planes = 1;
  fmt.pixel_format = VW_PIXEL_RGBA;
  fmt.channel_type = VW_CHANNEL_UINT8;

  typedef PixelRGBA<uint8> Px;
  typedef boost::shared_array<const uint8> Data;
  const Px src_[4] = {Px(1,2,3,4), Px(5,6,7,8), Px(9,10,11,12), Px(13,14,15,16)};
  const uint8* src = reinterpret_cast<const uint8*>(&src_[0]);

  SrcNoopResource src1(fmt, src, true),
                  src2(fmt, src, false);

  ASSERT_EQ(2 * 2 * 1 * 4 * 1, src1.native_size());
  ASSERT_EQ(2 * 2 * 1 * 4 * 1, src2.native_size());

  Data d1, d2;
  d1 = src1.native_ptr();
  d2 = src2.native_ptr();

  ASSERT_NE(src, d1.get());
  ASSERT_EQ(src, d2.get());

  EXPECT_RANGE_EQ(src, src+4, &d1[0], &d1[4]);
  EXPECT_RANGE_EQ(src, src+4, &d2[0], &d2[4]);
}

struct TestStream : public ::testing::Test {
  protected:
    static const size_t WIDTH = 2;
    static const size_t HEIGHT = 3;
    static const size_t SIZE = WIDTH * HEIGHT;
  private:
    static const char canary[SIZE];
  protected:
    static const char src_data[SIZE];

    char src_buf[SIZE];
    char dst_buf[SIZE];
    char back_buf[SIZE];
    ImageFormat fmt;

    size_t size() const {return SIZE;}

    void SetUp() {
      fmt.cols = WIDTH;
      fmt.rows = HEIGHT;
      fmt.planes = 1;
      fmt.pixel_format = VW_PIXEL_GRAY;
      fmt.channel_type = VW_CHANNEL_UINT8;

      std::copy(src_data, src_data+SIZE, src_buf);
    }

    void CHECK(const char* e, const char* a, size_t size = 0) const {
      if (size == 0) size = SIZE;
      EXPECT_RANGE_EQ(e, e+SIZE, a, a+SIZE);
    }
    void CANARY(char* loc, size_t size = 0) const {
      assert(size <= SIZE);
      if (size == 0) size = SIZE;
      std::copy(canary+0, canary+size, loc);
    }
};
const char TestStream::src_data[TestStream::SIZE] = {'P', 'a', 'n', 't', 's', '!'};
const char TestStream::canary[TestStream::SIZE] = {43, 47, 53, 59, 61, 67};

TEST_F(TestStream, Read) {
  const BBox2i box(0,0,WIDTH,HEIGHT);
  io::stream<io::array_source> s1(src_data, size());
  SrcImageResourceStream r(&s1);

  ImageBuffer dst(fmt, dst_buf);

  CANARY(dst_buf);
  ASSERT_NO_THROW(r.read(dst, box));
  CHECK(src_data, dst_buf);

  // reset the read pointer for another read
  ASSERT_NO_THROW(r.reset());

  // repeat read to make sure it's supported (after we seek).
  CANARY(dst_buf);
  ASSERT_NO_THROW(r.read(dst, box));
  CHECK(src_data, dst_buf);
}

TEST_F(TestStream, Write) {
  const BBox2i box(0,0,WIDTH,HEIGHT);
  io::stream<io::array_sink> s1(dst_buf, size());

  DstImageResourceStream r(&s1);

  ImageBuffer src(fmt, src_buf);

  // install canary
  CANARY(dst_buf);
  ASSERT_NO_THROW(r.write(src, box));
  ASSERT_NO_THROW(r.flush());
  CHECK(src_data, dst_buf);

  // reset the write pointer for another write
  ASSERT_NO_THROW(r.reset());

  // repeat write to make sure it's supported (after we reset)
  CANARY(dst_buf);
  ASSERT_NO_THROW(r.reset());
  ASSERT_NO_THROW(r.write(src, box));
  ASSERT_NO_THROW(r.flush());
  CHECK(src_data, dst_buf);
}

TEST_F(TestStream, ReadWrite) {
  const BBox2i box(0,0,WIDTH,HEIGHT);

  io::stream<io::array> ss(back_buf, back_buf+SIZE);
  ImageResourceStream r(&ss);
  EXPECT_THROW(r.cols(), LogicErr);

  ImageBuffer src(fmt, src_buf);
  ImageBuffer dst(fmt, dst_buf);

  CANARY(dst_buf);
  CANARY(back_buf);

  // Do a normal write, a reset, and then a normal read
  ASSERT_NO_THROW(r.write(src, box));
  ASSERT_NO_THROW(r.flush());
  CHECK(src_data, back_buf);

  // reset the pointer for the reread
  ASSERT_NO_THROW(r.reset());

  ASSERT_NO_THROW(r.read(dst, box));
  CHECK(src_data, dst_buf);
}

TEST_F(TestStream, WriteResize) {
  const BBox2i box(0,0,WIDTH,HEIGHT);

  UnlinkName name("resizewrite");
  std::ofstream f(name.c_str(), std::ios::binary | std::ios::trunc);
  ASSERT_TRUE(f.is_open());

  DstImageResourceStream r(&f);

  ImageBuffer src(fmt, src_buf);

  EXPECT_FALSE(f.bad());
  //ASSERT_NO_THROW(r.write(src, box));
  r.write(src, box);
  ASSERT_NO_THROW(r.flush());
  ASSERT_EQ(size(), fs::file_size(name));

  // reset the pointer for the rewrite
  ASSERT_NO_THROW(r.reset());

  // repeat write to make sure it doesn't resize even bigger (it should overwrite old)
  ASSERT_NO_THROW(r.write(src, box));
  ASSERT_NO_THROW(r.flush());
  ASSERT_EQ(size(), fs::file_size(name));

  // this should now append to the file
  ASSERT_NO_THROW(r.write(src, box));
  ASSERT_NO_THROW(r.flush());
  ASSERT_EQ(2*size(), fs::file_size(name));
}

#if 0
TEST_F(TestStream, ReadWriteConvert) {

#if 0
  ImageBuffer src2(fmt, src_buf);
  // Now that we've done a write, we can do a convert read
  src2.format.pixel_format = VW_PIXEL_SCALAR;
  // Set up canary
  std::copy(idata+0, idata+SIZE, dst_data);

  ASSERT_NO_THROW(r.read(dst, box));
  EXPECT_RANGE_EQ(rstr.begin(), rstr.end(), ss.str().begin(), ss.str().end());
  EXPECT_RANGE_EQ(rdata+0, rdata+SIZE, dst_data+0, dst_data+SIZE);
#endif
}
#endif


TEST_F(TestStream, StreamError) {
  const BBox2i box(0,0,WIDTH,HEIGHT);
  const BBox2i big_box(0,0,WIDTH+1,HEIGHT+1);
  const size_t BIG_SIZE = (WIDTH+1) * (HEIGHT+1);
  const char src_data_big[BIG_SIZE] = {'P', 'a', 'n', 't', 's', '!', 'H', 'e', 'l', 'l', 'o', '!'};
  char dst_buf_big[BIG_SIZE], src_buf_big[BIG_SIZE];

  io::stream<io::array> ss(back_buf+0, back_buf+SIZE);
  ImageResourceStream r(&ss);

  ImageFormat big_fmt = fmt;
  big_fmt.cols = WIDTH+1;
  big_fmt.rows = HEIGHT+1;

  // Set up canary in dst
  CANARY(back_buf);

  // Set up source data (since imagebuffer needs non-const void*)
  std::copy(src_data_big, src_data_big+BIG_SIZE, src_buf_big);

  ImageBuffer src(fmt, &src_buf);
  ImageBuffer big_src(big_fmt, &src_buf_big);

  // Try a write without enough room to hold it
  EXPECT_FALSE(ss.bad());
  EXPECT_THROW(r.write(big_src, big_box), IOErr);
  // Make sure that the error is cleared, we already "handled" it by throwing
  EXPECT_FALSE(ss.bad());

  // Now retry the write with a proper size
  ASSERT_NO_THROW(r.reset());
  ASSERT_NO_THROW(r.write(src, box));
  EXPECT_FALSE(ss.bad());
  CHECK(src_data, back_buf);

  // Now try a read without enough data
  CANARY(dst_buf);
  std::copy(src_data+0, src_data+SIZE, back_buf);

  ImageBuffer dst(fmt, &dst_buf);
  ImageBuffer big_dst(big_fmt, &dst_buf_big);

  EXPECT_FALSE(ss.fail());
  EXPECT_THROW(r.read(big_dst, big_box), IOErr);
  // Make sure that the error is cleared, we already "handled" it by throwing
  EXPECT_FALSE(ss.fail());

  ASSERT_NO_THROW(r.reset());

  // Now retry the read with a proper size
  ASSERT_NO_THROW(r.read(dst, box));
  EXPECT_FALSE(ss.fail());
  CHECK(src_data, dst_buf);
}

#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1

struct ImageResourceOpenCVTest : public ::testing::Test, private boost::noncopyable {
    typedef boost::shared_ptr<cv::Mat>       mat_t;
    typedef boost::shared_ptr<ImageResource> wir_t;
    typedef boost::shared_ptr<const ImageResource> rir_t;
    typedef uint8  in_px_t;
    typedef double out_px_t;
  private:
    static const size_t width   = 2;
    static const size_t height  = 4;
  public:
    static const size_t len = width * height;
  protected:
    static in_px_t  input_data[len];
    static out_px_t actual_data[len];
    static const out_px_t e_read_data[len];
    static const out_px_t e_write_data[len];

  void SetUp() {
    // these values chosen to be the same as the values in expected, below.
    // they act as canaries in case of input overrun.
    static const out_px_t pre[len] = {97,89,83,79,73,71,67,61};
    std::copy(pre, pre+len, actual_data);
  }

  template <typename T>
  void make_matrix(T data[], mat_t& matrix) const {
    matrix.reset(new cv::Mat_<T>(height, width, data, cv::Mat::AUTO_STEP));
    ASSERT_TRUE(matrix->isContinuous()) << "Assumption failure: Matrix should be continuous!";
  }

  void make_discontinuous(mat_t& matrix) const {
    matrix->flags &= ~cv::Mat::CONTINUOUS_FLAG;
    ASSERT_FALSE(matrix->isContinuous()) << "Assumption failure: Matrix should be discontinuous!";
  }

  BBox2i box(const ImageBuffer& buf) const {
    int cols = buf.format.cols, rows = buf.format.rows;
    return BBox2i(0,rows/2,cols,rows);
  }

  enum Mode {
    READ,
    WRITE
  };

  ImageBuffer buffer(Mode m) const {
    ImageFormat fmt;
    fmt.cols   = width;
    fmt.rows   = height/2;
    fmt.planes = 1;
    fmt.pixel_format      = m == READ ? PixelFormatEnum(PixelFormatID<in_px_t>::value)
                                      : PixelFormatEnum(PixelFormatID<out_px_t>::value);
    fmt.channel_type      = m == READ ? ChannelTypeEnum(ChannelTypeID<PixelChannelType<in_px_t>::type>::value)
                                      : ChannelTypeEnum(ChannelTypeID<PixelChannelType<out_px_t>::type>::value);
    return ImageBuffer(fmt, m == READ ? reinterpret_cast<void*>(input_data)
                                      : reinterpret_cast<void*>(actual_data));
  }
};

ImageResourceOpenCVTest::in_px_t ImageResourceOpenCVTest::input_data[ImageResourceOpenCVTest::len] = {
   2,  3,
  11, 13,
  23, 29,
  41, 43,
};

const ImageResourceOpenCVTest::out_px_t ImageResourceOpenCVTest::e_read_data[ImageResourceOpenCVTest::len]  = {
  11, 13,
  23, 29,
  73,71,
  67,61
};

const ImageResourceOpenCVTest::out_px_t ImageResourceOpenCVTest::e_write_data[ImageResourceOpenCVTest::len]  = {
  97,89,
   2, 3,
  11,13,
  67,61
};
ImageResourceOpenCVTest::out_px_t ImageResourceOpenCVTest::actual_data[ImageResourceOpenCVTest::len];


TEST_F(ImageResourceOpenCVTest, OpenCv_Read_Cont) {
  mat_t m;
  make_matrix<in_px_t>(input_data, m);
  rir_t r(new ImageResourceOpenCV(m));
  ImageBuffer buf = buffer(WRITE);

  ASSERT_NO_THROW(r->read(buf, box(buf)));
  EXPECT_RANGE_EQ(e_read_data+0, e_read_data+len, actual_data+0, actual_data+len);
}

TEST_F(ImageResourceOpenCVTest, OpenCv_Read_Discont) {
  mat_t m;
  make_matrix<in_px_t>(input_data, m);
  make_discontinuous(m);

  rir_t r(new ImageResourceOpenCV(m));
  ImageBuffer buf = buffer(WRITE);

  ASSERT_NO_THROW(r->read(buf, box(buf)));
  EXPECT_RANGE_EQ(e_read_data+0, e_read_data+len, actual_data+0, actual_data+len);
}

TEST_F(ImageResourceOpenCVTest, OpenCv_Write_Cont) {
  mat_t m;
  make_matrix<out_px_t>(actual_data, m);
  wir_t r(new ImageResourceOpenCV(m));
  ImageBuffer buf = buffer(READ);

  ASSERT_NO_THROW(r->write(buf, box(buf)));
  EXPECT_RANGE_EQ(e_write_data+0, e_write_data+len, actual_data+0, actual_data+len);
}

TEST_F(ImageResourceOpenCVTest, OpenCv_Write_Discont) {
  mat_t m;
  make_matrix<out_px_t>(actual_data, m);
  make_discontinuous(m);

  wir_t r(new ImageResourceOpenCV(m));
  ImageBuffer buf = buffer(READ);

  ASSERT_NO_THROW(r->write(buf, box(buf)));
  EXPECT_RANGE_EQ(e_write_data+0, e_write_data+len, actual_data+0, actual_data+len);
}

#endif
