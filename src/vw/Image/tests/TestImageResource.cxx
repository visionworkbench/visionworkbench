// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestImageResource.h
#include <gtest/gtest.h>

#include <vw/Image/ImageResource.h>
#include <vw/Image/ImageResourceImpl.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelTypes.h>

#include <test/Helpers.h>

#if defined(VW_HAVE_PKG_OPENCV) && (VW_HAVE_PKG_OPENCV==1)
# include <opencv/cxcore.h>
#endif

using namespace vw;

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

#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1

struct ImageResourceOpenCVTest : public ::testing::Test, private boost::noncopyable {
    typedef boost::shared_ptr<cv::Mat>             mat_t;
    typedef boost::shared_ptr<ImageResource>       wir_t;
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
    size_t cols = buf.format.cols, rows = buf.format.rows;
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
