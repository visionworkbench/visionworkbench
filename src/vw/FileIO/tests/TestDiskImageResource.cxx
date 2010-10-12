// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestDiskImageResource.h
#include <gtest/gtest.h>
#include <vw/FileIO.h>
#include <vw/FileIO/DiskImageResource_internal.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ImageView.h>
#include <vw/config.h>
#include <test/Helpers.h>

using namespace vw;
using namespace vw::internal;
using namespace vw::test;

using std::string;
using std::set;

template <class PixelT>
static void test_extension(string const& fn_base)
{
  ImageView<PixelT> img1(4,4), img2;
  UnlinkName fn(fn_base);

  try {
    write_image(fn, img1);
    read_image(img2, fn);
  } catch (vw::NoImplErr &e) {
    // this doesn't represent a test failure, just a lack of test coverage
    return;
  }

  EXPECT_EQ( img1.cols(),         img2.cols() ) << fn;
  EXPECT_EQ( img1.rows(),         img2.rows() ) << fn;
  EXPECT_EQ( img1.planes(),       img2.planes() ) << fn;
  EXPECT_EQ( img1.channels(),     img2.channels() ) << fn;
  EXPECT_EQ( img1.channel_type(), img2.channel_type() ) << fn;
}

TEST( DiskImageResource, WriteReadView ) {
  set<string> exclude;
  const char *ex_list[] = {"img", "lbl", "pds", "cub"}; // skip the ro PDS formats
  exclude.insert(ex_list, ex_list+4);

  foreach_ext("rwtest",
              test_extension<PixelRGB<float> >,   exclude);
  foreach_ext("rwtest",
              test_extension<PixelRGB<uint8> >,   exclude);
  foreach_ext("rwtest",
              test_extension<PixelRGBA<uint8> >,  exclude);
  foreach_ext("rwtest",
              test_extension<PixelGray<uint8> >,  exclude);
  foreach_ext("rwtest",
              test_extension<PixelGrayA<float> >, exclude);

  // there's no sane way to represent scalar images in ppm
  exclude.insert("ppm");
  foreach_ext("rwtest", test_extension<uint8>, exclude);
  foreach_ext("rwtest", test_extension<float>, exclude);
}

template <typename PixelT, int extension_i>
class ReadImage : public ::testing::Test {
protected:
  ReadImage() {}

  virtual void SetUp() {
    ImageView<PixelT> image;
    std::string filename;
    switch (extension_i) {
    case 0: // PNG
      filename = TEST_SRCDIR"/rgb2x2.png"; break;
    case 1: // TIF
      filename = TEST_SRCDIR"/rgb2x2.tif"; break;
    case 2: // JPG
    default:
      filename = TEST_SRCDIR"/rgb2x2.jpg"; break;
    }
    ASSERT_NO_THROW( read_image( image, filename ) );
    EXPECT_EQ( image.cols(), 2 );
    EXPECT_EQ( image.rows(), 2 );
    EXPECT_EQ( image.planes(), 1 );
    typedef typename CompoundChannelType<PixelT>::type ChannelT;
    ChannelT tol;
    if ( boost::is_floating_point<ChannelT>::value )
      tol = boost::numeric_cast<ChannelT>(extension_i == 2 ? 1e-2 : 1e-5);
    else
      tol = boost::numeric_cast<ChannelT>(extension_i == 2 ? 1 : 0);
    ChannelCastRescaleFunctor<ChannelT> conv;
    EXPECT_NEAR( image(0,0).r(), conv(uint8(128)), tol );
    EXPECT_NEAR( image(0,0).g(), conv(uint8(128)), tol );
    EXPECT_NEAR( image(0,0).b(), conv(uint8(128)), tol );
    EXPECT_NEAR( image(1,0).r(), conv(uint8(85)), tol );
    EXPECT_NEAR( image(1,0).g(), conv(uint8(0)), tol );
    EXPECT_NEAR( image(1,0).b(), conv(uint8(0)), tol );
    EXPECT_NEAR( image(0,1).r(), conv(uint8(0)), tol );
    EXPECT_NEAR( image(0,1).g(), conv(uint8(170)), tol );
    EXPECT_NEAR( image(0,1).b(), conv(uint8(0)), tol );
    EXPECT_NEAR( image(1,1).r(), conv(uint8(0)), tol );
    EXPECT_NEAR( image(1,1).g(), conv(uint8(0)), tol );
    EXPECT_NEAR( image(1,1).b(), conv(uint8(255)), tol );
  }
};

template <typename PixelT, int extension_i>
class WriteReadImage : public ::testing::Test {
protected:
  WriteReadImage() {}

  virtual void SetUp() {
    ImageView<PixelT> image(2,2);
    image(0,0) = pixel_cast_rescale<PixelT>(PixelRGB<uint8>(120,240,180));
    image(1,0) = pixel_cast_rescale<PixelT>(PixelRGB<uint8>(36,89,79));
    image(0,1) = pixel_cast_rescale<PixelT>(PixelRGB<uint8>(190,25,34));
    image(1,1) = pixel_cast_rescale<PixelT>(PixelRGB<uint8>(23,13,189));
    std::string base;
    switch ( extension_i ) {
    case 0: // PNG
      base = "tmp.png"; break;
    case 1: // TIF
      base = "tmp.tif"; break;
    case 2: // JPG
    default:
      DiskImageResourceJPEG::set_default_quality(1.0);
      DiskImageResourceJPEG::set_default_subsample_factor(1);
      base = "tmp.jpg"; break;
    }
    UnlinkName filename(base);
    ASSERT_NO_THROW( write_image( filename, image ) );
    ImageView<PixelT> image2;
    ASSERT_NO_THROW( read_image( image2, filename ) );
    EXPECT_EQ( image2.cols(), image.cols() );
    EXPECT_EQ( image2.rows(), image.rows() );
    EXPECT_EQ( image2.planes(), image.planes() );
    typedef typename CompoundChannelType<PixelT>::type ChannelT;
    ChannelT tol;
    if ( boost::is_floating_point<ChannelT>::value )
      tol = boost::numeric_cast<ChannelT>(extension_i == 2 ? 4e-2 : 1e-5);
    else
      tol = boost::numeric_cast<ChannelT>(extension_i == 2 ? 10 : 0);
    EXPECT_NEAR( image2(0,0).r(), image(0,0).r(), tol );
    EXPECT_NEAR( image2(0,0).g(), image(0,0).g(), tol );
    EXPECT_NEAR( image2(0,0).b(), image(0,0).b(), tol );
    EXPECT_NEAR( image2(1,0).r(), image(1,0).r(), tol );
    EXPECT_NEAR( image2(1,0).g(), image(1,0).g(), tol );
    EXPECT_NEAR( image2(1,0).b(), image(1,0).b(), tol );
    EXPECT_NEAR( image2(0,1).r(), image(0,1).r(), tol );
    EXPECT_NEAR( image2(0,1).g(), image(0,1).g(), tol );
    EXPECT_NEAR( image2(0,1).b(), image(0,1).b(), tol );
    EXPECT_NEAR( image2(1,1).r(), image(1,1).r(), tol );
    EXPECT_NEAR( image2(1,1).g(), image(1,1).g(), tol );
    EXPECT_NEAR( image2(1,1).b(), image(1,1).b(), tol );
  }
};

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
typedef ReadImage<PixelRGB<uint8>, 0> ReadImageRGBU8PNG;
TEST_F( ReadImageRGBU8PNG, RGB_U8_PNG ) {}

typedef ReadImage<PixelRGB<float>, 0> ReadImageRGBF32PNG;
TEST_F( ReadImageRGBF32PNG, RGB_F32_PNG ) {}

typedef WriteReadImage<PixelRGB<uint8>, 0> WriteReadImageRGBU8PNG;
TEST_F( WriteReadImageRGBU8PNG, RGB_U8_PNG ) {}

typedef WriteReadImage<PixelRGB<float>, 0> WriteReadImageRGBF32PNG;
TEST_F( WriteReadImageRGBF32PNG, RGB_F32_PNG ) {}
#endif

#if (defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1) || (defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1)
typedef ReadImage<PixelRGB<uint8>, 1> ReadImageRGBU8TIF;
TEST_F( ReadImageRGBU8TIF, RGB_U8_TIF ) {}

typedef ReadImage<PixelRGB<float>, 1> ReadImageRGBF32TIF;
TEST_F( ReadImageRGBF32TIF, RGB_F32_TIF ) {}

typedef WriteReadImage<PixelRGB<uint8>, 1> WriteReadImageRGBU8TIF;
TEST_F( WriteReadImageRGBU8TIF, RGB_U8_TIF ) {}

typedef WriteReadImage<PixelRGB<float>, 1> WriteReadImageRGBF32TIF;
TEST_F( WriteReadImageRGBF32TIF, RGB_F32_TIF ) {}
#endif

#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
typedef ReadImage<PixelRGB<uint8>, 2> ReadImageRGBU8JPG;
TEST_F( ReadImageRGBU8JPG, RGB_U8_JPG ) {}

typedef ReadImage<PixelRGB<float>, 2> ReadImageRGBF32JPG;
TEST_F( ReadImageRGBF32JPG, RGB_F32_JPG ) {}

// These test just don't work with JPG since these images are smaller
// than the tolerance of libjpeg.
typedef WriteReadImage<PixelRGB<uint8>, 2> WriteReadImageRGBU8JPG;
TEST_F( WriteReadImageRGBU8JPG, DISABLED_RGB_U8_JPG ) {}

typedef WriteReadImage<PixelRGB<float>, 2> WriteReadImageRGBF32JPG;
TEST_F( WriteReadImageRGBF32JPG, DISABLED_RGB_F32_JPG ) {}
#endif

TEST( DiskImageResource, TestPBM ) {
  const char
    //pr1[] = "P1 1 2 1 0",
    pr2[] = "P2 1 2 255 12 36",
    pr3[] = "P3 1 2 255 42 43 44 89 88 87",
    //pr4[] = "P4 1 2 \1\0",
    pr5[] = "P5 1 2 255 \xC\x24",
    pr6[] = "P6 1 2 255 \x2A\x2B\x2C\x59\x58\x57";

#define WF(x, ext) \
  UnlinkName fn ## x ("/test_p" #x "." ext); \
  do {                                                 \
    std::fstream f ## x(fn ## x.c_str(), std::fstream::out|std::fstream::binary); \
    f ## x << pr ## x;                                                  \
    f ## x.close();                                                     \
  } while (0);

  //WF(1, "pbm");
  WF(2, "pgm");
  WF(3, "ppm");
  //WF(4, "pbm");
  WF(5, "pgm");
  WF(6, "ppm");

  //ImageView<bool> p1, p4;
  ImageView<PixelGray<uint8> > p2, p5;
  ImageView<PixelRGB<uint8> >  p3, p6;

  //read_image( p1, fn1 );
  read_image( p2, fn2 );
  read_image( p3, fn3 );
  //read_image( p1, fn4 );
  read_image( p5, fn5 );
  read_image( p6, fn6 );

  EXPECT_EQ( p2.cols(),   1 );
  EXPECT_EQ( p2.rows(),   2 );
  EXPECT_EQ( p2.planes(), 1 );
  EXPECT_EQ( p2(0,0).v(), 12 );
  EXPECT_EQ( p2(0,1).v(), 36 );

  EXPECT_EQ( p5.cols(),   1 );
  EXPECT_EQ( p5.rows(),   2 );
  EXPECT_EQ( p5.planes(), 1 );
  EXPECT_EQ( p5(0,0).v(), 12 );
  EXPECT_EQ( p5(0,1).v(), 36 );

  EXPECT_EQ( p3.cols(),   1 );
  EXPECT_EQ( p3.rows(),   2 );
  EXPECT_EQ( p3.planes(), 1 );
  EXPECT_EQ( p3(0,0).r(), 42 );
  EXPECT_EQ( p3(0,0).g(), 43 );
  EXPECT_EQ( p3(0,0).b(), 44 );
  EXPECT_EQ( p3(0,1).r(), 89 );
  EXPECT_EQ( p3(0,1).g(), 88 );
  EXPECT_EQ( p3(0,1).b(), 87 );

  EXPECT_EQ( p6.cols(),   1 );
  EXPECT_EQ( p6.rows(),   2 );
  EXPECT_EQ( p6.planes(), 1 );
  EXPECT_EQ( p6(0,0).r(), 42 );
  EXPECT_EQ( p6(0,0).g(), 43 );
  EXPECT_EQ( p6(0,0).b(), 44 );
  EXPECT_EQ( p6(0,1).r(), 89 );
  EXPECT_EQ( p6(0,1).g(), 88 );
  EXPECT_EQ( p6(0,1).b(), 87 );
}

TEST( DiskImageResource, PBM_Case_Insentive ) {
  const char pr3[] = "P3 1 2 255 42 43 44 89 88 87";
  UnlinkName caps_p3("test_p3.PPM");

  WF(3, "ppm");
  ImageView<PixelRGB<uint8> > p3a, p3b;
  read_image(p3a, fn3);
  write_image(fn3, p3a);
  read_image(p3b, fn3);

  EXPECT_EQ( p3a.cols(),   p3b.cols() );
  EXPECT_EQ( p3a.rows(),   p3b.rows() );
  EXPECT_EQ( p3a.planes(), p3b.planes() );
  EXPECT_EQ( p3a(0,0).r(), p3b(0,0).r() );
  EXPECT_EQ( p3a(0,0).g(), p3b(0,0).g() );
  EXPECT_EQ( p3a(0,0).b(), p3b(0,0).b() );
  EXPECT_EQ( p3a(0,1).r(), p3b(0,1).r() );
  EXPECT_EQ( p3a(0,1).g(), p3b(0,1).g() );
  EXPECT_EQ( p3a(0,1).b(), p3b(0,1).b() );
}

#undef WF


TEST( DiskImageResource, NonExistentFiles ) {
  boost::scoped_ptr<DiskImageResource> r;

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  EXPECT_THROW(r.reset(DiskImageResourceGDAL::construct_open("nonfile.tif")),
               vw::ArgumentErr);
#endif
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
  EXPECT_THROW(r.reset(DiskImageResourcePNG::construct_open("nonfile.png")),
               vw::ArgumentErr);
#endif
#if defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1
  EXPECT_THROW(r.reset(DiskImageResourceTIFF::construct_open("nonfile.tif")),
               vw::ArgumentErr);
#endif
#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
  EXPECT_THROW(r.reset(DiskImageResourceJPEG::construct_open("nonfile.jpg")),
               vw::ArgumentErr);
#endif
#if defined(VW_HAVE_PKG_OPENEXR) && VW_HAVE_PKG_OPENEXR==1
  EXPECT_THROW(r.reset(DiskImageResourceOpenEXR::construct_open("nonfile.exr")),
               vw::ArgumentErr);
#endif
  EXPECT_THROW(r.reset(DiskImageResourcePDS::construct_open("nonfile.img")),
               vw::ArgumentErr);
  EXPECT_THROW(r.reset(DiskImageResourcePBM::construct_open("nonfile.pgm")),
               vw::ArgumentErr);
}

TEST( DiskImageResource, WrongFiles ) {
  boost::scoped_ptr<DiskImageResource> r;

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
  EXPECT_THROW(r.reset(DiskImageResourceGDAL::construct_open("TestDiskImageResource.cxx")),
               vw::ArgumentErr);
#endif
#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
  EXPECT_THROW(r.reset(DiskImageResourcePNG::construct_open("rgb2x2.jpg")),
               vw::ArgumentErr);
#endif
#if defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1
  EXPECT_THROW(r.reset(DiskImageResourceTIFF::construct_open("rgb2x2.jpg")),
               vw::ArgumentErr);
#endif
#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
  // OSX 10.5 fails to catch this error.
  //EXPECT_THROW(r.reset(DiskImageResourceJPEG::construct_open("rgb2x2.tif")),
  //             vw::ArgumentErr);
#endif
#if defined(VW_HAVE_PKG_OPENEXR) && VW_HAVE_PKG_OPENEXR==1
  EXPECT_THROW(r.reset(DiskImageResourceOpenEXR::construct_open("rgb2x2.tif")),
               vw::ArgumentErr);
#endif
  EXPECT_THROW(r.reset(DiskImageResourcePDS::construct_open("rgb2x2.tif")),
               vw::ArgumentErr);
  EXPECT_THROW(r.reset(DiskImageResourcePBM::construct_open("rgb2x2.tif")),
               vw::ArgumentErr);
}
