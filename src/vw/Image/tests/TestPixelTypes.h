// TestPixelTypes.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/PixelTypes.h>

using namespace std;
using namespace vw;

class TestPixelTypes : public CxxTest::TestSuite
{
public:

  void test_pixel_gray()
  {
    // Test default-construction and size with all supported channel types
    { PixelGray<int8> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==1 ); }
    { PixelGray<uint8> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==1 ); }
    { PixelGray<int16> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==2 ); }
    { PixelGray<uint16> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==2 ); }
    { PixelGray<int32> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelGray<uint32> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelGray<int64> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelGray<uint64> p; TS_ASSERT( p.v()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelGray<float32> p; TS_ASSERT( p.v()==0.0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelGray<float64> p; TS_ASSERT( p.v()==0.0 ); TS_ASSERT( sizeof(p)==8 ); }
    // Test channel-value-construction and accessors
    { PixelGray<int8> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<uint8> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<int16> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<uint16> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<int32> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<uint32> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<int64> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<uint64> p(1); TS_ASSERT( p.v()==1 ); TS_ASSERT( p[0]==1 ); TS_ASSERT( p(0)==1 ); }
    { PixelGray<float32> p(1.0); TS_ASSERT( p.v()==1.0 ); TS_ASSERT( p[0]==1.0 ); TS_ASSERT( p(0)==1.0 ); }
    { PixelGray<float64> p(1.0); TS_ASSERT( p.v()==1.0 ); TS_ASSERT( p[0]==1.0 ); TS_ASSERT( p(0)==1.0 ); }
  }

  void test_pixel_graya()
  {
    // Test default-construction and size with all supported channel types
    { PixelGrayA<int8> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==2 ); }
    { PixelGrayA<uint8> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==2 ); }
    { PixelGrayA<int16> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelGrayA<uint16> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelGrayA<int32> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelGrayA<uint32> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelGrayA<int64> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==16 ); }
    { PixelGrayA<uint64> p; TS_ASSERT( p.v()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==16 ); }
    { PixelGrayA<float32> p; TS_ASSERT( p.v()==0.0 && p.a()==0.0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelGrayA<float64> p; TS_ASSERT( p.v()==0.0 && p.a()==0.0 ); TS_ASSERT( sizeof(p)==16 ); }
    // Test channel-value-construction and accessors
    { PixelGrayA<int8> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<uint8> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<int16> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<uint16> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<int32> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<uint32> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<int64> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<uint64> p(1,2); TS_ASSERT( p.v()==1 && p.a()==2 ); TS_ASSERT( p[0]==1 && p[1]==2 ); TS_ASSERT( p(0)==1 && p(1)==2 ); }
    { PixelGrayA<float32> p(1.0,2.0); TS_ASSERT( p.v()==1.0 && p.a()==2.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 ); }
    { PixelGrayA<float64> p(1.0,2.0); TS_ASSERT( p.v()==1.0 && p.a()==2.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 ); }
  }

  void test_pixel_rgb()
  {
    // Test default-construction and size with all supported channel types
    { PixelRGB<int8> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==3 ); }
    { PixelRGB<uint8> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==3 ); }
    { PixelRGB<int16> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==6 ); }
    { PixelRGB<uint16> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==6 ); }
    { PixelRGB<int32> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelRGB<uint32> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelRGB<int64> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==24 ); }
    { PixelRGB<uint64> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 ); TS_ASSERT( sizeof(p)==24 ); }
    { PixelRGB<float32> p; TS_ASSERT( p.r()==0.0 && p.g()==0.0 && p.b()==0.0 ); TS_ASSERT( sizeof(p)==12 ); }
    { PixelRGB<float64> p; TS_ASSERT( p.r()==0.0 && p.g()==0.0 && p.b()==0.0 ); TS_ASSERT( sizeof(p)==24 ); }
    // Test channel-value-construction and accessors
    { PixelRGB<int8> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<uint8> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<int16> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<uint16> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<int32> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<uint32> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<int64> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<uint64> p(1,2,3); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 ); }
    { PixelRGB<float32> p(1.0,2.0,3.0); TS_ASSERT( p.r()==1.0 && p.g()==2.0 && p.b()==3.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 ); }
    { PixelRGB<float64> p(1.0,2.0,3.0); TS_ASSERT( p.r()==1.0 && p.g()==2.0 && p.b()==3.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 ); }
  }

  void test_pixel_rgba()
  {
    // Test default-construction and size with all supported channel types
    { PixelRGBA<int8> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelRGBA<uint8> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==4 ); }
    { PixelRGBA<int16> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelRGBA<uint16> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==8 ); }
    { PixelRGBA<int32> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==16 ); }
    { PixelRGBA<uint32> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==16 ); }
    { PixelRGBA<int64> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==32 ); }
    { PixelRGBA<uint64> p; TS_ASSERT( p.r()==0 && p.g()==0 && p.b()==0 && p.a()==0 ); TS_ASSERT( sizeof(p)==32 ); }
    { PixelRGBA<float32> p; TS_ASSERT( p.r()==0.0 && p.g()==0.0 && p.b()==0.0 && p.a()==0.0 ); TS_ASSERT( sizeof(p)==16 ); }
    { PixelRGBA<float64> p; TS_ASSERT( p.r()==0.0 && p.g()==0.0 && p.b()==0.0 && p.a()==0.0 ); TS_ASSERT( sizeof(p)==32 ); }
    // Test channel-value-construction and accessors
    { PixelRGBA<int8> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<uint8> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<int16> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<uint16> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<int32> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<uint32> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<int64> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<uint64> p(1,2,3,4); TS_ASSERT( p.r()==1 && p.g()==2 && p.b()==3 && p.a()==4 ); TS_ASSERT( p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 ); TS_ASSERT( p(0)==1 && p(1)==2 && p(2)==3 && p(3)==4 ); }
    { PixelRGBA<float32> p(1.0,2.0,3.0,4.0); TS_ASSERT( p.r()==1.0 && p.g()==2.0 && p.b()==3.0 && p.a()==4.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 && p[3]==4.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 && p(3)==4.0 ); }
    { PixelRGBA<float64> p(1.0,2.0,3.0,4.0); TS_ASSERT( p.r()==1.0 && p.g()==2.0 && p.b()==3.0 && p.a()==4.0 ); TS_ASSERT( p[0]==1.0 && p[1]==2.0 && p[2]==3.0 && p[3]==4.0 ); TS_ASSERT( p(0)==1.0 && p(1)==2.0 && p(2)==3.0 && p(3)==4.0 ); }
  }

  void test_rgb_to_gray()
  {
    int val0, val1, val2;
    val0 = val1 = val2 = 40;
   
    // Standard case

    PixelRGB<int8> test_rgb(val0, val1, val2);
    PixelGray<int8> test_gray(test_rgb);
    TS_ASSERT_EQUALS(test_gray.v(),
		     (test_rgb.r() + test_rgb.g() + test_rgb.b())/3 );

    // Case with type incompatibility

    PixelRGB<int16> test_rgb16(val0, val1, val2);
    PixelGray<int8> test_gray8(test_rgb);
    TS_ASSERT_EQUALS(test_gray8.v(),
		     (test_rgb16.r() + test_rgb16.g() + test_rgb16.b())/3 );
  }

  void test_rgb_to_hsv()
  {
    // Standard case (change types and try?)
    PixelRGB<float> input_rgb(100, 100, 100);
    PixelHSV<float> test_hsv(input_rgb);
    TS_ASSERT_EQUALS(test_hsv.h(), 1.0);
    TS_ASSERT_EQUALS(test_hsv.s(), 0);
    TS_ASSERT_EQUALS(test_hsv.v(), 100);
    
    //PixelRGB<uint8> input_rgb8(100, 100, 100);
    //PixelHSV<uint8> test_hsv8(input_rgb8);
    //TS_ASSERT_EQUALS(test_hsv8.h(), 0);
    //TS_ASSERT_EQUALS(test_hsv8.s(), 0);
    //TS_ASSERT_EQUALS(test_hsv8.v(), 100);

    //PixelRGB<uint16> input_rgb16(100, 100, 100);
    //PixelHSV<uint16> test_hsv16(input_rgb16);
    //TS_ASSERT_EQUALS(test_hsv16.h(), 0);
    //TS_ASSERT_EQUALS(test_hsv16.s(), 0);
    //TS_ASSERT_EQUALS(test_hsv16.v(), 100);
  }

  void test_hsv_to_rgb()
  {
    PixelHSV<float> input_hsv(0, 0, 100);
    PixelRGB<float> test_rgb(input_hsv);
    TS_ASSERT_EQUALS(test_rgb.r(), 100);
    TS_ASSERT_EQUALS(test_rgb.g(), 100);
    TS_ASSERT_EQUALS(test_rgb.b(), 100);

    PixelHSV<float> input_hsv_wrap_h(1, 0, 100);
    PixelRGB<float> test_rgb_wrap(input_hsv_wrap_h);
    TS_ASSERT_EQUALS(test_rgb_wrap.r(), 100);
    TS_ASSERT_EQUALS(test_rgb_wrap.g(), 100);
    TS_ASSERT_EQUALS(test_rgb_wrap.b(), 100);

    //PixelHSV<uint8> input_hsv8(0, 0, 100);
    //PixelRGB<uint8> test_rgb8(input_hsv8);
    //TS_ASSERT_EQUALS(test_rgb8.r(), 100);
    //TS_ASSERT_EQUALS(test_rgb8.g(), 100);
    //TS_ASSERT_EQUALS(test_rgb8.b(), 100);

    //PixelHSV<uint16> input_hsv16(0, 0, 100);
    //PixelRGB<uint16> test_rgb16(input_hsv16);
    //TS_ASSERT_EQUALS(test_rgb16.r(), 100);
    //TS_ASSERT_EQUALS(test_rgb16.g(), 100);
    //TS_ASSERT_EQUALS(test_rgb16.b(), 100);

    //PixelHSV<uint32> input_hsv32(0, 0, 100);
    //PixelRGB<uint32> test_rgb32(input_hsv32);
    //TS_ASSERT_EQUALS(test_rgb32.r(), 100);
    //TS_ASSERT_EQUALS(test_rgb32.g(), 100);
    //TS_ASSERT_EQUALS(test_rgb32.b(), 100);
  }

};


