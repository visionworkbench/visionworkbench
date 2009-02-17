#include <vw/Image.h>

using namespace std;
using namespace vw;

typedef int T;
typedef int T1;
typedef int T2;
typedef int DestT;
typedef int ValueT;

typedef int ChannelT;
typedef PixelRGBA<ChannelT> PixelT;
typedef PixelT Pixel1T;
typedef PixelT Pixel2T;

typedef ImageView<PixelT> ImageT;
typedef ImageT View1T;
typedef ImageT View3T;
typedef ImageT ViewT;
typedef ImageT Image1T;
typedef ImageT SrcT;
typedef ImageT KernelT;

typedef ImageView<PixelMask<PixelT> > MaskViewT;
typedef MaskViewT Image2T;
typedef MaskViewT View2T;

typedef char CharT;
typedef char _CharT;

class TestTransform : public TransformHelper<TestTransform,ConvexFunction,ConvexFunction> {
public:
  TestTransform( Matrix2x2 const& matrix ) {}
  inline Vector2 forward( Vector2 const& p ) const {}
  inline Vector2 reverse( Vector2 const& p ) const {}
  Vector2 operator() (Vector2 a) const {return a;}
};

typedef TestTransform TransformT;
typedef TransformT Tx1T;
typedef TransformT Tx2T;
typedef TransformT Tx3T;
typedef TransformT Tx4T;
typedef TransformT TxT;

typedef ZeroEdgeExtension EdgeT;
typedef EdgeT ExtensionT;
typedef EdgeT EdgeExtensionT;

typedef BilinearInterpolation InterpT;

typedef double ScalarT;
typedef BBox<double, 2> BBoxRealT;

typedef double HighT;
typedef double NewHighT;
typedef double NewLowT;
typedef double OldHighT;
typedef double OldLowT;
typedef double LowT;

typedef TestTransform ImplT;

typedef ViewT::pixel_accessor SrcAccessT;
typedef SrcAccessT KernelAccessT;

struct FuncT : ReturnFixedType<PixelT> {
  PixelT operator() (PixelT a, PixelT b = PixelT(), PixelT c = PixelT()) const {return PixelT();}
  PixelT operator() (MemoryStridingPixelAccessor<PixelT> a) const {return PixelT();}
  BBox2i operator() (BBox2i a) const {return a;}
  void operator() (PixelT& a, vw::PixelMask<PixelT>& b) const {}
  void operator() (PixelT& a, vw::PixelMask<PixelT>& b, PixelT& c) const {}
  BBox2i work_area() const { return BBox2i(); }
};

// ChildPixelT
// ChildT
// ChT
// FuncT
// ImplT
// KernelT
// KRangeT
// MaskViewT
// OutputT
// SourceT
// Src1T
// Src2T
// SrcAccessT
// SrcT
// ThreshT
// ValueT

#if 0
// TODO: HACK for Palette
namespace vw {
  void read_image(vw::ImageView<vw::PixelRGBA<int> >&, const std::string) {}
}

#include <vw/Image/Palette.h>

using namespace vw;

const int SizeN = 5;
const int ForwardType = ConvexFunction;
const int ReverseType = ConvexFunction;


typedef int ChannelT;
typedef ChannelT ChT;
typedef ChannelT NewChT;
typedef ChannelT OldChT;
typedef ChannelT DestT;

typedef PixelRGBA<ChannelT> PixelT;
typedef PixelT DerivedT;
typedef ImageView<PixelT> ViewT;
typedef ViewT ImageT;
typedef ViewT Image1T;
typedef ViewT Image2T;
typedef ViewT Image3T;
typedef PixelT ChildPixel1T;
typedef PixelT ChildPixel2T;
typedef ViewT KernelT;
typedef ZeroEdgeExtension EdgeT;
typedef EdgeT ExtensionT;
typedef BilinearInterpolation InterpT;
typedef int ArgsT;

typedef MinMaxAccumulator<int> AccumT;

typedef ViewT::pixel_accessor PixelAccessorT;
typedef PixelAccessorT ImageIterT;
typedef PixelAccessorT Image1IterT;
typedef PixelAccessorT Image2IterT;
typedef PixelAccessorT Image3IterT;
typedef PixelAccessorT ChildT;
typedef PixelAccessorT IterT;


class TestTransform : public TransformHelper<TestTransform,ConvexFunction,ConvexFunction> {
public:
  TestTransform( Matrix2x2 const& matrix ) {}
  inline Vector2 forward( Vector2 const& p ) const {}
  inline Vector2 reverse( Vector2 const& p ) const {}
  Vector2 operator() (Vector2 a) const {return a;}
};

typedef TestTransform TransformT;
typedef TransformT Tx1T;
typedef TransformT Tx2T;
typedef TransformT TxT;

typedef TransformT ImplT;

struct FuncT : ReturnFixedType<PixelT> {
  PixelT operator() (PixelT a, PixelT b = PixelT(), PixelT c = PixelT()) const {return PixelT();}
  PixelT operator() (MemoryStridingPixelAccessor<PixelT> a) const {return PixelT();}
  BBox2i operator() (BBox2i a) const {return a;}
  BBox2i work_area() const { return BBox2i(); }
};

#endif

#include "TestInstantiateFreeList.hh"

class TestInstantiateImageFree : public CxxTest::TestSuite
{
  public: void test_inst() {}
};
