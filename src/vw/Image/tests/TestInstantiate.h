#include <vw/Image.h>
#include <vw/Image/BlockCacheView.h>

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

typedef int T;

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


#include "TestInstantiateList.hh"

class TestInstantiateMath : public CxxTest::TestSuite
{
  public:
  void test_inst() {}
};
