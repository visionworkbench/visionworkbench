#include <vw/Cartography.h>
using namespace vw;

typedef ImageView<double> TerrainImageT;
typedef ImageView<PixelRGBA<double> > CameraImageT;
typedef ZeroEdgeExtension EdgeT;
typedef BilinearInterpolation InterpT;

#include "TestInstantiateRecordList.hh"

class TestInstantiateCartographyRecord : public CxxTest::TestSuite
{
  public: void test_inst() {}
};
