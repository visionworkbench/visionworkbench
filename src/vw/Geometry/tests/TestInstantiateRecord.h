#include <vw/Geometry.h>
using namespace vw;

const int DimN = 0;
const int IntegralN = 3;
const int SignedN = 3;

typedef double RealT;
typedef int ValT;

typedef Box<double, DimN> ShapeT;
typedef Sphere<double> SphereT;

#include "TestInstantiateRecordList.hh"

class TestInstantiateGeometryRecord : public CxxTest::TestSuite
{
  public: void test_inst() {}
};
