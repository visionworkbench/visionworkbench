// TestCameraCurve.h
#include <cxxtest/TestSuite.h>
#include <vw/Image/ImageView.h>
#include <vw/FileIO.h>
#include <vw/HDR/CameraCurve.h>

using namespace std;
using namespace vw;
using namespace vw::hdr;

class TestCameraCurve : public CxxTest::TestSuite
{
 public:

  void test_curve_generation()
  {
    vnl_matrix<double> pairs = read_matrix("pair_list.exr");
    int degree = 9;

    vnl_real_polynomial poly(degree);

    estimate_camera_curve(pairs, poly, degree);

    TS_TRACE("The resulting coefficients are:");
    TS_TRACE(stringify(poly.coefficients()));

    vnl_matrix<double> graph(255,1);

    for (unsigned int i = 0; i < 255; i++) {
      graph(i,0) = poly.evaluate((double) i / 255.0);
    }

    write_matrix(graph, "graph.exr");
  }
};
