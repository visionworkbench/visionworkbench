#include <vw/Camera/LensDistortion.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Math/LevenbergMarquardt.h>

namespace {
  // Optimization functor for computing the undistorted coordinates using levenberg marquardt.
  struct UndistortOptimizeFunctor : public vw::math::LeastSquaresModelBase<UndistortOptimizeFunctor> {
    typedef vw::Vector2 result_type;
    typedef vw::Vector2 domain_type;
    typedef vw::Matrix<double> jacobian_type;

    const vw::camera::PinholeModel& m_cam;
    const vw::camera::LensDistortion &m_distort;
    UndistortOptimizeFunctor(const vw::camera::PinholeModel& cam, const vw::camera::LensDistortion& d) : m_cam(cam), m_distort(d) {}

    inline result_type operator()( domain_type const& x ) const {
      return m_distort.distorted_coordinates(m_cam, x);
    }
  };
}

std::ostream & vw::camera::operator<<(std::ostream & os, const vw::camera::LensDistortion& ld) {
  ld.write(os);
  return os;
}

vw::Vector2 vw::camera::LensDistortion::undistorted_coordinates(const vw::camera::PinholeModel& cam, vw::Vector2 const& v) const {
  UndistortOptimizeFunctor model(cam, *this);
  int status;
  vw::Vector2 solution = vw::math::levenberg_marquardt( model, v, v, status, 0.1, 0.1 ); // tol = 0.1 pixels
  VW_DEBUG_ASSERT( status != vw::math::optimization::eConvergedRelTolerance, PixelToRayErr() << "undistorted_coordinates: failed to converge." );
  return solution;
}

vw::Vector2 vw::camera::TsaiLensDistortion::distorted_coordinates(const vw::camera::PinholeModel& cam, vw::Vector2 const& p) const {

  double fu, fv, cu, cv;
  cam.intrinsic_parameters(fu, fv, cu, cv);

  double du = p[0] - cu;
  double dv = p[1] - cv;

  if (fu < 1e-300 || fv < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  double x = du / fu;
  double y = dv / fv;

  double x1 = m_distortion[3] / x;
  double y1 = m_distortion[2] / y;

  double r2 = x * x + y * y;

  double x3 = 2.0 * m_distortion[3] * x;
  double y3 = 2.0 * m_distortion[2] * y;

  double b = r2 * (m_distortion[0] + r2 * m_distortion[1]) + x3 + y3;

  double bx = b + r2 * x1;
  double by = b + r2 * y1;

  // Prevent divide by zero at the origin or along the x and y center line
  vw::Vector2 result(p[0] + bx * du, p[1] + by * dv);
  if (p[0] == cu)
    result[0] =  p[0];
  if (p[1] == cv)
    result[1] =  p[1];

  return result;
}
