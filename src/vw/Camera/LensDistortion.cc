// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Camera/LensDistortion.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Math/LevenbergMarquardt.h>

using namespace vw;

// Special LMA Models to figure out foward and backward ---------

// Optimization functor for computing the undistorted coordinates
// using levenberg marquardt.
struct UndistortOptimizeFunctor : public math::LeastSquaresModelBase<UndistortOptimizeFunctor> {
  typedef Vector2 result_type;
  typedef Vector2 domain_type;
  typedef Matrix<double> jacobian_type;

  const camera::PinholeModel& m_cam;
  const camera::LensDistortion &m_distort;
  UndistortOptimizeFunctor(const camera::PinholeModel& cam, const camera::LensDistortion& d) : m_cam(cam), m_distort(d) {}

  inline result_type operator()( domain_type const& x ) const {
    return m_distort.distorted_coordinates(m_cam, x);
  }
};

struct DistortOptimizeFunctor :  public math::LeastSquaresModelBase<DistortOptimizeFunctor> {
  typedef Vector2 result_type;
  typedef Vector2 domain_type;
  typedef Matrix<double> jacobian_type;

  const camera::PinholeModel& m_cam;
  const camera::LensDistortion &m_distort;
  DistortOptimizeFunctor(const camera::PinholeModel& cam, const camera::LensDistortion& d) : m_cam(cam), m_distort(d) {}
  inline result_type operator()( domain_type const& x ) const {
    Vector2 result = m_distort.undistorted_coordinates(m_cam, x);
    return result;
  }
};

// Backup implemenations for Lens Distortion -------------------

Vector2
vw::camera::LensDistortion::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  UndistortOptimizeFunctor model(cam, *this);
  int status;
  Vector2 solution =
    math::levenberg_marquardt( model, v, v, status, 1e-6, 1e-6 );
  VW_DEBUG_ASSERT( status != math::optimization::eConvergedRelTolerance, PixelToRayErr() << "undistorted_coordinates: failed to converge." );
  return solution;
}

vw::Vector2
vw::camera::LensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  DistortOptimizeFunctor model(cam, *this);
  int status;
  Vector2 solution =
    math::levenberg_marquardt( model, v, v, status, 1e-6, 1e-6 );
  VW_DEBUG_ASSERT( status != math::optimization::eConvergedRelTolerance, PixelToRayErr() << "distorted_coordinates: failed to converge." );
  return solution;
}

// Specific Implementations -------------------------------------

Vector2
vw::camera::TsaiLensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {

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
  Vector2 result(p[0] + bx * du, p[1] + by * dv);
  if (p[0] == cu)
    result[0] =  p[0];
  if (p[1] == cv)
    result[1] =  p[1];

  return result;
}

Vector2
vw::camera::BrownConradyDistortion::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {
  double fu, fv, cu, cv;
  cam.intrinsic_parameters(fu, fv, cu, cv);
  Vector2 offset(cu,cv);
  Vector2 intermediate = p - m_principal_point - offset;
  double r2 = norm_2_sqr(intermediate);
  double radial = 1 + m_radial_distortion[0]*r2 +
    m_radial_distortion[1]*r2*r2 + m_radial_distortion[2]*r2*r2*r2;
  double tangental = m_centering_distortion[0]*r2 + m_centering_distortion[1]*r2*r2;
  intermediate *= radial;
  intermediate[0] -= tangental*sin(m_centering_angle);
  intermediate[1] += tangental*cos(m_centering_angle);
  return intermediate+offset;
}

std::ostream& vw::camera::operator<<(std::ostream & os,
                                     const camera::LensDistortion& ld) {
  ld.write(os);
  return os;
}

