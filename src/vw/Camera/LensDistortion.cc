// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
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
    return m_distort.undistorted_coordinates(m_cam, x);
  }
};

// Backup implemenations for Lens Distortion -------------------

Vector2
vw::camera::LensDistortion::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  UndistortOptimizeFunctor model(cam, *this);
  int status;
  Vector2 solution =
    math::levenberg_marquardt( model, v, v, status, 1e-6, 1e-6, 50 );
  VW_DEBUG_ASSERT( status != math::optimization::eConvergedRelTolerance, PixelToRayErr() << "undistorted_coordinates: failed to converge." );
  return solution;
}

vw::Vector2
vw::camera::LensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& v) const {
  DistortOptimizeFunctor model(cam, *this);
  int status;
  Vector2 solution =
    math::levenberg_marquardt( model, v, v, status, 1e-6, 1e-6, 50 );
  VW_DEBUG_ASSERT( status != math::optimization::eConvergedRelTolerance, PixelToRayErr() << "distorted_coordinates: failed to converge." );
  return solution;
}

// Specific Implementations -------------------------------------

Vector2
vw::camera::TsaiLensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {

  Vector2 focal = cam.focal_length();
  Vector2 offset = cam.point_offset();

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  Vector2 p_0 = elem_quot(p - offset, focal); // represents x and y
  double r2 = norm_2_sqr( p_0 );
  Vector2 distortion( m_distortion[3], m_distortion[2] );
  Vector2 p_1 = elem_quot(distortion, p_0);
  Vector2 p_3 = 2*elem_prod(distortion,p_0);

  Vector2 b =  elem_prod(r2,p_1);
  b = elem_sum(b,r2*(m_distortion[0] + r2 * m_distortion[1]) + sum(p_3));

  // Prevent divide by zero at the origin or along the x and y center line
  Vector2 result = p + elem_prod(b,p - offset);
  if (p[0] == offset[0])
    result[0] = p[0];
  if (p[1] == offset[1])
    result[1] = p[1];

  return result;
}

Vector2
vw::camera::BrownConradyDistortion::undistorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p) const {
  Vector2 offset = cam.point_offset();
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

Vector2
vw::camera::AdjustableTsaiLensDistortion::distorted_coordinates(const camera::PinholeModel& cam, Vector2 const& p )  const {
  Vector2 focal = cam.focal_length();
  Vector2 offset = cam.point_offset();

  if (focal[0] < 1e-300 || focal[1] < 1e-300)
    return Vector2(HUGE_VAL, HUGE_VAL);

  // Create normalized coordinates
  Vector2 p_0 = elem_quot(p - offset, focal); // represents x and y
  double r2 = norm_2_sqr( p_0 );

  // Calculating Radial effects
  double r_n = 1, radial = 0;
  for ( unsigned i = 0; i < m_distortion.size()-3; i++ ) {
    r_n *= r2;
    radial += m_distortion[i]*r_n;
  }

  // Calculating Tangential effects
  Vector2 tangent;
  Vector2 swap_coeff(m_distortion[m_distortion.size()-2],
                     m_distortion[m_distortion.size()-3]);
  tangent += elem_prod(swap_coeff,elem_sum(r2,2*elem_prod(p_0,p_0)));
  tangent += 2*prod(p_0)*subvector(m_distortion,m_distortion.size()-3,2);

  // Final normalized result
  Vector2 result = p_0 + tangent + radial*p_0;

  // Running back through intrinsic matrix (with alpha or skew)
  return elem_prod(result+Vector2(m_distortion[m_distortion.size()-1]*result.y(),0),focal)+offset;
}

std::ostream& vw::camera::operator<<(std::ostream & os,
                                     const camera::LensDistortion& ld) {
  ld.write(os);
  return os;
}

