// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Core/Log.h>
#include <vw/Camera/Extrinsics.h>

using namespace vw;
using namespace vw::camera;

Curve3DPositionInterpolation::Curve3DPositionInterpolation(
                                std::vector<Vector3> const& position_samples,
                                double t0, double dt) {
  Matrix<double> Z(position_samples.size()*3, 9);

  Vector<double> p(position_samples.size() * 3);
  // Reshape the position_samples matrix into a column vector
  for (size_t i = 0; i < position_samples.size(); i++) {
    p(3*i) = position_samples[i][0];
    p(3*i+1) = position_samples[i][1];
    p(3*i+2) = position_samples[i][2];
  }

  Vector<double> t(position_samples.size());
  for (size_t i = 0; i < t.size(); i++) {
    t(i) = t0 + dt*i;
  }

  // Populate the Z matrix
  for (size_t i = 0; i < position_samples.size(); i++) {
    Z(3*i  , 0) = 1.0;
    Z(3*i  , 1) = t(i);
    Z(3*i  , 2) = t(i)*t(i);
    Z(3*i+1, 3) = 1.0;
    Z(3*i+1, 4) = t(i);
    Z(3*i+1, 5) = t(i)*t(i);
    Z(3*i+2, 6) = 1.0;
    Z(3*i+2, 7) = t(i);
    Z(3*i+2, 8) = t(i)*t(i);
  }

  Vector<double> x = least_squares(Z,p);
  Matrix3x3 coeff;
  coeff(0,0) = x(0);  coeff(0,1) = x(1); coeff(0,2) = x(2);
  coeff(1,0) = x(3);  coeff(1,1) = x(4); coeff(1,2) = x(5);
  coeff(2,0) = x(6);  coeff(2,1) = x(7); coeff(2,2) = x(8);

  m_cached_fit = coeff;
}

Quat SLERPPoseInterpolation::slerp(double alpha, Quat const& a,
                                   Quat const& b, int spin) const {
  const double SLERP_EPSILON = 1.0E-6;              // a tiny number
  double beta;                      // complementary interp parameter
  double theta;                     // angle between A and B
  double sin_t, cos_t;              // sine, cosine of theta
  double phi;                       // theta plus spins
  int bflip;                        // use negation of B?

  // cosine theta = dot product of A and B
  cos_t = a(1)*b(1) + a(2)*b(2) + a(3)*b(3) + a(0)*b(0);

  // if B is on opposite hemisphere from A, use -B instead
  if (cos_t < 0.0) {
    cos_t = -cos_t;
    bflip = true;
  } else {
    bflip = false;
  }

  // if B is (within precision limits) the same as A,
  // just linear interpolate between A and B.
  // Can't do spins, since we don't know what direction to spin.
  if (1.0 - cos_t < SLERP_EPSILON) {
    beta = 1.0 - alpha;
  } else {                          /* normal case */
    theta = acos(cos_t);
    phi = theta + spin * M_PI;
    sin_t = sin(theta);
    beta = sin(theta - alpha*phi) / sin_t;
    alpha = sin(alpha*phi) / sin_t;
  }

  if (bflip)
    alpha = -alpha;

  // interpolate
  return Quat( beta*a(0) + alpha*b(0),
               beta*a(1) + alpha*b(1),
               beta*a(2) + alpha*b(2),
               beta*a(3) + alpha*b(3) );
}

Quat SLERPPoseInterpolation::operator()(double t) const {
  // Make sure that t lies within the range [t0, t0+dt*length(points)]
  if ((t < m_t0) || (t > m_t0+m_dt*m_pose_samples.size())) {
    vw_out() << "Time: " << t << "   min: " << m_t0
             << "   max: " << (m_t0+m_dt*m_pose_samples.size()) <<"\n";
    vw_throw( ArgumentErr() << "Cannot extrapolate point for time "
              << t << ". Out of valid range." );
  }

  size_t low_ind = (size_t)floor( (t-m_t0) / m_dt );
  size_t high_ind = (size_t)ceil( (t-m_t0) / m_dt );

  // If there are not enough points to interpolate at the end, we
  // will limit the high_ind here.
  if ( high_ind > m_pose_samples.size() ) {
    vw_throw( ArgumentErr() << "Attempted to interpolate a quaternion past the last available control point." );
  } else if (high_ind == m_pose_samples.size()) {
    high_ind = m_pose_samples.size() - 1;
  }

  double low_t =  m_t0 + m_dt * low_ind;
  double norm_t = (t - low_t)/m_dt;

  return this->slerp(norm_t, m_pose_samples[low_ind],
                     m_pose_samples[high_ind], 0);
}
