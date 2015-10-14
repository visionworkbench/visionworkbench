// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <vw/Core/Exception.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Quaternion.h>
#include <vw/Math/Vector.h>
#include <vw/Camera/Extrinsics.h>

#include <map>
#include <utility>
#include <vector>

using namespace vw;
using namespace vw::camera;

//======================================================================
// LinearPositionInterpolation class

LinearPositionInterpolation::LinearPositionInterpolation(Vector3 const& initial_position, Vector3 const& velocity_vector) :
  m_initial_position(initial_position), m_velocity_vector(velocity_vector) {}

Vector3
LinearPositionInterpolation::operator()(double t) const {
  return m_initial_position + t * m_velocity_vector;
}

//======================================================================
// LinearPiecewisePositionInterpolation class

LinearPiecewisePositionInterpolation::LinearPiecewisePositionInterpolation( std::vector<Vector3> const& position_samples,
									    double t0, double dt ) :
  m_position(position_samples), m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_position.size() - 1)) {}

Vector3 LinearPiecewisePositionInterpolation::operator()( double t ) const {
  VW_ASSERT( t >= m_t0 && t <= m_tend,
	     ArgumentErr() << "Cannot extrapolate position for time "
	     << t << ". Out of valid range. Expecting " << m_t0 << " <= " << t << " <= " << m_tend << "\n" );

  // Get bounding indices
  size_t low_i  = (size_t) floor( ( t - m_t0 ) / m_dt );
  size_t high_i = (size_t) ceil ( ( t - m_t0 ) / m_dt );

  double low_t  = m_t0 + m_dt * low_i;
  double norm_t = ( t - low_t) / m_dt; // t as fraction of time between points

  Vector3 result = m_position[low_i] + norm_t * ( m_position[high_i] - m_position[low_i] );

  return result;
}

//======================================================================
// PiecewiseAPositionInterpolation class

PiecewiseAPositionInterpolation::PiecewiseAPositionInterpolation( std::vector<Vector3> const& position_samples,
								  std::vector<Vector3> const& velocity_samples,
								  double t0, double dt ) :
  m_position( position_samples ), m_velocity( velocity_samples ),
  m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_position.size() - 1)) {}

Vector3 PiecewiseAPositionInterpolation::operator()( double t ) const {
  VW_ASSERT( t >= m_t0 && t < m_tend,
	     ArgumentErr() << "Cannot extrapolate position for time "
	     << t << ". Out of valid range. Expecting " << m_t0 << " <= " << t << " < " << m_tend << "\n" );

  // Get the bounding indices and the distance from the time at the lower index
  size_t low_i    = (size_t) floor( ( t - m_t0 ) / m_dt );
  size_t high_i   = low_i + 1;
  double offset_t = t - (m_t0 + m_dt * low_i);

  Vector3 a = ( m_velocity[high_i] - m_velocity[low_i] ) / m_dt; // Mean acceleration across the range
  return m_position[low_i] + m_velocity[low_i] * offset_t + a * offset_t * offset_t / 2;
}

//======================================================================
// Curve3DPositionInterpolation class

Curve3DPositionInterpolation::Curve3DPositionInterpolation(std::vector<Vector3> const& position_samples,
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

Vector3 Curve3DPositionInterpolation::operator()( double t ) const {
  Vector3 T(1, t, t*t);
  return m_cached_fit * T;
}

//======================================================================
// HermitePositionInterpolation class

HermitePositionInterpolation::HermitePositionInterpolation( std::vector<Vector3> const& position_samples,
							    std::vector<Vector3> const& velocity_samples,
							    double t0, double dt ) :
  m_position( position_samples ), m_velocity( velocity_samples ),
  m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_position.size() - 1)) {}

Vector3 HermitePositionInterpolation::operator()( double t ) const {
  VW_ASSERT( t >= m_t0 && t < m_tend,
	     ArgumentErr() << "Cannot extrapolate position for time "
	     << t << ". Out of valid range. Expecting " << m_t0 << " <= " << t << " < " << m_tend << "\n");

  size_t low_i = (size_t) floor( ( t - m_t0 ) / m_dt );
  size_t high_i = low_i + 1;

  double low_t = m_t0 + m_dt * low_i;
  double norm_t = ( t - low_t) / m_dt;
  Vector4 poly(1,0,0,0);
  for ( size_t i = 0; i < 3; i++ )
    poly[i+1] = norm_t * poly[i];

  return dot_prod(Vector4(1,0,-3,2), poly) * m_position[low_i] +
    dot_prod(Vector4(0,1,-2,1), poly) * ( m_velocity[low_i] * m_dt ) +
    dot_prod(Vector4(0,0,3,-2), poly) * m_position[high_i] +
    dot_prod(Vector4(0,0,-1,1), poly) * ( m_velocity[high_i] * m_dt );
}

//======================================================================
// ConstantPoseInterpolation class

ConstantPoseInterpolation::ConstantPoseInterpolation(Quat const& pose) : m_pose(pose) {}

//======================================================================
// SLERPPoseInterpolation class

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

SLERPPoseInterpolation::SLERPPoseInterpolation(std::vector<Quat > const& pose_samples, double t0, double dt) :
  m_pose_samples(pose_samples), m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_pose_samples.size() - 1)) {}

Quat SLERPPoseInterpolation::operator()(double t) const {
  // Make sure that t lies within the range [t0, t0+dt*length(points)]
  VW_ASSERT( t >= m_t0 && t <= m_tend,
	     ArgumentErr() << "Cannot extrapolate point for time "
	     << t << ". Out of valid range. Expecting: " << m_t0 << " <= " << t << " <= " << m_tend << "\n");

  size_t low_ind  = (size_t)floor( (t-m_t0) / m_dt );
  size_t high_ind = (size_t)ceil ( (t-m_t0) / m_dt );

  // If there are not enough points to interpolate at the end, we
  // will limit the high_ind here.
  if ( high_ind >= m_pose_samples.size() ) {
    vw_throw( ArgumentErr() << "Attempted to interpolate a quaternion past the last available control point." );
  }

  double low_t =  m_t0 + m_dt * low_ind;
  double norm_t = (t - low_t)/m_dt;

  return this->slerp(norm_t, m_pose_samples[low_ind],
		     m_pose_samples[high_ind], 0);
}

//======================================================================
// LinearTimeInterpolation class

LinearTimeInterpolation::LinearTimeInterpolation( double initial_time, double time_per_line ) :
  m_t0( initial_time ), m_dt(time_per_line) {}

double LinearTimeInterpolation::operator()( double line ) const {
  return m_dt * line + m_t0;
}

//======================================================================
// TLCTimeInterpolation class

TLCTimeInterpolation::TLCTimeInterpolation(std::vector<std::pair<double, double> > const& tlc,
					   double time_offset ) {
  for ( size_t i = 0; i < tlc.size() - 1; i++ ) {
    double t = time_offset + tlc[i].second;
    m_m[tlc[i].first] = ( tlc[i].second - tlc[i+1].second ) / ( tlc[i].first - tlc[i+1].first );
    m_b[tlc[i].first] = t - m_m[tlc[i].first] * tlc[i].first;
  }
}

double TLCTimeInterpolation::operator()( double line ) const {
  map_type::const_iterator m = m_m.lower_bound( line );
  map_type::const_iterator b = m_b.lower_bound( line );
  if ( m != m_m.begin() ) {
    m--; b--;
  }

  return line  * m->second + b->second;
}
