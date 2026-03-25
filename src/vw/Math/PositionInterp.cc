// __BEGIN_LICENSE__
//  Copyright (c) 2006-2026, United States Government as represented by the
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
#include <vw/Math/Vector.h>
#include <vw/Math/PositionInterp.h>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

using namespace vw;

namespace vw {

  typedef std::pair<double, int> val_index;

  struct sort_descending_by_val {
    bool operator()(const val_index &left, const val_index &right) {
      return left.first > right.first;
    }
  };

  // Given the values ti = t0 + i*dt for i = 0, 1, ..., num_samples - 1,
  // and the value t, compute the weights
  // w_i=exp(-sigma*(t-ti)^2/dt^2). Keep the largest m_num_wts weights.
  // Normalize the weights to add to 1 and return them. We'll use them
  // later for interpolation. This works even when t is out of range.
  void get_wts_and_indices(double t0, double dt, double t, int num_samples,
			   int num_wts, double sigma,
			   std::vector<double> & wts, std::vector<int> & indices){

    std::vector<val_index> V(num_samples);
    for (int i = 0; i < num_samples; i++) {
      double ratio = (t0 + dt * i - t)/dt;
      V[i] = val_index(exp(-sigma*ratio*ratio), i);
    }

    std::sort(V.begin(), V.end(), sort_descending_by_val());
    V.resize(std::min(num_samples, num_wts));

    int len = V.size();
    double sum = 0.0;
    wts.resize(len);
    indices.resize(len);
    for (int i = 0; i < len; i++) {
      wts[i] = V[i].first;
      sum += wts[i];
      indices[i] = V[i].second;
    }

    // Normalize the weights
    for (int i = 0; i < len; i++) {
      wts[i] /= sum;
    }

  }

} // namespace vw

//======================================================================
// LinearPositionInterpolation class

LinearPositionInterpolation::LinearPositionInterpolation(Vector3 const& initial_position,
                                                         Vector3 const& initial_velocity) :
  m_initial_position(initial_position), m_initial_velocity(initial_velocity) {}

Vector3
LinearPositionInterpolation::operator()(double t) const {
  return m_initial_position + t * m_initial_velocity;
}

//======================================================================
// LinearPiecewisePositionInterpolation class

LinearPiecewisePositionInterpolation::LinearPiecewisePositionInterpolation
(std::vector<Vector3> const& position_samples,
 double t0, double dt) :
  m_position_samples(position_samples), m_t0(t0), m_dt(dt),
  m_tend(m_t0 + m_dt * (m_position_samples.size() - 1)) {}

Vector3 LinearPiecewisePositionInterpolation::operator()(double t) const {

  VW_ASSERT(t >= m_t0 && t <= m_tend + 1.0e-8,
            ArgumentErr() << "Cannot extrapolate position for time "
            << t << ". Out of valid range. Expecting " << m_t0 << " <= " << t << " <= "
            << m_tend << "\n");

  // Get bounding indices
  int low_i  = (int) floor((t - m_t0) / m_dt);
  int high_i = (int) ceil ((t - m_t0) / m_dt);

  // Adjustment for when t equals to m_tend within numerical
  // precision, when the result should be the last element.  If in
  // doubt, extrapolate a bit rather than fail, so t can be allowed to
  // be a tiny bit more above t_end.
  if (std::abs(t - m_tend) < 1.0e-8 && high_i == (int)m_position_samples.size())
    high_i = low_i;

  VW_ASSERT(low_i >= 0 && high_i < (int)m_position_samples.size(),
	     ArgumentErr() << "Out of bounds in LinearPiecewisePositionInterpolation.\n");

  double low_t  = m_t0 + m_dt * low_i;
  double norm_t = (t - low_t) / m_dt; // t as fraction of time between points

  Vector3 result = m_position_samples[low_i]
    + norm_t * (m_position_samples[high_i] - m_position_samples[low_i]);

  return result;
}

//======================================================================
// SmoothPiecewisePositionInterpolation class

SmoothPiecewisePositionInterpolation::SmoothPiecewisePositionInterpolation
(std::vector<Vector3> const& position_samples, double t0, double dt, int num_wts, double sigma):
  m_position_samples(position_samples), m_t0(t0), m_dt(dt),
  m_tend(m_t0 + m_dt * (m_position_samples.size() - 1)), m_num_wts(num_wts), m_sigma(sigma) {

  VW_ASSERT(m_position_samples.size() > 1,
	    ArgumentErr() << "Expecting at least two position samples.\n");
}

Vector3 SmoothPiecewisePositionInterpolation::operator()(double t) const {

  // Make sure that t lies within the range [t0, t0+dt*length(points)].
  VW_ASSERT(t >= m_t0 && t <= m_tend,
	     ArgumentErr() << "Cannot extrapolate point for time "
	     << t << ". Out of valid range. Expecting: "
	     << m_t0 << " <= " << t << " <= " << m_tend << "\n");

  std::vector<double> wts;
  std::vector<int> indices;
  vw::get_wts_and_indices(m_t0, m_dt, t, m_position_samples.size(),
			  m_num_wts, m_sigma,
			  wts, indices);

  Vector3 ans;
  for (size_t i = 0; i < indices.size(); i++)
    ans += wts[i]*m_position_samples[indices[i]];

  return ans;
}

// Get the indices corresponding to the largest weights
std::vector<int> SmoothPiecewisePositionInterpolation::get_indices_of_largest_weights(double t)
  const {

  // Make sure that t lies within the range [t0, t0+dt*length(points)]
  VW_ASSERT(t >= m_t0 && t <= m_tend,
	     ArgumentErr() << "Cannot extrapolate point for time "
	     << t << ". Out of valid range. Expecting: "
	     << m_t0 << " <= " << t << " <= " << m_tend << "\n");

  std::vector<double> wts;
  std::vector<int> indices;
  vw::get_wts_and_indices(m_t0, m_dt, t, m_position_samples.size(),
			  m_num_wts, m_sigma,
			  wts, indices);

  return indices;
}

//======================================================================
// LagrangianInterpolationVarTime class
LagrangianInterpolationVarTime::LagrangianInterpolationVarTime
(std::vector<Vector3> const& samples, std::vector<double> const& times,
 int radius):
  m_samples(samples), m_times(times), m_radius(radius) {

  VW_ASSERT(m_samples.size() > 1,
            ArgumentErr() << "Expecting at least two samples.\n");
  VW_ASSERT(m_samples.size() == m_times.size(),
            ArgumentErr() << "The number of samples and times must be equal.\n");
  VW_ASSERT(m_radius > 1, ArgumentErr() << "Radius must be > 0.\n");
}

Vector3 LagrangianInterpolationVarTime::operator()(double t) const {

  // Find where t lies in our list of samples
  const int num_samples = static_cast<int>(m_times.size());
  int next = -1;
  for (int i = 0; i < num_samples; ++i) {
    if (m_times[i] > t) {
      next = i;
      break;
    }
  }
  // If t is at or past the last sample, set next to the last index
  if (next < 0)
    next = num_samples - 1;

  // Clamp the interpolation window to stay within bounds
  int start = next - m_radius;
  int end   = next + m_radius; // the last index used is end - 1
  if (start < 0) {
    end   -= start;
    start  = 0;
  }
  if (end > num_samples) {
    start -= (end - num_samples);
    end    = num_samples;
  }
  VW_ASSERT((start >= 0) && (end <= num_samples),
            ArgumentErr() << "Not enough samples to interpolate time "
                          << t << "\n");

  // Perform the interpolation
  Vector3 ans;
  for (int j = start; j < end; ++j) {
    double num_part = 1.0, denominator = 1.0;
    for (int i = start; i < end; ++i) {
      if (i == j)
        continue;
      num_part *= (t - m_times[i]);
    }
    // TODO(oalexan1): the denominator could be cached
    for (int i = start; i < end; ++i) {
      if (i == j)
        continue;
      denominator *= (m_times[j] - m_times[i]);
    }
    ans += m_samples[j] * (num_part / denominator);
  }

  return ans;
}

// LagrangianInterpolation class
// interpolation order = 2 * m_radius
LagrangianInterpolation::LagrangianInterpolation
(std::vector<Vector3> const& samples, double start_time, double time_delta,
 double last_time, int radius):
  m_samples(samples), m_start_time(start_time), m_time_delta(time_delta),
  m_last_time(last_time), m_radius(radius) {

  // Perform a bunch of checks on construction
  VW_ASSERT(m_samples.size() > 1,            ArgumentErr() << "Expecting at least two samples.\n");
  VW_ASSERT(m_radius         > 1,            ArgumentErr() << "Radius must be > 0.\n"            );
  VW_ASSERT(m_time_delta     > 0,            ArgumentErr() << "Time delta must be > 0.\n"        );
  VW_ASSERT(m_last_time      > m_start_time, ArgumentErr() << "Last time must be > start time.\n");

  size_t num_times = round((m_last_time - m_start_time) / m_time_delta)+1;
  VW_ASSERT(m_samples.size() == num_times,
	    ArgumentErr() << "The number of samples and times must be equal.\n");

	// We can precalculate the denominator of the equation here
	//  since the time intervals between the data points are constant.
	const int num_points = 2*radius; // Number of points used in each calculation
	m_denoms.resize(num_points);
	m_times_temp.resize(num_points);

  for (int j=0; j<num_points; ++j) {

    double denominator = 1.0;
    for (int i=0; i<num_points; ++i){
      if (i == j)
        continue;
      denominator *= (j-i)*m_time_delta;
    }
    m_denoms[j] = denominator;
  } // End outer loop
}

Vector3 LagrangianInterpolation::operator()(double t) const {

  // Get the bounding indices
  int    low_i    = static_cast<int>(floor((t - m_start_time) / m_time_delta));
  int    high_i   = low_i + 1;

  VW_ASSERT((low_i >= 0) && (high_i < static_cast<int>(m_samples.size())),
	    ArgumentErr() << "Out of bounds in LagrangianInterpolation for time " << t << ".\n");

  // Check that we have enough bordering points to interpolate
  int start = low_i  - (m_radius-1);
  int end   = high_i + (m_radius-1);

  if (start < 0) {
    // Have to use points more on the right
    int shift = -start;
    start += shift;
    end   += shift;
  }

  if (end >= static_cast<int>(m_samples.size())) {
    // Have to use more points on the left
    int shift = end - static_cast<int>(m_samples.size()) + 1;
    start -= shift;
    end -= shift;
  }

  VW_ASSERT((start >= 0) && (end < static_cast<int>(m_samples.size())),
	    ArgumentErr() << "Not enough samples to interpolate time " << t << ".\n");

  // Compute the times of the points being used for interpolation
  m_times_temp[0] = m_start_time + start*m_time_delta;
  for (size_t k=1; k<m_times_temp.size(); ++k)
    m_times_temp[k] = m_times_temp[k-1] + m_time_delta;

  // Perform the interpolation. The interval [start, end] has 2 *
  // m_radius values.  We end up multiplying 2 * m_radius - 1 values
  // for each numerator.
  Vector3 ans(0,0,0);
  for (int j=start; j<=end; ++j) {

    double numerator = 1.0;
    for (int i=start; i<=end; ++i){
      if (i == j)
        continue;
      numerator *= (t - m_times_temp[i-start]);
    }

    double denominator = m_denoms[j-start];
    ans += m_samples[j] * (numerator/denominator);
  }

  return ans;
}

//======================================================================
// PiecewiseAPositionInterpolation class

PiecewiseAPositionInterpolation::PiecewiseAPositionInterpolation
(std::vector<Vector3> const& position_samples,
 std::vector<Vector3> const& velocity_samples,
 double t0, double dt):
  m_position_samples(position_samples), m_velocity(velocity_samples),
  m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_position_samples.size() - 1)) {}

Vector3 PiecewiseAPositionInterpolation::operator()(double t) const {

  VW_ASSERT(t >= m_t0 && t <= m_tend + 1.0e-8,
            ArgumentErr() << "Cannot extrapolate position for time "
            << t << ". Out of valid range. Expecting " << m_t0 << " <= "
            << t << " < " << m_tend << ".\n");

  // Get the bounding indices and the distance from the time at the lower index
  int low_i    = (int) floor((t - m_t0) / m_dt);
  int high_i   = low_i + 1;
  double offset_t = t - (m_t0 + m_dt * low_i);

  // Adjustment for when t equals to m_tend within numerical
  // precision, when the result should be the last element.  If in
  // doubt, extrapolate a bit rather than fail, so t can be allowed to
  // be a tiny bit more above t_end.
  if (std::abs(t - m_tend) < 1.0e-8 && high_i == (int)m_position_samples.size())
    high_i = low_i;

  VW_ASSERT(low_i >= 0 && high_i < (int)m_position_samples.size(),
	     ArgumentErr() << "Out of bounds in PiecewiseAPositionInterpolation.\n");

  // Mean acceleration across the range
  Vector3 a = (m_velocity[high_i] - m_velocity[low_i]) / m_dt;

  return m_position_samples[low_i] + m_velocity[low_i] * offset_t + a * offset_t * offset_t / 2;
}

//======================================================================
// Curve3DPositionInterpolation class

Curve3DPositionInterpolation::Curve3DPositionInterpolation
(std::vector<Vector3> const& position_samples,
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

Vector3 Curve3DPositionInterpolation::operator()(double t) const {
  Vector3 T(1, t, t*t);
  return m_cached_fit * T;
}

//======================================================================
// HermitePositionInterpolation class

HermitePositionInterpolation::HermitePositionInterpolation
(std::vector<Vector3> const& position_samples,
 std::vector<Vector3> const& velocity_samples,
 double t0, double dt) :
  m_position_samples(position_samples), m_velocity(velocity_samples),
  m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_position_samples.size() - 1)) {}

Vector3 HermitePositionInterpolation::operator()(double t) const {

  VW_ASSERT(t >= m_t0 && t < m_tend,
            ArgumentErr() << "Cannot extrapolate position for time "
            << t << ". Out of valid range. Expecting "
            << m_t0 << " <= " << t << " < " << m_tend << "\n");

  int low_i = (int) floor((t - m_t0) / m_dt);
  int high_i = low_i + 1;

  VW_ASSERT(low_i >= 0 && high_i < (int)m_position_samples.size(),
	     ArgumentErr() << "Out of bounds in HermitePositionInterpolation.\n");

  double low_t = m_t0 + m_dt * low_i;
  double norm_t = (t - low_t) / m_dt;
  Vector4 poly(1,0,0,0);
  for (size_t i = 0; i < 3; i++)
    poly[i+1] = norm_t * poly[i];

  return dot_prod(Vector4(1,0,-3,2), poly) * m_position_samples[low_i] +
    dot_prod(Vector4(0,1,-2,1), poly) * (m_velocity[low_i] * m_dt) +
    dot_prod(Vector4(0,0,3,-2), poly) * m_position_samples[high_i] +
    dot_prod(Vector4(0,0,-1,1), poly) * (m_velocity[high_i] * m_dt);
}
