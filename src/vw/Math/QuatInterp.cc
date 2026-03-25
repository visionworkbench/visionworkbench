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
#include <vw/Math/Quaternion.h>
#include <vw/Math/QuatInterp.h>
#include <vw/Math/PositionInterp.h>

#include <cmath>
#include <vector>

using namespace vw;

// QuatLagrangianInterpolationVarTime class.
// Same algorithm as LagrangianInterpolationVarTime but for quaternions.
// The result is normalized after interpolation.
QuatLagrangianInterpolationVarTime::QuatLagrangianInterpolationVarTime
(std::vector<Quat> const& samples, std::vector<double> const& times,
 int radius):
  m_samples(samples), m_times(times), m_radius(radius) {

  VW_ASSERT(m_samples.size() > 1,
            ArgumentErr() << "Expecting at least two samples.\n");
  VW_ASSERT(m_samples.size() == m_times.size(),
            ArgumentErr() << "The number of samples and times must be equal.\n");
  VW_ASSERT(m_radius > 1, ArgumentErr() << "Radius must be > 0.\n");
}

Quat QuatLagrangianInterpolationVarTime::operator()(double t) const {

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

  // Perform the interpolation on each quaternion component
  double q0 = 0, q1 = 0, q2 = 0, q3 = 0;
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
    double w = num_part / denominator;
    q0 += w * m_samples[j].w();
    q1 += w * m_samples[j].x();
    q2 += w * m_samples[j].y();
    q3 += w * m_samples[j].z();
  }

  // Normalize the result
  Quat ans(q0, q1, q2, q3);
  double norm = sqrt(ans.w() * ans.w() + ans.x() * ans.x() +
                     ans.y() * ans.y() + ans.z() * ans.z());
  if (norm > 0)
    ans = Quat(ans.w() / norm, ans.x() / norm,
               ans.y() / norm, ans.z() / norm);

  return ans;
}

//======================================================================
// QuatLagrangianInterpolation class
// interpolation order = 2 * m_radius

// TODO(oalexan1): This needs to be tested and compared with SLERPPoseInterpolation,
// with and without use of splines in the latter.

QuatLagrangianInterpolation::QuatLagrangianInterpolation
(std::vector<Quat> const& samples, double start_time, double time_delta,
 double last_time, int radius):
  m_samples(samples), m_start_time(start_time), m_time_delta(time_delta),
  m_last_time(last_time), m_radius(radius) {

  // Perform a bunch of checks on construction
  VW_ASSERT(m_samples.size() > 1,            ArgumentErr() << "Expecting at least two samples.\n");
  VW_ASSERT(m_radius         > 1,            ArgumentErr() << "Radius must be > 0.\n"            );
  VW_ASSERT(m_time_delta     > 0,            ArgumentErr() << "Time delta must be > 0.\n"        );
  VW_ASSERT(m_last_time      > m_start_time, ArgumentErr() << "Last time must be > start time.\n");

  size_t num_times = round((m_last_time - m_start_time) / m_time_delta) + 1;
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

Quat QuatLagrangianInterpolation::operator()(double t) const {

  // Get the bounding indices
  int    low_i    = static_cast<int>(floor((t - m_start_time) / m_time_delta));
  int    high_i   = low_i + 1;

  VW_ASSERT((low_i >= 0) && (high_i < static_cast<int>(m_samples.size())),
	    ArgumentErr()
            << "XOut of bounds in QuatLagrangianInterpolation for time "
            << t << ".\n");

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
  Quat ans(0, 0, 0, 0);
  for (int j=start; j<=end; ++j) {

    double numerator = 1.0;
    for (int i=start; i<=end; ++i){
      if (i == j)
        continue;
      numerator *= (t - m_times_temp[i-start]);
    }

    double denominator = m_denoms[j-start];
    ans = ans + m_samples[j] * (numerator/denominator);
  }

  return normalize(ans);
}

//======================================================================
// ConstantPoseInterpolation class

ConstantPoseInterpolation::ConstantPoseInterpolation(Quat const& pose) : m_pose(pose) {}

//======================================================================
// SLERPPoseInterpolation class
SLERPPoseInterpolation::SLERPPoseInterpolation(std::vector<Quat> const& pose_samples,
                                               double t0, double dt):
  m_pose_samples(pose_samples), m_t0(t0), m_dt(dt),
  m_tend(m_t0 + m_dt * (m_pose_samples.size() - 1)) {
}

Quat SLERPPoseInterpolation::operator()(double t) const {

  // Make sure that t lies within the range [t0, t0+dt*length(points)]
  VW_ASSERT(t >= m_t0 && t <= m_tend,
	     ArgumentErr() << "Cannot extrapolate point for time "
	     << t << ". Out of valid range. Expecting: "
	     << m_t0 << " <= " << t << " <= " << m_tend << "\n");


  double ratio = (t - m_t0) / m_dt;

  // When t is m_t0 + i * m_dt, this ratio is very close to i, but due to limited
  // numerical precision when doing the floor operation one can get i - 1,
  // or the ceil can be i + 1 which can go out of bounds even though i is within
  // bounds. Hence do a bit of a tune-up here.

  if (std::abs(ratio - round(ratio)) < 1.0e-10)
    ratio = round(ratio);

  int low_i  = (int)floor(ratio);
  int high_i = (int)ceil (ratio);

  double low_t =  m_t0 + m_dt * low_i;
  double norm_t = (t - low_t)/m_dt;

  // norm_t is between 0 and 1 unless it is barely beyond these due to
  // numerical precision issues. The slerp interpolation expects it to
  // not to exceed these bounds.
  if (norm_t < 0.0)
    norm_t = 0.0;
  if (norm_t > 1.0)
    norm_t = 1.0;

  VW_ASSERT(low_i >= 0 && high_i < (int)m_pose_samples.size(),
            ArgumentErr() << "Out of bounds in SLERPPoseInterpolation.\n");

  return vw::math::slerp(norm_t, m_pose_samples[low_i], m_pose_samples[high_i], 0);
}

/// Simple slerp interpolation between a table of pointing directions arranged on a grid.
SlerpGridPointingInterpolation
::SlerpGridPointingInterpolation(std::vector<std::vector<vw::Vector3>> const& directions,
                                 double row0, double drow, double col0, double dcol):
  m_directions(directions), m_row0(row0), m_drow(drow), m_col0(col0), m_dcol(dcol){


  VW_ASSERT(!m_directions.empty() && !m_directions.front().empty(),
	     ArgumentErr() << "Empty input table in SlerpGridPointingInterpolation.\n");

  m_row_end = m_row0 + m_drow * (m_directions.size() - 1);
  m_col_end = m_col0 + m_dcol * (m_directions.front().size() - 1);
}

// Careful here, pix[0] is a column, and pix[1] is a row, so we'll
// access directions(pix[1], pix[0]).
Vector3 SlerpGridPointingInterpolation::operator()(vw::Vector2 const& pix) const {

  double row = pix[1], col = pix[0];
  VW_ASSERT(row >= m_row0 && row <= m_row_end,
    ArgumentErr()
      << "Cannot interpolate for pixel row "
      << row << ". Out of valid range. Expecting "
      << m_row0 << " <= " << row << " <= " << m_row_end << "\n");

  // Calculations for the row
  int low_irow  = (int) floor((row - m_row0) / m_drow);
  int high_irow = (int) ceil ((row - m_row0) / m_drow);
  VW_ASSERT(low_irow >= 0 && high_irow < (int)m_directions.size(),
	     ArgumentErr() << "Out of bounds in SlerpGridPointingInterpolation.\n");
  double low_row  = m_row0 + m_drow * low_irow;
  double norm_row = (row - low_row) / m_drow; // row as fraction of time between points

  // Calculations for the col
  int low_icol  = (int) floor((col - m_col0) / m_dcol);
  int high_icol = (int) ceil ((col - m_col0) / m_dcol);

  // It is convenient to extrapolate for out-of-range columns. This helps guide
  // the solver making use of this towards the solution. This is a bugfix.
  if (col < m_col0) {
    low_icol = 0;
    high_icol = 1;
  } else if (col > m_col_end) {
    low_icol = m_col_end - 1;
    high_icol = m_col_end;
  }

  // We will not use these when out of range in col
  double low_col  = m_col0 + m_dcol * low_icol;
  double norm_col = (col - low_col) / m_dcol;

  vw::Vector3 p;
  Quat L, H;

  {
    // Interpolate the pointing vector for low col
    p = m_directions[low_irow][low_icol];
    Quat ll(0, p[0], p[1], p[2]);
    p = m_directions[high_irow][low_icol];
    Quat hl(0, p[0], p[1], p[2]);
    L = vw::math::slerp(norm_row, ll, hl, 0);
  }

  {
    // Interpolate the pointing vector for high col
    p = m_directions[low_irow][high_icol];
    Quat lh(0, p[0], p[1], p[2]);
    p = m_directions[high_irow][high_icol];
    Quat hh(0, p[0], p[1], p[2]);
    H = vw::math::slerp(norm_row, lh, hh, 0);
  }

  Vector3 result;
  if (m_col0 <= col && col <= m_col_end) {
    // Within range, interpolate
    Quat Res = vw::math::slerp(norm_col, L, H, 0);
    result = Vector3(Res.x(), Res.y(), Res.z());
  } else {
    // Extrapolate. Have to spell out the math per coordinate,
    // as the Quat class does not support this.
    double t = (col - m_col0) / (m_col_end - m_col0);
    result = Vector3((1-t) * L.x() + t * H.x(),
                     (1-t) * L.y() + t * H.y(),
                     (1-t) * L.z() + t * H.z());
  }

  return result;
}

//======================================================================
// SmoothSLERPPoseInterpolation class
// Very speculative. Use with a lot of care.
SmoothSLERPPoseInterpolation::SmoothSLERPPoseInterpolation(std::vector<Quat> const& pose_samples, double t0, double dt, int num_wts, double sigma):
  m_pose_samples(pose_samples), m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_pose_samples.size() - 1)), m_num_wts(num_wts), m_sigma(sigma) {

  VW_ASSERT(m_pose_samples.size() > 1,
	    ArgumentErr() << "Expecting at least two pose samples.\n");
}

Quat SmoothSLERPPoseInterpolation::operator()(double t) const {

  // Make sure that t lies within the range [t0, t0+dt*length(points)]
  VW_ASSERT(t >= m_t0 && t <= m_tend,
	     ArgumentErr() << "Cannot extrapolate point for time "
	     << t << ". Out of valid range. Expecting: "
	     << m_t0 << " <= " << t << " <= " << m_tend << "\n");

  std::vector<double> wts;
  std::vector<int> indices;
  vw::get_wts_and_indices(m_t0, m_dt, t, m_pose_samples.size(),
			  m_num_wts,  m_sigma,
			  wts, indices);

  std::vector<Quat> Q;
  for (size_t i = 0; i < indices.size(); i++)
    Q.push_back(m_pose_samples[indices[i]]);

  return vw::math::slerp_n(wts, Q, 0);
}
