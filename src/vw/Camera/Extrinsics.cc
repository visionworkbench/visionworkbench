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


namespace vw { namespace camera {

  // Given the values ti = t0 + i*dt for i = 0, 1, ..., num_samples - 1,
  // and the value t, compute the weights
  // w_i=exp(-sigma*(t-ti)^2/dt^2). Keep the largest m_num_wts weights.
  // Normalize the weights to add to 1 and return them. We'll use them
  // later for interpolation. This works even when t is out of range.
  
  // TODO: Move this somewhere else!

  typedef std::pair<double, int> val_index;

  struct sort_descending_by_val {
    bool operator()(const val_index &left, const val_index &right) {
      return left.first > right.first;
    }
  };

  void get_wts_and_indices(double t0, double dt, double t, int num_samples,
			   int num_wts, double sigma,
			   std::vector<double> & wts, std::vector<int> & indices){

    std::vector<val_index> V(num_samples);
    for (int i = 0; i < num_samples; i++) {
      double ratio = (t0 + dt * i - t)/dt;
      V[i] = val_index( exp(-sigma*ratio*ratio), i );
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

}}

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

LinearPiecewisePositionInterpolation::LinearPiecewisePositionInterpolation( std::vector<Vector3> const& position_samples,
									    double t0, double dt ) :
  m_position_samples(position_samples), m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_position_samples.size() - 1)) {}

Vector3 LinearPiecewisePositionInterpolation::operator()( double t ) const {

  // Make sure that t lies within the range [t0, t0+dt*length(points)]
  VW_ASSERT( t >= m_t0 && t <= m_tend,
	     ArgumentErr() << "Cannot extrapolate position for time "
	     << t << ". Out of valid range. Expecting " << m_t0 << " <= " << t << " <= " << m_tend << "\n" );

  // Get bounding indices
  int low_i  = (int) floor( ( t - m_t0 ) / m_dt );
  int high_i = (int) ceil ( ( t - m_t0 ) / m_dt );

  VW_ASSERT( low_i >= 0 && high_i < (int)m_position_samples.size(),
	     ArgumentErr() << "Out of bounds in LinearPiecewisePositionInterpolation.\n" );

  double low_t  = m_t0 + m_dt * low_i;
  double norm_t = ( t - low_t) / m_dt; // t as fraction of time between points

  Vector3 result = m_position_samples[low_i] + norm_t * ( m_position_samples[high_i] - m_position_samples[low_i] );

  return result;
}

//======================================================================
// SmoothPiecewisePositionInterpolation class

SmoothPiecewisePositionInterpolation::SmoothPiecewisePositionInterpolation
(std::vector<Vector3> const& position_samples, double t0, double dt, int num_wts, double sigma):
  m_position_samples(position_samples), m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_position_samples.size() - 1)), m_num_wts(num_wts), m_sigma(sigma) {

  VW_ASSERT(m_position_samples.size() > 1,
	    ArgumentErr() << "Expecting at least two position samples.\n" );
}

Vector3 SmoothPiecewisePositionInterpolation::operator()( double t ) const {

  // Make sure that t lies within the range [t0, t0+dt*length(points)]
  VW_ASSERT( t >= m_t0 && t <= m_tend,
	     ArgumentErr() << "Cannot extrapolate point for time "
	     << t << ". Out of valid range. Expecting: "
	     << m_t0 << " <= " << t << " <= " << m_tend << "\n");

  std::vector<double> wts;
  std::vector<int> indices;
  vw::camera::get_wts_and_indices(m_t0, m_dt, t, m_position_samples.size(),
				  m_num_wts, m_sigma,
				  wts, indices);

  Vector3 ans;
  for (size_t i = 0; i < indices.size(); i++)
    ans += wts[i]*m_position_samples[indices[i]];

  return ans;
}

// Get the indices corresponding to the largest weights
std::vector<int> SmoothPiecewisePositionInterpolation::get_indices_of_largest_weights(double t) const {

  // Make sure that t lies within the range [t0, t0+dt*length(points)]
  VW_ASSERT( t >= m_t0 && t <= m_tend,
	     ArgumentErr() << "Cannot extrapolate point for time "
	     << t << ". Out of valid range. Expecting: "
	     << m_t0 << " <= " << t << " <= " << m_tend << "\n");

  std::vector<double> wts;
  std::vector<int> indices;
  vw::camera::get_wts_and_indices(m_t0, m_dt, t, m_position_samples.size(),
				  m_num_wts, m_sigma,
				  wts, indices);

  return indices;
}

//======================================================================
// LagrangianInterpolationVarTime class

LagrangianInterpolationVarTime::LagrangianInterpolationVarTime
(std::vector<Vector3> const& samples, std::vector<double> const& times, int radius):
  m_samples(samples), m_times(times), m_radius(radius) {

  VW_ASSERT(m_samples.size() > 1, ArgumentErr() << "Expecting at least two samples.\n" );
  VW_ASSERT(m_samples.size() == m_times.size(),
	    ArgumentErr() << "The number of samples and times must be equal.\n" );
  VW_ASSERT(m_radius > 1, ArgumentErr() << "Radius must be > 0.\n" );
}

Vector3 LagrangianInterpolationVarTime::operator()( double t ) const {

  // Find where t lies in our list of samples
  const int num_samples = static_cast<int>(m_times.size());
  int next = -1;
  for (int i=0; i<num_samples; ++i) {
    if (m_times[i] > t) {
      next = i;
      break;
    }
  }
  
  // Check that we have enough bordering points to interpolate
  int start = next - m_radius;
  int end   = next + m_radius; // Note: The last index we use is end-1!
  VW_ASSERT((start >= 0) && (end <= num_samples),
	    ArgumentErr() << "Not enough samples to interpolate time " << t << "\n" );
    
  // Perform the interpolation
  Vector3 ans;
  
  for (int j=start; j<end; ++j) {
    double  num_part=1.0, denominator=1.0;
    // Numerator
    for (int i=start; i<end; ++i){
      if (i == j)
        continue;
      num_part *= (t - m_times[i]);
    }

    // Denominator
    for (int i=start; i<end; ++i){
      if (i == j)
        continue;
      denominator *= (m_times[j] - m_times[i]);
    }
    
    ans += m_samples[j] * (num_part/denominator);
  }

  return ans;
}

//======================================================================
// LagrangianInterpolation class

LagrangianInterpolation::LagrangianInterpolation
(std::vector<Vector3> const& samples, double start_time, double time_delta, double last_time, int radius):
  m_samples(samples), m_start_time(start_time), m_time_delta(time_delta), 
  m_last_time(last_time), m_radius(radius) {

  // Perform a bunch of checks on construction
  VW_ASSERT(m_samples.size() > 1,            ArgumentErr() << "Expecting at least two samples.\n" );
  VW_ASSERT(m_radius         > 1,            ArgumentErr() << "Radius must be > 0.\n"             );
  VW_ASSERT(m_time_delta     > 0,            ArgumentErr() << "Time delta must be > 0.\n"         );
  VW_ASSERT(m_last_time      > m_start_time, ArgumentErr() << "Last time must be > start time.\n" );
  
  size_t num_times = round((m_last_time - m_start_time) / m_time_delta)+1;
  VW_ASSERT(m_samples.size() == num_times,
	    ArgumentErr() << "The number of samples and times must be equal.\n" );
	
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

Vector3 LagrangianInterpolation::operator()( double t ) const {

  // Get the bounding indices
  int    low_i    = static_cast<int>(floor( (t - m_start_time) / m_time_delta ));
  int    high_i   = low_i + 1;

  // Check that we have enough bordering points to interpolate
  int start = low_i  - (m_radius-1);
  int end   = high_i + (m_radius-1);
  VW_ASSERT((start >= 0) && (end < static_cast<int>(m_samples.size())),
	    ArgumentErr() << "Not enough samples to interpolate time " << t << "\n" );
  
  // Compute the times of the points being used for interpolation
  m_times_temp[0] = m_start_time + start*m_time_delta;
  for (size_t k=1; k<m_times_temp.size(); ++k)
    m_times_temp[k] = m_times_temp[k-1] + m_time_delta;
    
  // Perform the interpolation
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

PiecewiseAPositionInterpolation::PiecewiseAPositionInterpolation( std::vector<Vector3> const& position_samples,
								  std::vector<Vector3> const& velocity_samples,
								  double t0, double dt ) :
  m_position_samples( position_samples ), m_velocity( velocity_samples ),
  m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_position_samples.size() - 1)) {}

Vector3 PiecewiseAPositionInterpolation::operator()( double t ) const {

  VW_ASSERT( t >= m_t0 && t < m_tend,
	     ArgumentErr() << "Cannot extrapolate position for time "
	     << t << ". Out of valid range. Expecting " << m_t0 << " <= " << t << " < " << m_tend << "\n" );

  // Get the bounding indices and the distance from the time at the lower index
  int low_i    = (int) floor( ( t - m_t0 ) / m_dt );
  int high_i   = low_i + 1;
  double offset_t = t - (m_t0 + m_dt * low_i);

  VW_ASSERT( low_i >= 0 && high_i < (int)m_position_samples.size(),
	     ArgumentErr() << "Out of bounds in PiecewiseAPositionInterpolation.\n" );

  Vector3 a = ( m_velocity[high_i] - m_velocity[low_i] ) / m_dt; // Mean acceleration across the range
  return m_position_samples[low_i] + m_velocity[low_i] * offset_t + a * offset_t * offset_t / 2;
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
  m_position_samples( position_samples ), m_velocity( velocity_samples ),
  m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_position_samples.size() - 1)) {}

Vector3 HermitePositionInterpolation::operator()( double t ) const {

  VW_ASSERT( t >= m_t0 && t < m_tend,
	     ArgumentErr() << "Cannot extrapolate position for time "
	     << t << ". Out of valid range. Expecting " << m_t0 << " <= " << t << " < " << m_tend << "\n");

  int low_i = (int) floor( ( t - m_t0 ) / m_dt );
  int high_i = low_i + 1;

  VW_ASSERT( low_i >= 0 && high_i < (int)m_position_samples.size(),
	     ArgumentErr() << "Out of bounds in HermitePositionInterpolation.\n" );

  double low_t = m_t0 + m_dt * low_i;
  double norm_t = ( t - low_t) / m_dt;
  Vector4 poly(1,0,0,0);
  for ( size_t i = 0; i < 3; i++ )
    poly[i+1] = norm_t * poly[i];

  return dot_prod(Vector4(1,0,-3,2), poly) * m_position_samples[low_i] +
    dot_prod(Vector4(0,1,-2,1), poly) * ( m_velocity[low_i] * m_dt ) +
    dot_prod(Vector4(0,0,3,-2), poly) * m_position_samples[high_i] +
    dot_prod(Vector4(0,0,-1,1), poly) * ( m_velocity[high_i] * m_dt );
}

//======================================================================
// ConstantPoseInterpolation class

ConstantPoseInterpolation::ConstantPoseInterpolation(Quat const& pose) : m_pose(pose) {}



//======================================================================
// SLERPPoseInterpolation class
SLERPPoseInterpolation::SLERPPoseInterpolation(std::vector<Quat> const& pose_samples,
                                               double t0, double dt) :
  m_pose_samples(pose_samples), m_t0(t0), m_dt(dt),
  m_tend(m_t0 + m_dt * (m_pose_samples.size() - 1)) {}

Quat SLERPPoseInterpolation::operator()(double t) const {

  // Make sure that t lies within the range [t0, t0+dt*length(points)]
  VW_ASSERT( t >= m_t0 && t <= m_tend,
	     ArgumentErr() << "Cannot extrapolate point for time "
	     << t << ". Out of valid range. Expecting: "
	     << m_t0 << " <= " << t << " <= " << m_tend << "\n");

  int low_i  = (int)floor( (t-m_t0) / m_dt );
  int high_i = (int)ceil ( (t-m_t0) / m_dt );

  VW_ASSERT( low_i >= 0 && high_i < (int)m_pose_samples.size(),
	     ArgumentErr() << "Out of bounds in SLERPPoseInterpolation.\n" );

  double low_t =  m_t0 + m_dt * low_i;
  double norm_t = (t - low_t)/m_dt;

  return vw::math::slerp(norm_t, m_pose_samples[low_i], m_pose_samples[high_i], 0);
}

/// Simple slerp interpolation between a table of pointing directions arranged on a grid.
SlerpGridPointingInterpolation
::SlerpGridPointingInterpolation(std::vector< std::vector<vw::Vector3> > const& directions,
                                 double row0, double drow, double col0, double dcol):
  m_directions(directions), m_row0(row0), m_drow(drow), m_col0(col0), m_dcol(dcol){
  

  VW_ASSERT( !m_directions.empty() && !m_directions.front().empty(),
	     ArgumentErr() << "Empty input table in SlerpGridPointingInterpolation.\n" );

  m_row_end = m_row0 + m_drow * (m_directions.size() - 1);
  m_col_end = m_col0 + m_dcol * (m_directions.front().size() - 1);
}

// Careful here, pix[0] is a column, and pix[1] is a row, so we'll
// access directions(pix[1], pix[0]).
Vector3 SlerpGridPointingInterpolation::operator()(vw::Vector2 const& pix) const {

  double row = pix[1], col = pix[0];
  VW_ASSERT( row >= m_row0 && row <= m_row_end,
	     ArgumentErr() << "Cannot interpolate for pixel row "
	     << row << ". Out of valid range. Expecting "
             << m_row0 << " <= " << row << " <= " << m_row_end << "\n" );
  VW_ASSERT( col >= m_col0 && col <= m_col_end,
	     ArgumentErr() << "Cannot interpolate for pixel col "
	     << col << ". Out of valid range. Expecting "
             << m_col0 << " <= " << col << " <= " << m_col_end << "\n" );

  // Calculations for the row
  int low_irow  = (int) floor( ( row - m_row0 ) / m_drow );
  int high_irow = (int) ceil ( ( row - m_row0 ) / m_drow );
  VW_ASSERT( low_irow >= 0 && high_irow < (int)m_directions.size(),
	     ArgumentErr() << "Out of bounds in SlerpGridPointingInterpolation.\n" );
  double low_row  = m_row0 + m_drow * low_irow;
  double norm_row = ( row - low_row) / m_drow; // row as fraction of time between points

  // Calculations for the col
  int low_icol  = (int) floor( ( col - m_col0 ) / m_dcol );
  int high_icol = (int) ceil ( ( col - m_col0 ) / m_dcol );
  VW_ASSERT( low_icol >= 0 && high_icol < (int)m_directions.front().size(),
	     ArgumentErr() << "Out of bounds in SlerpGridPointingInterpolation.\n" );
  double low_col  = m_col0 + m_dcol * low_icol;
  double norm_col = ( col - low_col) / m_dcol; // col as fraction of time between points

  vw::Vector3 p;
  Quat L, H;

  {
    // Interpolate the pointing vector for low_col
    p = m_directions[low_irow][low_icol];
    Quat ll(0, p[0], p[1], p[2]);
    p = m_directions[high_irow][low_icol];
    Quat hl(0, p[0], p[1], p[2]);
    L = vw::math::slerp(norm_row, ll, hl, 0);
  }

  {
    // Interpolate the pointing vector for high_col
    p = m_directions[low_irow][high_icol];
    Quat lh(0, p[0], p[1], p[2]);
    p = m_directions[high_irow][high_icol];
    Quat hh(0, p[0], p[1], p[2]);
    H = vw::math::slerp(norm_row, lh, hh, 0);
  }

  // Now interpolate between low_col and high_col
  Quat Res = vw::math::slerp(norm_col, L, H, 0);

  Vector3 result( Res.x(), Res.y(), Res.z());
  return result;
}

//======================================================================
// SmoothSLERPPoseInterpolation class
// Very speculative. Use with a lot of care. 
SmoothSLERPPoseInterpolation::SmoothSLERPPoseInterpolation(std::vector<Quat> const& pose_samples, double t0, double dt, int num_wts, double sigma):
  m_pose_samples(pose_samples), m_t0(t0), m_dt(dt), m_tend(m_t0 + m_dt * (m_pose_samples.size() - 1)), m_num_wts(num_wts), m_sigma(sigma) {

  VW_ASSERT(m_pose_samples.size() > 1,
	    ArgumentErr() << "Expecting at least two pose samples.\n" );
}

Quat SmoothSLERPPoseInterpolation::operator()(double t) const {

  // Make sure that t lies within the range [t0, t0+dt*length(points)]
  VW_ASSERT( t >= m_t0 && t <= m_tend,
	     ArgumentErr() << "Cannot extrapolate point for time "
	     << t << ". Out of valid range. Expecting: "
	     << m_t0 << " <= " << t << " <= " << m_tend << "\n");

  std::vector<double> wts;
  std::vector<int> indices;
  vw::camera::get_wts_and_indices(m_t0, m_dt, t, m_pose_samples.size(),
				  m_num_wts,  m_sigma,
				  wts, indices);

  std::vector<Quat> Q;
  for (size_t i = 0; i < indices.size(); i++)
    Q.push_back(m_pose_samples[indices[i]]);

  return vw::math::slerp_n(wts, Q, 0);
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
  // Loop until next-to-last entry
  for ( size_t i = 0; i < tlc.size() - 1; i++ ) {
    const double this_line = tlc[i].first;
    const double t         = time_offset + tlc[i].second; // The time for this entry
    
    // Compute instantaneous slope at this time = (time diff) / (line diff)
    m_m[this_line] = ( tlc[i].second - tlc[i+1].second ) / ( tlc[i].first - tlc[i+1].first );
    // ?
    m_b[this_line] = t - m_m[this_line] * this_line;
  }
}

double TLCTimeInterpolation::operator()( double line ) const {
  map_type::const_iterator m = m_m.lower_bound( line );
  map_type::const_iterator b = m_b.lower_bound( line );
  if ( m != m_m.begin() ) {
    m--; b--;
  }
  // ?
  return line  * m->second + b->second;
}
