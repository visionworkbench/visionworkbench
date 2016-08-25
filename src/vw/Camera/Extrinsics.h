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


/// \file Extrinsics.h
///
/// Utilities for describing extrinsic camera parameters (position and pose).


// TODO: Most of these are generic interpolation algorithms, they ought to 
//       live somewhere else.

#ifndef __VW_CAMERA_EXTRINSICS_H__
#define __VW_CAMERA_EXTRINSICS_H__

#include <vw/Math/Quaternion.h>
#include <vw/Math/LinearAlgebra.h>

#include <iostream>
#include <map>
#include <vector>

namespace vw {
namespace camera {

  // --------------------------------------------------------------------------
  // POSITION INTERPOLATION - can also be used for velocity.
  // --------------------------------------------------------------------------

  /// Simple position interpolation when only an initial point and velocity are known.
  class LinearPositionInterpolation {
    Vector3 m_initial_position, m_initial_velocity;

  public:
    LinearPositionInterpolation(Vector3 const& initial_position,
                        				Vector3 const& initial_velocity);

    Vector3 operator()(double t) const;
  };

  /// Simple linear interpolation between a series of positions.
  class LinearPiecewisePositionInterpolation {
    std::vector<Vector3> m_position_samples;
    double m_t0, m_dt, m_tend;
  public:
    LinearPiecewisePositionInterpolation(){}
    LinearPiecewisePositionInterpolation( std::vector<Vector3> const& position_samples,
					  double t0, double dt );

    Vector3 operator()( double t ) const;
    double get_t0  () const { return m_t0;  }
    double get_dt  () const { return m_dt;  }
    double get_tend() const { return m_tend;}
  };

  /// Performs smooth interpolation between sparse pose data points using gaussians.
  class SmoothPiecewisePositionInterpolation {
    std::vector<Vector3> m_position_samples;
    double m_t0, m_dt, m_tend;
    int m_num_wts;
    double m_sigma;

  public:
    SmoothPiecewisePositionInterpolation(){}
    SmoothPiecewisePositionInterpolation(std::vector<Vector3> const& pose_samples,
					 double t0, double dt,
					 int num_wts, double sigma);

    /// Compute the pose at a given time t.  The pose will be an interpolated value
    Vector3 operator()(double t) const;

    double get_t0  () const { return m_t0;  }
    double get_dt  () const { return m_dt;  }
    double get_tend() const { return m_tend;}
    std::vector<int> get_indices_of_largest_weights(double t) const;
  };
  
  /// Performs Lagrangian interpolation between data points with flexible times.
  /// - If it becomes useful this can be combined with the fixed interval Lagrangian
  ///   class using a boolean check.
  class LagrangianInterpolationVarTime {
    std::vector<Vector3> m_samples;
    std::vector<double > m_times;
    int m_radius;
  public:
    /// Construct with a set of data samples and times.
    /// - The radius is the number of points before and after time t used for interpolation.
    LagrangianInterpolationVarTime(std::vector<Vector3> const& samples, 
                                   std::vector<double > const& times, int radius=4);
  
    /// Compute the interpolated value at a given time t.
    Vector3 operator()(double t) const;
  };

  /// Performs Lagrangian interpolation between data points with constant time intervals.
  class LagrangianInterpolation {
    std::vector<Vector3> m_samples;
    double m_start_time, m_time_delta, m_last_time;
    int m_radius;
    std::vector<double> m_denoms;
    mutable std::vector<double> m_times_temp;
  public:
    /// Construct with a set of data samples and times.
    /// - The radius is the number of points before and after time t used for interpolation.
    LagrangianInterpolation(std::vector<Vector3> const& samples, 
                            double start_time, double time_delta, double last_time,
                            int radius=4);
  
    /// Compute the interpolated value at a given time t.
    Vector3 operator()(double t) const;
  };


  /// Interpolation between a series of positions incorporating accelleration information.
  class PiecewiseAPositionInterpolation {
    std::vector<Vector3> m_position_samples, m_velocity;
    double m_t0, m_dt, m_tend;
  public:
    PiecewiseAPositionInterpolation( std::vector<Vector3> const& position_samples,
				     std::vector<Vector3> const& velocity_samples,
				     double t0, double dt );

    Vector3 operator()( double t ) const;
    double get_t0  () const { return m_t0;  }
    double get_dt  () const { return m_dt;  }
    double get_tend() const { return m_tend;}
  };


  class Curve3DPositionInterpolation {
    Matrix<double,3,3> m_cached_fit;

  public:
    /// The constructor fits a 3D curve to the provided points and
    /// caches the fit so that it can be quickly evaluated using the
    /// call operator.
    ///
    /// This routine will fit a 3d parametric curve (2nd degree
    /// polynomial) to a series of N points in the Nx3 matrix points.
    ///
    /// The curve is specified by the coefficient
    /// matrix, which has the form:
    ///
    ///          | a b c |
    ///   A  =   | d e f |
    ///          | g h i |
    ///
    /// and the curve is is evaluatied using the equation:
    ///
    ///              |  1  |
    ///   x(t) = A * |  t  |
    ///              | t^2 |
    ///
    Curve3DPositionInterpolation(std::vector<Vector3> const& position_samples,
				 double t0, double dt);

    /// The call operator evaluates the curve at the given time t.
    /// See the constructor documentation for a summary of the
    /// mathematical operations that are carried out in this function.
    Vector3 operator()(double t) const;
  };

  /// Cubic Hermite Spline interpolation or CSpline. Assumes you
  /// provide the velocity measurements.
  class HermitePositionInterpolation {
    std::vector<Vector3> m_position_samples, m_velocity;
    double m_t0, m_dt, m_tend;
  public:
    HermitePositionInterpolation( std::vector<Vector3> const& position_samples,
				  std::vector<Vector3> const& velocity_samples,
				  double t0, double dt );

    Vector3 operator()( double t ) const;
    double get_t0  () const { return m_t0;  }
    double get_dt  () const { return m_dt;  }
    double get_tend() const { return m_tend;}
  };

  // --------------------------------------------------------------------------
  // POSE INTERPOLATION
  // --------------------------------------------------------------------------

  /// Always returns the input pose.
  class ConstantPoseInterpolation {
    Quat m_pose;
  public:
    ConstantPoseInterpolation(Quat const& pose);

    inline Quat operator()(double /*t*/) const {
      return m_pose;
    }
  };

  /// Performs interpolation between sparse pose data points using the
  /// spherical linear interpolation algorithm.
  class SLERPPoseInterpolation {
    std::vector<Quat> m_pose_samples;
    double m_t0, m_dt, m_tend;

  public:
    SLERPPoseInterpolation(){}
    SLERPPoseInterpolation(std::vector<Quat> const& pose_samples, double t0, double dt);

    /// Compute the pose at a given time t.  The pose will be an interpolated value
    Quat operator()(double t) const;

    double get_t0  () const { return m_t0;  }
    double get_dt  () const { return m_dt;  }
    double get_tend() const { return m_tend;}
  };

  /// Simple slerp interpolation between a table of pointing directions arranged on a grid.
  class SlerpGridPointingInterpolation {
    std::vector< std::vector<vw::Vector3> > m_directions;
    double m_row0, m_drow, m_row_end;
    double m_col0, m_dcol, m_col_end;
  public:
    SlerpGridPointingInterpolation(){}
    SlerpGridPointingInterpolation(std::vector< std::vector<vw::Vector3> > const& directions,
                                   double row0, double drow, double col0, double dcol);

    // Careful here, pix[0] is a column, and pix[1] is a row,
    // so we'll access directions(pix[1], pix[0]).
    Vector3 operator()( vw::Vector2 const& pix ) const;
    double get_row0    () const { return m_row0;    }
    double get_drow    () const { return m_drow;    }
    double get_row_end () const { return m_row_end; }
    double get_col0    () const { return m_col0;    }
    double get_dcol    () const { return m_dcol;    }
    double get_col_end () const { return m_col_end; }
  };

  /// Performs smooth interpolation between sparse pose data points using the
  /// spherical linear interpolation algorithm. Highly experimental and not tested.
  class SmoothSLERPPoseInterpolation {
    std::vector<Quat> m_pose_samples;
    double m_t0, m_dt, m_tend;
    int m_num_wts;
    double m_sigma;

  public:
    SmoothSLERPPoseInterpolation(){}
    SmoothSLERPPoseInterpolation(std::vector<Quat> const& pose_samples, double t0, double dt,
				 int num_wts, double sigma);

    /// Compute the pose at a given time t.  The pose will be an interpolated value
    Quat operator()(double t) const;

    double get_t0  () const { return m_t0;  }
    double get_dt  () const { return m_dt;  }
    double get_tend() const { return m_tend;}
  };

  // --------------------------------------------------------------------------
  // TIME INTERPOLATION
  // --------------------------------------------------------------------------

  /// Simple linear interpolation of the time at a given line.
  class LinearTimeInterpolation {
    double m_t0, m_dt;

  public:
    LinearTimeInterpolation( double initial_time, double time_per_line );

    double operator()( double line ) const;
  };

  ///
  class TLCTimeInterpolation {
    typedef std::map<double, double> map_type;
    map_type m_m, m_b; // Tables keyed on line: time = m * line + b;

  public:
    // TLC is straight from the IMG XML tag from Digital Globe
    // products. The pairings are expected to be <Line, Time>.
    TLCTimeInterpolation( std::vector<std::pair<double, double> > const& tlc, double time_offset = 0 );

    double operator()( double line ) const;
  };

  
}} // namespace vw::camera

#endif // __VW_CAMERA_EXTRINSICS_H__
