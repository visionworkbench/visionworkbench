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

/// \file PositionInterp.h
///
/// Position interpolation algorithms (linear, piecewise, Lagrangian,
/// Hermite, etc.). Can also be used for velocity.

#ifndef __VW_MATH_POSITION_INTERP_H__
#define __VW_MATH_POSITION_INTERP_H__

#include <vw/Math/Vector.h>
#include <vw/Math/LinearAlgebra.h>

#include <vector>

namespace vw {

  // Given the values ti = t0 + i*dt for i = 0, 1, ..., num_samples - 1,
  // and the value t, compute the weights
  // w_i=exp(-sigma*(t-ti)^2/dt^2). Keep the largest num_wts weights.
  // Normalize the weights to add to 1 and return them. Used for
  // Gaussian-weighted interpolation. This works even when t is out of range.
  void get_wts_and_indices(double t0, double dt, double t, int num_samples,
                           int num_wts, double sigma,
                           std::vector<double> & wts, std::vector<int> & indices);

  /// Simple position interpolation when only an initial point and velocity are known.
  class LinearPositionInterpolation {
  public:
    Vector3 m_initial_position, m_initial_velocity;

    LinearPositionInterpolation(Vector3 const& initial_position,
                                Vector3 const& initial_velocity);

    Vector3 operator()(double t) const;
  };

  /// Simple linear interpolation between a series of positions.
  class LinearPiecewisePositionInterpolation {
  public:
    std::vector<Vector3> m_position_samples;
    double m_t0, m_dt, m_tend;
    LinearPiecewisePositionInterpolation(){}
    LinearPiecewisePositionInterpolation(std::vector<Vector3> const& position_samples,
                                         double t0, double dt);

    Vector3 operator()(double t) const;
    double get_t0  () const { return m_t0;  }
    double get_dt  () const { return m_dt;  }
    double get_tend() const { return m_tend;}
  };

  /// Performs smooth interpolation between sparse pose data points using gaussians.
  class SmoothPiecewisePositionInterpolation {
  public:
    std::vector<Vector3> m_position_samples;
    double m_t0, m_dt, m_tend;
    int m_num_wts;
    double m_sigma;

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
  /// - See QuatLagrangianInterpolationVarTime for the quaternion version.
  class LagrangianInterpolationVarTime {
  public:
    std::vector<Vector3> m_samples;
    std::vector<double> m_times;
    int m_radius; // interpolation order = 2 * m_radius

    /// Construct with a set of data samples and times.
    /// - The radius is the number of points before and after time t
    ///   used for interpolation.
    LagrangianInterpolationVarTime(std::vector<Vector3> const& samples,
                                   std::vector<double> const& times,
                                   int radius = 4);

    /// Compute the interpolated value at a given time t.
    Vector3 operator()(double t) const;
  };

  /// Performs Lagrangian interpolation between data points with
  /// constant time intervals. Given the input radius, form the
  /// polynomial of degree 2*radius - 1 which goes through the
  /// 2*radius samples whose times are closest to the time at which
  /// interpolation must happen.
  /// TODO(oalexan1): What if is desired to use an odd number of
  /// samples?
  class LagrangianInterpolation {
  public:
    double m_start_time, m_time_delta, m_last_time;
    int m_radius;
    std::vector<double> m_denoms;
    mutable std::vector<double> m_times_temp;
    std::vector<Vector3> m_samples;

    /// Construct with a set of data samples and times.
    /// - The radius is the number of points before and after time t used for interpolation.
    LagrangianInterpolation(std::vector<Vector3> const& samples,
                            double start_time, double time_delta, double last_time,
                            int radius=4);

    double get_t0  () const { return m_start_time; }
    double get_dt  () const { return m_time_delta; }
    double get_tend() const { return m_last_time;  }

    /// Compute the interpolated value at a given time t.
    Vector3 operator()(double t) const;
  };

  /// Interpolation between a series of positions using velocity information.
  class PiecewiseAPositionInterpolation {
  public:
    std::vector<Vector3> m_position_samples, m_velocity;
    double m_t0, m_dt, m_tend;
    PiecewiseAPositionInterpolation(std::vector<Vector3> const& position_samples,
                                    std::vector<Vector3> const& velocity_samples,
                                    double t0, double dt);

    Vector3 operator()(double t) const;
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

  /// Cubic Hermite Spline interpolation or CSpline. Assumes the velocity measurements are provided.
  class HermitePositionInterpolation {
  public:
    std::vector<Vector3> m_position_samples, m_velocity;
    double m_t0, m_dt, m_tend;
    HermitePositionInterpolation(std::vector<Vector3> const& position_samples,
                                 std::vector<Vector3> const& velocity_samples,
                                 double t0, double dt);

    Vector3 operator()(double t) const;
    double get_t0  () const { return m_t0;  }
    double get_dt  () const { return m_dt;  }
    double get_tend() const { return m_tend;}
  };

} // namespace vw

#endif // __VW_MATH_POSITION_INTERP_H__
