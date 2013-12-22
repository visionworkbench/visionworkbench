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
/// Utilities for describing extrinsic camera parameters (position and
/// pose).
///
///
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
  // POSITION INTERPOLATION
  // --------------------------------------------------------------------------

  class LinearPositionInterpolation {
    Vector3 m_initial_position, m_velocity_vector;

  public:
    LinearPositionInterpolation(Vector3 const& initial_position,
                                Vector3 const& velocity_vector);

    Vector3 operator()(double t) const;
  };

  class Curve3DPositionInterpolation {
    Matrix<double,3,3> m_cached_fit;

  public:
    /// The constructor fits a 3D curve to the provided points and
    /// caches the fit so that it can be quickly evaluated using the
    /// call operator.
    ///
    /// This routine will fit a 3d parametric curve (2nd degree
    /// polynomial) to a serios of N points in the Nx3 matrix points.
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

  // Cubic Hermite Spline interpolation or CSpline. Assumes you
  // provide the velocity measurements.
  class HermitePositionInterpolation {
    std::vector<Vector3> m_position, m_velocity;
    double m_t0, m_dt;
  public:
    HermitePositionInterpolation( std::vector<Vector3> const& position_samples,
                                  std::vector<Vector3> const& velocity_samples,
                                  double t0, double dt );

    Vector3 operator()( double t ) const;
  };

  class PiecewiseAPositionInterpolation {
    std::vector<Vector3> m_position, m_velocity;
    double m_t0, m_dt;
  public:
    PiecewiseAPositionInterpolation( std::vector<Vector3> const& position_samples,
                                     std::vector<Vector3> const& velocity_samples,
                                     double t0, double dt );

    Vector3 operator()( double t ) const;
  };

  class LinearPiecewisePositionInterpolation {
    std::vector<Vector3> m_position;
    double m_t0, m_dt;
  public:
    LinearPiecewisePositionInterpolation( std::vector<Vector3> const& position_samples,
                                          double t0, double dt );

    Vector3 operator()( double t ) const;
  };


  // --------------------------------------------------------------------------
  // POSE INTERPOLATION
  // --------------------------------------------------------------------------

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
    std::vector<Quat > m_pose_samples;
    double m_t0, m_dt;

    Quat slerp(double alpha, Quat const& a, Quat const& b, int spin) const;

  public:
    SLERPPoseInterpolation(std::vector<Quat > const& pose_samples, double t0, double dt);

    /// Compute the pose at a given time t.  The pose will be an interpolated value
    Quat operator()(double t) const;

  };

  // --------------------------------------------------------------------------
  // POSE INTERPOLATION
  // --------------------------------------------------------------------------

  class LinearTimeInterpolation {
    double m_t0, m_dt;

  public:
    LinearTimeInterpolation( double initial_time, double time_per_line );

    double operator()( double line ) const;
  };

  class TLCTimeInterpolation {
    typedef std::map<double, double> map_type;
    map_type m_m, m_b; // Tables keyed on line: time = m * line + b;

  public:
    // TLC is straight from the IMG XML tag from Digital Globe
    // products. The pairings are expected to be <Line, Time>.
    TLCTimeInterpolation( std::vector<std::pair<double, double> > const& tlc,
                          double time_offset = 0 );

    double operator()( double line ) const;
  };

}} // namespace vw::camera

#endif // __VW_CAMERA_EXTRINSICS_H__
