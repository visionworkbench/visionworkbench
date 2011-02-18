// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

namespace vw {
namespace camera {

  // --------------------------------------------------------------------------
  // POSITION INTERPOLATION
  // --------------------------------------------------------------------------

  class LinearPositionInterpolation {
    Vector3 m_initial_position, m_velocity_vector;

  public:
    LinearPositionInterpolation(Vector3 initial_position, Vector3 velocity_vector) :
      m_initial_position(initial_position), m_velocity_vector(velocity_vector) {}

    Vector3 operator()(double t) const {
      return m_initial_position + t * m_velocity_vector;
    }

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
    Vector3 operator()(double t) const {
      Vector3 T(1, t, t*t);
      return m_cached_fit * T;
    }
  };


  // --------------------------------------------------------------------------
  // POSE INTERPOLATION
  // --------------------------------------------------------------------------

  class ConstantPoseInterpolation {
    Quat m_pose;
  public:
    ConstantPoseInterpolation(Quat const& pose) : m_pose(pose) {}

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
    SLERPPoseInterpolation(std::vector<Quat > const& pose_samples, double t0, double dt) :
      m_pose_samples(pose_samples), m_t0(t0), m_dt(dt) {}

    /// Compute the pose at a given time t.  The pose will be an interpolated value
    Quat operator()(double t) const;

  };

}} // namespace vw::camera

#endif // __VW_CAMERA_EXTRINSICS_H__
