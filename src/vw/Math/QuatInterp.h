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

/// \file QuatInterp.h
///
/// Quaternion interpolation algorithms (SLERP, Lagrangian, etc.).

#ifndef __VW_MATH_QUAT_INTERP_H__
#define __VW_MATH_QUAT_INTERP_H__

#include <vw/Math/Quaternion.h>
#include <vw/Math/Vector.h>

#include <boost/shared_ptr.hpp>

#include <vector>

namespace vw {
  namespace math {
    template <class Point>
    class catmull_rom;
  }
}

namespace vw {

  /// Performs Lagrangian interpolation between quaternion data points with
  /// flexible times. Same algorithm as LagrangianInterpolationVarTime but
  /// operates on Quat (4 components) and normalizes the result.
  class QuatLagrangianInterpolationVarTime {
  public:
    std::vector<Quat>   m_samples;
    std::vector<double> m_times;
    int m_radius; // interpolation order = 2 * m_radius

    /// Construct with a set of quaternion samples and times.
    /// - The radius is the number of points before and after time t
    ///   used for interpolation.
    QuatLagrangianInterpolationVarTime(std::vector<Quat> const& samples,
                                       std::vector<double> const& times,
                                       int radius = 4);

    /// Compute the interpolated quaternion at a given time t.
    /// The result is normalized.
    Quat operator()(double t) const;
  };

  /// Performs Lagrangian interpolation between quaternions with
  /// constant time intervals. Given the input radius, form the
  /// polynomial of degree 2*radius - 1 which goes through the
  /// 2*radius samples whose times are closest to the time at which
  /// interpolation must happen.
  /// TODO(oalexan1): What if is desired to use an odd number of
  /// samples?

  // TODO(oalexan1): This needs to be tested and compared with SLERPPoseInterpolation,
  // with and without use of splines in the latter.
  class QuatLagrangianInterpolation {
  public:
    std::vector<Quat> m_samples;
    double m_start_time, m_time_delta, m_last_time;
    int m_radius;
    std::vector<double> m_denoms;
    mutable std::vector<double> m_times_temp;
    /// Construct with a set of data samples and times.
    /// - The radius is the number of points before and after time t used for interpolation.
    QuatLagrangianInterpolation(std::vector<Quat> const& samples,
                            double start_time, double time_delta, double last_time,
                            int radius=4);

    double get_t0  () const { return m_start_time; }
    double get_dt  () const { return m_time_delta; }
    double get_tend() const { return m_last_time;  }

    /// Compute the interpolated value at a given time t.
    Quat operator()(double t) const;
  };

  /// Always returns the input pose.
  class ConstantPoseInterpolation {
  public:
    Quat m_pose;
    ConstantPoseInterpolation(Quat const& pose);

    inline Quat operator()(double /*t*/) const {
      return m_pose;
    }
  };

  /// Performs interpolation between sparse pose data points using the
  /// spherical linear interpolation algorithm. Using splines is not
  /// enabled by default and is experimental.
  class SLERPPoseInterpolation {
public:
    SLERPPoseInterpolation(){}
    SLERPPoseInterpolation(std::vector<Quat> const& pose_samples, double t0, double dt,
                           bool use_splines = false);

    /// Compute the pose at a given time t.  The pose will be an interpolated value
    Quat operator()(double t) const;

    double get_t0  () const { return m_t0;  }
    double get_dt  () const { return m_dt;  }
    double get_tend() const { return m_tend;}
    std::vector<Quat> m_pose_samples;
    double m_t0, m_dt, m_tend;
    bool m_use_splines;
    boost::shared_ptr<vw::math::catmull_rom<std::array<double, 4>>> m_spline_ptr;
  };

  /// Simple slerp interpolation between a table of pointing directions arranged on a grid.
  class SlerpGridPointingInterpolation {
  public:
    std::vector<std::vector<vw::Vector3>> m_directions;
    double m_row0, m_drow, m_row_end;
    double m_col0, m_dcol, m_col_end;
    SlerpGridPointingInterpolation(){}
    SlerpGridPointingInterpolation(std::vector< std::vector<vw::Vector3>> const& directions,
                                   double row0, double drow, double col0, double dcol);

    // Careful here, pix[0] is a column, and pix[1] is a row,
    // so we'll access directions(pix[1], pix[0]).
    Vector3 operator()(vw::Vector2 const& pix) const;
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
  public:
    std::vector<Quat> m_pose_samples;
    double m_t0, m_dt, m_tend;
    int m_num_wts;
    double m_sigma;

    SmoothSLERPPoseInterpolation(){}
    SmoothSLERPPoseInterpolation(std::vector<Quat> const& pose_samples, double t0, double dt,
                                 int num_wts, double sigma);

    /// Compute the pose at a given time t.  The pose will be an interpolated value
    Quat operator()(double t) const;

    double get_t0  () const { return m_t0;  }
    double get_dt  () const { return m_dt;  }
    double get_tend() const { return m_tend;}
  };

} // namespace vw

#endif // __VW_MATH_QUAT_INTERP_H__
