// __BEGIN_LICENSE__
//  Copyright (c) 2006-2025, United States Government as represented by the
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

/// \file TimeInterp.h
///
/// Time interpolation algorithms for linescan cameras.

#ifndef __VW_CAMERA_TIME_INTERP_H__
#define __VW_CAMERA_TIME_INTERP_H__

#include <map>
#include <utility>
#include <vector>

namespace vw {

  /// Simple linear interpolation of the time at a given line.
  class LinearTimeInterpolation {
  public:
    LinearTimeInterpolation(double initial_time, double time_per_line);
    double operator()(double line) const;
    double m_t0, m_dt;
  };

  // Compute the time at given line using piecewise linear interpolation.
  // TLC is straight from the IMG XML tag from Digital Globe
  // products. The pairings are expected to be <Line, Time>.
  class TLCTimeInterpolation {
  public:

    typedef std::map<double, double> map_type;

    // Tables keyed on line: time = m * line + b;
    map_type m_m, m_b;

    std::vector<std::pair<double, double>> m_tlc;
    double m_time_offset;

    TLCTimeInterpolation(std::vector<std::pair<double, double>> const& tlc,
                         double time_offset = 0.0);

    double operator()(double line) const;
  };

} // namespace vw

#endif // __VW_CAMERA_TIME_INTERP_H__
