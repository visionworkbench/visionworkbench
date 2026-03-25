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

#include <vw/Camera/TimeInterp.h>

using namespace vw;

//======================================================================
// LinearTimeInterpolation class

LinearTimeInterpolation::LinearTimeInterpolation(double initial_time, double time_per_line) :
  m_t0(initial_time), m_dt(time_per_line) {}

double LinearTimeInterpolation::operator()(double line) const {
  return m_dt * line + m_t0;
}

//======================================================================
// TLCTimeInterpolation class

TLCTimeInterpolation::TLCTimeInterpolation(std::vector<std::pair<double, double>> const& tlc,
                                           double time_offset) {

  // These original inputs can be retrieved later if the precomputed data as below
  // are not sufficient.
  m_tlc = tlc;
  m_time_offset = time_offset;

  // Loop until next-to-last entry
  for (size_t i = 0; i + 1 < tlc.size(); i++) {
    const double this_line = tlc[i].first;
    const double t         = time_offset + tlc[i].second; // The time for this entry

    // Compute the slope between this time instance and the next
    m_m[this_line] = (tlc[i+1].second - tlc[i].second) / (tlc[i+1].first - tlc[i].first);

    // Compute the intercept
    m_b[this_line] = t - m_m[this_line] * this_line;
  }
}

double TLCTimeInterpolation::operator()(double line) const {
  map_type::const_iterator m = m_m.lower_bound(line);
  map_type::const_iterator b = m_b.lower_bound(line);
  if (m != m_m.begin()) {
    m--; b--;
  }
  // Find time at given line
  return line  * m->second + b->second;
}
