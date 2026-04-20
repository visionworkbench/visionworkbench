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


/// \file Statistics.cc
///
/// Non-template Histogram methods. The rest of Statistics is templated
/// and lives in Statistics.tcc / Statistics.h.

#include <vw/Math/Statistics.h>
#include <vw/Core/Exception.h>

#include <cmath>
#include <fstream>
#include <string>

namespace vw {
namespace math {

void Histogram::initialize(size_t num_bins, double min_value, double max_value) {

  if (num_bins == 0)
    vw_throw(ArgumentErr() << "Can't create a Histogram with zero bins!");

  m_num_bins    = num_bins;
  m_max_bin     = num_bins-1;
  m_num_values  = 0;
  m_min_value   = min_value;
  m_max_value   = max_value;
  m_bin_centers.resize(num_bins);
  m_bin_values.resize (num_bins);

  // Compute the bin centers
  m_range     = max_value - min_value;
  m_bin_width = m_range / num_bins;

  for (size_t i=0; i<num_bins; ++i) {
    m_bin_centers[i] = min_value + i*m_bin_width + m_bin_width/2.0;
    m_bin_values [i] = 0;
  }
}

void Histogram::add_value(double value, bool saturate) {
  // Compute the bin
  int bin = static_cast<int>(round(m_max_bin * ((value - m_min_value)/m_range)));
  // Handle out of range input values
  if (bin < 0) {
    if (saturate)
      bin = 0;
    else
      return;
  }
  if (bin >= m_num_bins) {
    if (saturate)
      bin = m_max_bin;
    else
      return;
  }
  m_bin_values[bin] += 1.0;

  ++m_num_values;
}

void Histogram::add_value_no_check(double value) {
  // Compute the bin
  int bin = static_cast<int>(round(m_max_bin * ((value - m_min_value)/m_range)));
  m_bin_values[bin] += 1.0;
  ++m_num_values;
}

size_t Histogram::get_percentile(double percentile) const {

  // Verify the input percentile is in the legal range
  if ((percentile < 0) || (percentile > 1.0)) {
    vw_throw(ArgumentErr() << "get_histogram_percentile: illegal percentile request: "
             << percentile << "\n");
  }

  double sum = static_cast<double>(m_num_values);

  double running_percentile = 0;
  for (int i=0; i<m_num_bins; ++i) {
    double this_percent = get_bin_value(i) / sum;
    running_percentile += this_percent;
    if (running_percentile >= percentile)
      return i;
  }
  vw_throw(LogicErr() << "get_histogram_percentile: Illegal histogram encountered!");
  return 0;
}

void Histogram::write_to_disk(std::string const& path) const {
  std::ofstream f(path.c_str());
  f << "# Histogram_Value, Histogram_Count";
  for (int i=0; i<m_num_bins; ++i) {
    f << get_bin_center(i) << ", " << get_bin_value(i) << "\n";
  }
  f.close();
}

}} // namespace vw::math
