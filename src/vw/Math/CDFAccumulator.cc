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


/// \file CDFAccumulator.cc

#include <vw/Math/CDFAccumulator.h>
#include <vw/Core/Exception.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace vw {
namespace math {

template <class ValT>
void CDFAccumulator<ValT>::resize(size_t buffersize, size_t quantiles) {
  VW_ASSERT(quantiles > 0, LogicErr() << "Cannot have 0 quantiles");
  m_buffer_idx = m_num_samples = 0;
  m_sample_buf.resize(buffersize);

  m_q0 =  std::numeric_limits<double>::max();
  m_qm = -std::numeric_limits<double>::max();

  m_num_quantiles = quantiles;
  if (!(quantiles%2))
    m_num_quantiles++;
  m_quantile.resize(m_num_quantiles);
  m_cdf.resize(m_num_quantiles);

  // Setting a generic cdf to start things off, where 80% of the
  // distribution is in the middle third.
  size_t third = m_num_quantiles/3;
  size_t third2 = third*2;
  double slope = 10.0 / double(third);
  double first_tertile_gain = 1.0 - slope;

  // Filling middle
  for (size_t j = third; j <= third2; j++)
    m_cdf[j] = 0.8*(double(j-third)/double(third2-third))+0.1;
  // Filling first tertile
  for (ssize_t j = third-1; j >= 0; j--)
    m_cdf[j] = first_tertile_gain*m_cdf[j+1];
  // Filling third tertile
  for (size_t j = third2+1; j < m_num_quantiles; j++)
    m_cdf[j] = 1.0 - first_tertile_gain*(1.0-m_cdf[j-1]);
}

template <class ValT>
void CDFAccumulator<ValT>::update() {
  // Early exit if an update already happened
  if (!m_buffer_idx)
    return;

  size_t jd=0, jq=1;
  double target, told=0, tnew=0, qold, qnew;
  std::vector<double> m_new_quantile(m_num_quantiles);
  std::sort(m_sample_buf.begin(),
            m_sample_buf.begin()+m_buffer_idx); // For partial updates
  // Setting to global min and max;
  qold = qnew = m_quantile[0] = m_new_quantile[0] = m_q0;
  m_quantile.back() = m_new_quantile.back() = m_qm;
  // .. then setting comparable probabilities
  m_cdf[0] = std::min(0.5/(m_buffer_idx+m_num_samples),
                      0.5*m_cdf[1]);
  m_cdf.back() = std::max(1-0.5/(m_buffer_idx+m_num_samples),
                          0.5*(1+m_cdf[m_num_quantiles-2]));
  // Looping over target probability values for interpolation
  for (size_t iq = 1; iq < m_num_quantiles-1; iq++) {
    target = (m_num_samples+m_buffer_idx)*m_cdf[iq];
    if (tnew < target) {
      while (1) {
        // Locating a succession of abscissa-ordinate pairs
        // (qnew,tnew) that are the discontinuities of value or
        // slope, breaking to perform an interpolation as we cross
        // each target.
        if (jq < m_num_quantiles &&
            (jd >= m_buffer_idx ||
             m_quantile[jq] < m_sample_buf[jd])) {
          // Found slope discontinuity from old CDF.
          qnew = m_quantile[jq];
          tnew = jd + m_num_samples*m_cdf[jq++];
          if (tnew >= target) break;
        } else {
          // Found value discontinuity from batch data CDF.
          qnew = m_sample_buf[jd];
          tnew = told;
          if (m_quantile[jq] > m_quantile[jq-1])
            tnew += m_num_samples*(m_cdf[jq]-m_cdf[jq-1])*
              (qnew-qold)/(m_quantile[jq]-m_quantile[jq-1]);
          jd++;
          if (tnew >= target) break;
          told = tnew++;
          qold = qnew;
          if (tnew >= target) break;
        }
        told = tnew;
        qold = qnew;
      }
    }
    // Performing new interpolation
    if (tnew == told)
      m_new_quantile[iq] = 0.5*(qold+qnew);
    else
      m_new_quantile[iq] = qold + (qnew-qold)*(target-told)/(tnew-told);
    told = tnew;
    qold = qnew;
  }
  // Reset'n
  m_quantile = m_new_quantile;
  m_num_samples += m_buffer_idx;
  m_buffer_idx = 0;
}

template <class ValT>
void CDFAccumulator<ValT>::operator()(ValT const& arg) {
  // Assimilate, We are the Borg, your data is my data!
  m_sample_buf[m_buffer_idx++] = arg;
  if (arg < m_q0)
    m_q0 = arg; // stretch cdf?
  if (arg > m_qm)
    m_qm = arg;
  if (m_buffer_idx == m_sample_buf.size())
    update(); // merge cdf?
}

template <class ValT>
ValT CDFAccumulator<ValT>::quantile(double const& arg) const {
  // Defend against querying before any data has been fed. Previously callers
  // (e.g. DisparityCdfFunctor) had to guard externally; now the class does it.
  if (m_num_samples + m_buffer_idx == 0)
    return ValT(0);

  double q;
  size_t jl=0, jh=m_num_quantiles-1, j;
  while (jh - jl > 1) {
    j = (jh+jl)>>1;
    if (arg > m_cdf[j]) jl=j;
    else jh=j;
  }
  j = jl;
  q = m_quantile[j]+(m_quantile[j+1]-m_quantile[j])*(arg-m_cdf[j])/(m_cdf[j+1]-m_cdf[j]);

  // Keeping estimate in CDF
  return std::max(m_quantile[0], std::min(m_quantile.back(), q));
}

template <class ValT>
ValT CDFAccumulator<ValT>::approximate_mean(float const& stepping) const {
  ValT   mean = 0;
  size_t count = 0;
  for (float i = stepping; i < 1+stepping; i+=stepping) {
    count++;
    mean += (quantile(i) + quantile(i-stepping)) / 2.0;
  }
  return mean / ValT(count);
}

template <class ValT>
ValT CDFAccumulator<ValT>::approximate_stddev(float const& stepping) const {
  ValT mean = approximate_mean(stepping);
  ValT stddev = 0;
  size_t count = 0;
  for (float i = stepping; i < 1+stepping; i+=stepping) {
    count++;
    stddev += pow((quantile(i) + quantile(i-stepping))/2-mean, 2);
  }
  return sqrt(stddev/ValT(count));
}

// Explicit instantiations. Add more here if a new ValT is needed.
template class CDFAccumulator<float>;
template class CDFAccumulator<double>;

}} // namespace vw::math
