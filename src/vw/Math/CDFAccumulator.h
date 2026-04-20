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


/// \file CDFAccumulator.h
///
/// Streaming CDF approximation for memory-efficient quantile estimation.
/// Instantiated for float and double in CDFAccumulator.cc.

#ifndef __VW_MATH_CDFACCUMULATOR_H__
#define __VW_MATH_CDFACCUMULATOR_H__

#include <vw/Core/Functors.h>

#include <vector>

namespace vw {
namespace math {

/// CDF (Cumulative Distribution Function) Accumulator
/// Actually it's an approximation. It allows for a more memory efficient
/// calculation of any quantile. Probably most importantly the median.
/// - Use the quantile() function buried way down below to obtain percentile
///   values of the image, useful for intensity stretching of images.
/// - Warning: This class may not function properly for very small numbers of inputs!
///
/// Taken from Numerical Recipes (3rd E) pg 435
template <class ValT>
class CDFAccumulator: public ReturnFixedType<void> {

public:
  CDFAccumulator(size_t buffersize = 1000, size_t quantiles = 251) {
    this->resize(buffersize, quantiles);
  }

  /// Allow user to change post constructor (see ChannelAccumulator)
  void resize(size_t buffersize, size_t quantiles);

  /// Merge in Bundles
  void update();

  /// User update function. (Bundles Data)
  void operator()(ValT const& arg);

  /// Function to merge to CDFs
  void operator()(CDFAccumulator<ValT>& other);

  /// Make this object an exact copy of the other object
  void duplicate(CDFAccumulator<ValT> const& other);

  /// Extract a percentile. Returns ValT(0) if no samples have been fed yet.
  ValT quantile(double const& arg) const;

  // Predefine functions
  ValT median        () const { return quantile(0.5 ); }
  ValT first_quartile() const { return quantile(0.25); }
  ValT third_quartile() const { return quantile(0.75); }

  ValT approximate_mean  (float const& stepping = 0.1) const;
  ValT approximate_stddev(float const& stepping = 0.1) const;


private: // Variables
  size_t m_num_quantiles, m_buffer_idx;
  size_t m_num_samples; // nq, nd, nt
  std::vector<double> m_cdf, m_sample_buf, m_quantile;
  double m_q0, m_qm;  // quantile min and max;

private: // Functions
  // Trapezoidal Rule functor for numeric integration
  template <class InputIterator1, class InputIterator2>
  double trapezoidal_rule(InputIterator1 first1, InputIterator1 last1,
                          InputIterator2 first2);

  // PDF differentiation. Requires an output vector that is qual
  // in length to the CDF but doesn't touch index 0 or the last element.
  void pdf_differentiation(std::vector<double> const& quantile,
                           std::vector<double> const& cdf,
                           std::vector<double> & output_pdf);

  void correct_pdf_integral_to_1(std::vector<double> const& quantile,
                                 std::vector<double> & pdf,
                                 double additional_scalar);
}; // End class CDFAccumulator

}} // namespace vw::math

#endif // __VW_MATH_CDFACCUMULATOR_H__
