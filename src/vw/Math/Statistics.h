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


/// \file Statistics.h
///
/// Assorted useful statistical routines and functors.
///
#ifndef __MATH_STATISTICS_H__
#define __MATH_STATISTICS_H__

#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>

#include <vw/Core/CompoundTypes.h>
#include <vw/Core/TypeDeduction.h>
#include <vw/Core/Functors.h>
#include <vw/Core/Exception.h>

namespace vw {
namespace math {

//---------------------------------------------------------------------------------

/// Compute the mean of an std container.
/// - It us up to the user to verify values.size() > 0
template <typename T>
double mean(T const& values);

/// Compute the standard deviation of an std container given the mean.
template <typename T>
double standard_deviation(T const& values, double mean_value);

/// Find the median of an std container
template <typename T>
double median(T const& values);



//---------------------------------------------------------------------------------


/// Using the quartile range, determine the values b and e so that
/// all elements in p outsize of [b, e] are outliers.
/* If Q1 and Q3 are the percentiles for pct and 1-pct, the outlier brackets are:
   b = Q1 - outlier_factor*(Q3-Q1)
   e = Q3 + outlier_factor*(Q3-Q1)
   This algorithm works best if the data is distributed rather uniformly.

   Suggested values for the inputs:
   pct = 0.25, outlier_factor: 1.5.
   We expect 0 < pct < 0.5.
*/
template<class T>
void find_outlier_brackets(std::vector<T> const& p,
                           double pct, double outlier_factor,
                           double & b, double & e);


/// CDF (Cumulative Distribution Function) Accumulator
/// Actually it's an approximation. It allows for a more memory efficient
/// calculation of any quantile. Probably most importantly the median.
/// - Use the quantile() function buried way down below to obtain percentile
///   values of the image, useful for intensity stretching of images.
/// - Warning: This class may not function properly for very small numbers of inputs!
///
/// Taken from Numerical Recipes (3rd E) pg 435
template <class ValT>
class CDFAccumulator : public ReturnFixedType<void> {

public:
  CDFAccumulator( size_t buffersize = 1000, size_t quantiles = 251) {
    this->resize( buffersize, quantiles );
  }

  /// Allow user to change post constructor (see ChannelAccumulator)
  void resize( size_t buffersize, size_t quantiles );

  /// Merge in Bundles
  void update();

  /// User update function. (Bundles Data)
  void operator()( ValT const& arg );

  /// Function to merge to CDFs
  void operator()( CDFAccumulator<ValT>& other );

  /// Make this object an exact copy of the other object
  void duplicate(CDFAccumulator<ValT> const& other);
  
  /// Extract a percentile
  ValT quantile( double const& arg ) const;

  // Predefine functions
  ValT median        () const { return quantile(0.5 ); }
  ValT first_quartile() const { return quantile(0.25); }
  ValT third_quartile() const { return quantile(0.75); }

  ValT approximate_mean  ( float const& stepping = 0.1 ) const;
  ValT approximate_stddev( float const& stepping = 0.1 ) const;
  
  
private: // Variables
  size_t m_num_quantiles, m_buffer_idx;
  size_t m_num_samples; // nq, nd, nt
  std::vector<double> m_cdf, m_sample_buf, m_quantile;
  double m_q0, m_qm;  // quantile min and max;

private: // Functions
  // Trapezoidal Rule functor for numeric integration
  template <class InputIterator1, class InputIterator2>
  double trapezoidal_rule( InputIterator1 first1, InputIterator1 last1,
                           InputIterator2 first2 );

  // PDF differentiation. Requires an output vector that is qual
  // in length to the CDF but doesn't touch index 0 or the last element.
  void pdf_differentiation( std::vector<double> const& quantile,
                            std::vector<double> const& cdf,
                            std::vector<double> & output_pdf );

  void correct_pdf_integral_to_1( std::vector<double> const& quantile,
                                  std::vector<double> & pdf,
                                  double additional_scalar );
}; // End class CDFAccumulator



/// Simple histogram management class
class Histogram {

public: // Functions

  /// Default constructor
  Histogram() : m_num_bins(0) {}

  /// Constructor
  Histogram(size_t num_bins, double min_value, double max_value) {
    initialize(num_bins, min_value, max_value);
  }

  // Set up the histogram so values can be added
  void initialize(size_t num_bins, double min_value, double max_value);

  double get_bin_width ()        const { return m_bin_width;        }
  double get_bin_center(int bin) const { return m_bin_centers[bin]; }
  double get_bin_value (int bin) const { return m_bin_values [bin]; }
  size_t get_total_num_values()  const { return m_num_values;       }

  /// Add a single input value to the histogram.
  /// - If saturate is not set, values outside the specified input range are ignored.
  void add_value(double value, bool saturate=true);
  void operator()(double value) { add_value(value); }

  /// Add a value with no bounds checking on the input!
  void add_value_no_check(double value);

  /// Return the bin index containing the specified histogram percentile
  size_t get_percentile(double percentile) const;

  /// Write a simple text file to disk containing the histogram
  void write_to_disk(std::string const& path) const;

private: // Variables

  int    m_num_bins, m_max_bin;
  size_t m_num_values;
  double m_min_value, m_max_value, m_range, m_bin_width;
  std::vector<double> m_bin_centers;
  std::vector<double> m_bin_values;

}; // End class histogram







//---------------------------------------------------------------------------
// The two functors in this class should not be confused with the similar
// classes in Functors.h.  These two are only used in the Geometry.h file.

/// Finds the mean of a set of points.
class MeanFunctor {
  bool m_homogeneous;
public:
  /// If this functor is going to be applied to points in a
  /// projective space (i.e. homogeneous coordinates), you should
  /// set this flag to to true. The resulting mean will be in
  /// the same coordinates.
  MeanFunctor(bool homogeneous_points = false)
    : m_homogeneous(homogeneous_points) {}

  /// This function can use points in any container that supports
  /// the size() and operator[] methods.  The container is usually a
  /// vw::Vector<>, but you could substitute other classes here as well.
  template <class ContainerT>
  ContainerT operator() (std::vector<ContainerT> const& points) const;
};

/// Finds the standard deviation of a set of points.
class StandardDeviationFunctor {
  bool m_homogeneous;
public:
  /// If this functor is going to be applied to points in a
  /// projective space (i.e. homogeneous coordinates), you should
  /// set this flag to to true. The resulting standard deviation
  /// will be in the same coordinates.
  StandardDeviationFunctor(bool homogeneous_points = false)
    : m_homogeneous(homogeneous_points) {}

  /// This function can use points in any container that supports
  /// the size() and operator[] methods.  The container is usually a
  /// vw::Vector<>, but you could substitute other classes here as well.
  template <class ContainerT>
  ContainerT operator() (std::vector<ContainerT> const& points) const;
};


#include <vw/Math/Statistics.tcc>

}} // namespace vw::math

#endif // __MATH_STATISTICS_H__
