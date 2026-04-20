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
/* If Q1 and Q3 are the percentiles for pct_factor and 1-pct_factor, the outlier brackets are:
   b = Q1 - outlier_factor*(Q3-Q1)
   e = Q3 + outlier_factor*(Q3-Q1)
   This algorithm works best if the data is distributed rather uniformly.

   Suggested values for the inputs:
   pct_factor = 0.25, outlier_factor: 1.5.
   We expect 0 < pct_factor < 0.5.
   
   Return true on success.
*/
template<class T>
bool find_outlier_brackets(std::vector<T> const& vec,
                           double pct_factor, double outlier_factor,
                           double & b, double & e);


// CDFAccumulator moved to vw/Math/CDFAccumulator.h.

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


//---------------------------------------------------------------------------------
// Template implementations

/// Compute the mean of an std container.
/// - It us up to the user to verify values.size() > 0
template <typename T>
double mean(T const& values) {

  double result = 0;
  typename T::const_iterator iter;
  for (iter = values.begin(); iter != values.end(); ++iter)
    result += static_cast<double>(*iter);

  return result / static_cast<double>(values.size());
}

/// Compute the standard deviation of an std container given the mean.
template <typename T>
double standard_deviation(T const& values, double mean_value) {

  double result = 0, diff = 0;
  typename T::const_iterator iter;
  for (iter = values.begin(); iter != values.end(); ++iter) {
    diff = static_cast<double>(*iter) - mean_value;
    result += diff * diff;
  }

  // No error checking done here as
  return sqrt(result / static_cast<double>(values.size()-1));
}

/// Find the median of an std container
template <typename T>
double median(T const& values) {

  // Copy the input list to a vector
  const size_t num_elements = values.size();
  std::vector<double> v(num_elements);
  typename T::const_iterator iter = values.begin();
  for (size_t i=0; i<num_elements; ++i) {
    v[i] = *iter;
    ++iter;
  }

  // A copy is needed so we can preserve the order of the input list.
  std::sort(v.begin(), v.end());
  size_t center = num_elements/2;
  if (num_elements % 2 == 1) // Odd length case
    return v[center];
  else  // Even length case
    return (v[center] + v[center-1]) / 2.0;
}

//---------------------------------------------------------------------------------

template <class ContainerT>
ContainerT MeanFunctor::operator() (std::vector<ContainerT> const& points) const {
  ContainerT result = points[0]; // to resize container if necessary
  size_t num_points = points.size();
  size_t dimensions = points[0].size();
  size_t last       = points[0].size() - 1;
  if (m_homogeneous)
    dimensions--;

  for (size_t i = 0; i < dimensions; ++i)
    result[i] = 0;
  if (m_homogeneous)
    result[last] = 1;

  if (m_homogeneous) {
    for (size_t i = 0; i < num_points; ++i)
      for (size_t j = 0; j < dimensions; ++j)
        result[j] += points[i][j] / points[i][last];
  }
  else {
    for (size_t i = 0; i < num_points; ++i)
      for (size_t j = 0; j < dimensions; ++j)
        result[j] += points[i][j];
  }

  for (size_t i = 0; i < dimensions; ++i)
    result[i] /= num_points;

  return result;
}

template <class ContainerT>
ContainerT StandardDeviationFunctor::operator() (std::vector<ContainerT> const& points) const {
  ContainerT result = points[0]; // to resize container if necessary
  ContainerT temp = points[0]; // to resize container if necessary
  MeanFunctor mean_func(m_homogeneous);
  ContainerT mean = mean_func(points);
  unsigned num_points = points.size();
  unsigned dimensions = points[0].size();
  unsigned last = points[0].size() - 1;
  if (m_homogeneous)
    dimensions--;

  for (unsigned int i = 0; i < dimensions; ++i)
    result[i] = 0;
  if (m_homogeneous)
    result[last] = 1;

  if (m_homogeneous) {
    for (unsigned i = 0; i < num_points; ++i)
      for (unsigned int j = 0; j < dimensions; ++j) {
        temp[j] = points[i][j] / points[i][last] - mean[j];
        result[j] += temp[j] * temp[j];
      }
  }
  else {
    for (unsigned i = 0; i < num_points; ++i)
      for (unsigned int j = 0; j < dimensions; ++j) {
        temp[j] = points[i][j] - mean[j];
        result[j] += temp[j] * temp[j];
      }
  }

  for (unsigned int i = 0; i < dimensions; ++i) {
    result[i] /= num_points;
    result[i] = sqrt(result[i]);
  }

  return result;
}

// See the header file for the doc
template<class T>
bool find_outlier_brackets(std::vector<T> const& vec,
                           double pct_factor, double outlier_factor,
                           double & b, double & e) {
  if (pct_factor > 0.5)
    pct_factor = 1.0 - pct_factor;
  if (pct_factor < 0)
    pct_factor = 0;

  b = 0.0; e = 0.0; // initialize
  std::vector<T> q = vec;
  std::sort(q.begin(), q.end());

  int len = q.size();
  if (len <= 0) return false;

  b = q[0]; e = q[len-1];
  if (len <= 3) return false; // too few points for analysis

  int bn = int(round(pct_factor*len));
  int en = int(round((1.0-pct_factor)*len))-1;

  // Adjust the bounds just in case
  bn = std::min(len-1, std::max(bn, 0));
  en = std::min(len-1, std::max(en, 0));

  double Q1 = q[bn];
  double Q3 = q[en];
  b = Q1 - outlier_factor*(Q3 - Q1);
  e = Q3 + outlier_factor*(Q3 - Q1);

  return true;
}


// Histogram methods are non-template and live in Statistics.cc.

}} // namespace vw::math

#endif // __MATH_STATISTICS_H__
