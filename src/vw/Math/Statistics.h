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

namespace vw {
namespace math {

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
    ContainerT operator() (std::vector<ContainerT> const& points) const {
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
    ContainerT operator() (std::vector<ContainerT> const& points) const {
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
  };


  template<class T>
  void find_outlier_brackets(std::vector<T> const& p,
                             double pct, double outlier_factor,
                             double & b, double & e){
  
    // Using the quartile range, determine the values b and e so that
    // all elements in p outsize of [b, e] are outliers.

    // If Q1 and Q3 are the percentiles for pct and 1-pct,
    // the outlier brackets are
    // b = Q1 - outlier_factor*(Q3-Q1)
    // e = Q3 + outlier_factor*(Q3-Q1)

    // This algorithm works best if the data is distributed rather
    // uniformly.
  
    // Suggested values for the inputs:
    // pct = 0.25, outlier_factor: 1.5.
    // We expect 0 < pct < 0.5.
  
    b = 0.0; e = 0.0; // initialize
    std::vector<T> q = p;
    std::sort(q.begin(), q.end());
    int len = q.size();
    if (len <= 0) return;
    b = q[0]; e = q[len-1];
    if (len <= 3) return; // too few points for analysis
      
    int bn = int(round(pct*len));
    int en = int(round((1.0-pct)*len))-1;
      
    double Q1 = q[bn];
    double Q3 = q[en];
    b = Q1 - outlier_factor*(Q3-Q1);
    e = Q3 + outlier_factor*(Q3-Q1);

    return;
  }
  
}} // namespace vw::math

#endif // __MATH_STATISTICS_H__
