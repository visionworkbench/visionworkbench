// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    ContainerT operator() (std::vector<ContainerT> const& points) const {
      ContainerT result = points[0]; // to resize container if necessary
      size_t num_points = points.size();
      size_t dimensions = points[0].size();
      if (m_homogeneous)
        dimensions--;

      for (size_t i = 0; i < dimensions; ++i)
        result[i] = 0;
      if (m_homogeneous)
        result[dimensions] = 1;

      if (m_homogeneous) {
        for (size_t i = 0; i < num_points; ++i)
          for (size_t j = 0; j < dimensions; ++j)
            result[j] += points[i][j] / points[i][dimensions];
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
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    ContainerT operator() (std::vector<ContainerT> const& points) const {
      ContainerT result = points[0]; // to resize container if necessary
      ContainerT temp = points[0]; // to resize container if necessary
      MeanFunctor mean_func(m_homogeneous);
      ContainerT mean = mean_func(points);
      unsigned num_points = points.size();
      unsigned dimensions = points[0].size();
      if (m_homogeneous)
        dimensions--;

      for (unsigned int i = 0; i < dimensions; ++i)
        result[i] = 0;
      if (m_homogeneous)
        result[dimensions] = 1;

      if (m_homogeneous) {
        for (unsigned i = 0; i < num_points; ++i)
          for (unsigned int j = 0; j < dimensions; ++j) {
            temp[j] = points[i][j] / points[i][dimensions] - mean[j];
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

}} // namespace vw::math

#endif // __MATH_STATISTICS_H__
