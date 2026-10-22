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


/// \file Matcher.h
///
/// Classes and functions for matching image interest points.
///
#ifndef _INTERESTPOINT_MATCHER_H_
#define _INTERESTPOINT_MATCHER_H_

#include <vw/Core/Log.h>
#include <vw/InterestPoint/Descriptor.h>
#include <vw/InterestPoint/InterestPoint.h>
#include <vector>
#include <boost/foreach.hpp>

#ifdef VW_HAVE_PKG_FLANN
#include <vw/Math/FLANNTree.h>
#else
#error the FLANN library is now required!
#endif

namespace vw {
namespace ip {

  //======================================================================
  // Interest Point distance metrics

  /// Interest point metrics must be of form:
  ///
  /// float operator() (const InterestPoint& ip1, const InterestPoint &ip2, float maxdist = DBL_MAX)
  ///
  /// --> This one is for interoperability with our FLANNTRee class which does all our heavy-duty matching.
  /// static const math::FLANN_DistType flann_type=FLANN_DistType;

  /// L2 Norm: returns Euclidean distance squared between pair of
  /// interest point descriptors.  Optional argument "maxdist" to
  /// provide early termination of computation if result will exceed maxdist.
  struct L2NormMetric {
    float operator() (InterestPoint const& ip1, InterestPoint const& ip2,
                      float maxdist = std::numeric_limits<float>::max()) const;
    static const math::FLANN_DistType flann_type = math::FLANN_DistType_L2;
  };

  /// Hamming distance for binary descriptors
  struct HammingMetric {
    float operator() (InterestPoint const& ip1, InterestPoint const& ip2,
                      float maxdist = std::numeric_limits<float>::max()) const;
    static const math::FLANN_DistType flann_type = math::FLANN_DistType_Hamming;
  };

  // ---------------------------------------------------------------------------
  //                         Interest Point Matcher
  // ---------------------------------------------------------------------------

  /// Interest point matcher class
  template < class MetricT >
  class InterestPointMatcher {
    std::string m_flann_method;
    MetricT     m_distance_metric;
    double      m_threshold;

  public:

    InterestPointMatcher(std::string const& flann_method,
                         double threshold = 0.5,
                         MetricT metric = MetricT()):
    m_flann_method(flann_method), m_distance_metric(metric),
    m_threshold(threshold) {}

    /// Write to index_list the matching index in ip2 for each point in ip1.
    /// index_list is the same length as ip1, with the max value of size_t
    /// where no match was found.
    void operator()(std::vector<InterestPoint> const& ip1, std::vector<InterestPoint> const& ip2,
                    std::vector<size_t>& index_list,
                    const ProgressCallback &progress_callback
                      = ProgressCallback::dummy_instance(),
                    bool quiet = false) const;

    /// Return the two lists of matching interest points based on the metric.
    void operator()( std::vector<InterestPoint> const& ip1, std::vector<InterestPoint> const& ip2,
                     std::vector<InterestPoint>& matched_ip1, std::vector<InterestPoint>& matched_ip2,
                     const ProgressCallback &progress_callback = ProgressCallback::dummy_instance(),
                     bool quiet = false) const;
  };


  /// A even more basic interest point matcher that doesn't rely on
  /// KDTree as it sometimes produces incorrect results
  template < class MetricT >
  class InterestPointMatcherSimple {
    MetricT     m_distance_metric;
    double      m_threshold;

  public:

  InterestPointMatcherSimple(double threshold = 0.5, MetricT metric = MetricT())
    : m_distance_metric(metric), m_threshold(threshold) { }

    /// Return the two lists of matching interest points based on the metric.
    void operator()( std::vector<InterestPoint> const& ip1, std::vector<InterestPoint> const& ip2,
                     std::vector<InterestPoint>& matched_ip1, std::vector<InterestPoint>& matched_ip2,
                     const ProgressCallback &progress_callback = ProgressCallback::dummy_instance(),
                     bool quiet = false ) const;
  };


  //////////////////////////////////////////////////////////////////////////////
  // Convenience typedef
  typedef InterestPointMatcher<L2NormMetric> DefaultMatcher;

  /// Matching doesn't constraint a point to being matched to only one
  /// other point. Here's a way to remove duplicates and have only pairwise points.
  void remove_duplicates(std::vector<InterestPoint>& ip1,
                         std::vector<InterestPoint>& ip2);


}} // namespace vw::ip

#endif // _INTEREST_POINT_MATCHER_H_
