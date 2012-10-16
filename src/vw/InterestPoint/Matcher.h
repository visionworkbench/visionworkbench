// __BEGIN_LICENSE__
//  Copyright (c) 2006-2012, United States Government as represented by the
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

#include <algorithm>

#include <vw/Core/Log.h>
#include <vw/InterestPoint/Descriptor.h>
#include <vector>
#include <boost/foreach.hpp>

#if VW_HAVE_PKG_FLANN
#include <vw/Math/FLANNTree.h>
#else
#include <vw/Math/KDTree.h>
#endif

namespace vw {
namespace ip {

  /// Interest point metrics must be of form:
  ///
  /// float operator() (const InterestPoint& ip1, const InterestPoint &ip2, float maxdist = DBL_MAX)

  /// L2 Norm: returns Euclidian distance squared between pair of
  /// interest point descriptors.  Optional argument "maxdist" to
  /// provide early termination of computation if result will exceed
  /// maxdist.
  struct L2NormMetric {
    float operator() (InterestPoint const& ip1, InterestPoint const& ip2,
                      float maxdist = std::numeric_limits<float>::max()) const;
  };

  /// KL distance (relative entropy) between interest point
  /// descriptors. Optional argument "maxdist" to provide early
  /// termination of computation if result will exceed maxdist.
  struct RelativeEntropyMetric {
    float operator() (InterestPoint const& ip1, InterestPoint const& ip2,
                      float maxdist = std::numeric_limits<float>::max()) const;
  };

  template <class ListT>
  inline void sort_interest_points(ListT const& ip1, ListT const& ip2,
                                   std::vector<ip::InterestPoint> & ip1_sorted,
                                   std::vector<ip::InterestPoint> & ip2_sorted){
    // The interest points may be in random order due to the fact that
    // they are obtained using multiple threads. Must sort them to
    // restore a unique order for the matcher to give a unique result.
    ip1_sorted.reserve( ip1.size() ); ip1_sorted.clear();
    ip2_sorted.reserve( ip2.size() ); ip2_sorted.clear();
    BOOST_FOREACH( ip::InterestPoint const& ip, ip1 )
      ip1_sorted.push_back(ip);
    BOOST_FOREACH( ip::InterestPoint const& ip, ip2 )
      ip2_sorted.push_back(ip);
    std::sort(ip1_sorted.begin(), ip1_sorted.end(), InterestPointLessThan);
    std::sort(ip2_sorted.begin(), ip2_sorted.end(), InterestPointLessThan);
  }
  
  /// Interest Point Match contraints functors to return a list of
  /// allowed match candidates to an interest point.
  ///
  /// Must have be of form: bool operator()( const InterestPoint&
  ///     baseline_ip, const InterestPoint& test_ip); returning true
  ///     or false depending on whether the test_ip is satisfies the
  ///     constraints relative to basline_ip..
  ///
  /// To use bimatch, must also include the inverse function: bool
  ///     inverse ()( const InterestPoint& baseline_ip, const
  ///     InterestPoint& test_ip) returning true or false under the
  ///     inverse constraint mapping.
  ///
  /// This default matching constraint functor providing NO
  /// constraints on putative matches.
  class NullConstraint {
  public:
    bool operator()(InterestPoint const& /*baseline_ip*/, InterestPoint const& /*test_ip*/) const {
      return true;
    }

    bool inverse(InterestPoint const& /*baseline_ip*/, InterestPoint const& /*test_ip*/) const {
      return true;
    }
  };

  /// Constraint on interest point orientation and scale differences.
  class ScaleOrientationConstraint {
  public:
    double scale_ratio_min;
    double scale_ratio_max;
    double ori_diff_min;
    double ori_diff_max;

    ScaleOrientationConstraint(double srmin = 0.9, double srmax = 1.1,
                               double odmin = -0.1, double odmax = 0.1) :
      scale_ratio_min(srmin), scale_ratio_max(srmax),
      ori_diff_min(odmin), ori_diff_max(odmax) {}

    bool operator()(InterestPoint const& baseline_ip, InterestPoint const& test_ip) const;

    bool inverse(InterestPoint const& baseline_ip, InterestPoint const& test_ip) const {
      return this->operator()(test_ip, baseline_ip);
    }
  };

  /// Constraint on interest point location differences.
  class PositionConstraint {
  public:
    double min_x;
    double max_x;
    double min_y;
    double max_y;

    PositionConstraint(double _min_x = -10.0, double _max_x = 10.0,
                       double _min_y = -10.0, double _max_y = 10.0) :
      min_x(_min_x), max_x(_max_x), min_y(_min_y), max_y(_max_y) {}

    bool operator()(InterestPoint const& baseline_ip, InterestPoint const& test_ip) const;

    bool inverse(InterestPoint const& baseline_ip, InterestPoint const& test_ip) const {
      return this->operator()(test_ip, baseline_ip);
    }
  };

  // ---------------------------------------------------------------------------
  //                         Interest Point Matcher
  // ---------------------------------------------------------------------------

  /// Interest point matcher class
  template < class MetricT, class ConstraintT >
  class InterestPointMatcher {
    ConstraintT m_constraint;
    MetricT m_distance_metric;
    double m_threshold;

  public:

    InterestPointMatcher(double threshold = 0.5, MetricT metric = MetricT(), ConstraintT constraint = ConstraintT())
      : m_constraint(constraint), m_distance_metric(metric), m_threshold(threshold) { }

    /// Given two lists of interest points, this routine returns the two lists
    /// of matching interest points based on the Metric and Constraints
    /// provided by the user.
    template <class ListT, class MatchListT>
    void operator()( ListT const& ip1, ListT const& ip2,
                     MatchListT& matched_ip1, MatchListT& matched_ip2,
                     bool bidirectional = false,
                     const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) const {
      typedef typename ListT::const_iterator IterT;

      Timer *total = new Timer("Total elapsed time", DebugMessage, "interest_point");

      std::vector<ip::InterestPoint> ip1_sorted, ip2_sorted;
      sort_interest_points(ip1, ip2, ip1_sorted, ip2_sorted);

      matched_ip1.clear(); matched_ip2.clear();
      if (!ip1_sorted.size() || !ip2_sorted.size()) {
        vw_out(InfoMessage,"interest_point") << "KD-Tree: no points to match, exiting\n";
        progress_callback.report_finished();
        return;
      }

      size_t size = ip1_sorted.size();
      float inc_amt = 1.0f/float(size);

#if VW_HAVE_PKG_FLANN
      Matrix<float> ip2_matrix( ip2_sorted.size(), ip2_sorted.begin()->size() );
      Matrix<float>::iterator ip2_matrix_it = ip2_matrix.begin();
      BOOST_FOREACH( InterestPoint const& ip, ip2_sorted )
        ip2_matrix_it = std::copy( ip.begin(), ip.end(), ip2_matrix_it );

      math::FLANNTree<flann::L2<float> > kd( ip2_matrix );
      vw_out(InfoMessage,"interest_point") << "FLANN-Tree created. Searching...\n";

      Vector<int> indices(2);
      Vector<float> distances(2);
#else
      math::KDTree<ListT> kd(ip2_sorted.begin()->size(), ip2_sorted);
      vw_out(InfoMessage,"interest_point") << "KD-Tree created with " << kd.size() << " nodes and depth ranging from " << kd.min_depth() << " to " << kd.max_depth() << ".  Searching...\n";
#endif
      progress_callback.report_progress(0);

      BOOST_FOREACH( InterestPoint ip, ip1_sorted ) {
        if (progress_callback.abort_requested())
          vw_throw( Aborted() << "Aborted by ProgressCallback" );
        progress_callback.report_incremental_progress(inc_amt);

        std::vector<InterestPoint> nearest_records(2);
#if VW_HAVE_PKG_FLANN
        kd.knn_search( ip.descriptor, indices, distances, 2 );
        nearest_records[0] = ip2_sorted[indices[0]];
        nearest_records[1] = ip2_sorted[indices[1]];
#else
        int num_records = kd.m_nearest_neighbors(ip, nearest_records, 2);
        if (num_records != 2)
          continue; // Ignore if there are no matches
#endif

        bool constraint_satisfied = false;
        if (bidirectional) {
          if (m_constraint(nearest_records[0], ip) &&
              m_constraint(ip, nearest_records[0]))
            constraint_satisfied = true;
        } else {
          if (m_constraint(nearest_records[0], ip))
            constraint_satisfied = true;
        }
        if (constraint_satisfied) {
          double dist0 = m_distance_metric(nearest_records[0], ip);
          double dist1 = m_distance_metric(nearest_records[1], ip);

          if (dist0 < m_threshold * dist1) {
            matched_ip1.push_back(ip);
            matched_ip2.push_back(nearest_records[0]);
          }
        }
      }

      progress_callback.report_finished();
      delete total;
    }
  };

  // Specialization to capture more speed.
  template <>
  class InterestPointMatcher<L2NormMetric, NullConstraint> {
    double m_threshold;
    L2NormMetric m_metric;

  public:

    InterestPointMatcher(double threshold = 0.5, L2NormMetric metric = L2NormMetric(), NullConstraint /*constraint*/ = NullConstraint())
      : m_threshold(threshold), m_metric(metric) { }

    /// Given two lists of interest points, this routine returns the two lists
    /// of matching interest points based on the Metric and Constraints
    /// provided by the user.
    template <class ListT, class MatchListT>
    void operator()( ListT const& ip1, ListT const& ip2,
                     MatchListT& matched_ip1, MatchListT& matched_ip2,
                     bool /*bidirectional*/ = false,
                     const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) const {
      typedef typename ListT::const_iterator IterT;

      Timer *total = new Timer("Total elapsed time", DebugMessage, "interest_point");

      std::vector<ip::InterestPoint> ip1_sorted, ip2_sorted;
      sort_interest_points(ip1, ip2, ip1_sorted, ip2_sorted);
      
      matched_ip1.clear(); matched_ip2.clear();
      if (!ip1_sorted.size() || !ip2_sorted.size()) {
        vw_out(InfoMessage,"interest_point") << "KD-Tree: no points to match, exiting\n";
        progress_callback.report_finished();
        return;
      }

      size_t size = ip1_sorted.size();
      float inc_amt = 1.0f/float(size);

#if VW_HAVE_PKG_FLANN
      Matrix<float> ip2_matrix( ip2_sorted.size(), ip2_sorted.begin()->size() );
      Matrix<float>::iterator ip2_matrix_it = ip2_matrix.begin();
      BOOST_FOREACH( InterestPoint const& ip, ip2_sorted )
        ip2_matrix_it = std::copy( ip.begin(), ip.end(), ip2_matrix_it );

      math::FLANNTree<flann::L2<float> > kd( ip2_matrix );
      vw_out(InfoMessage,"interest_point") << "FLANN-Tree created. Searching...\n";

      Vector<int> indices(2);
      Vector<float> distances(2);
#else
      math::KDTree<ListT> kd(ip2_sorted.begin()->size(), ip2_sorted);
      vw_out(InfoMessage,"interest_point") << "KD-Tree created with " << kd.size() << " nodes and depth ranging from " << kd.min_depth() << " to " << kd.max_depth() << ".  Searching...\n";
#endif
      progress_callback.report_progress(0);

      BOOST_FOREACH( InterestPoint ip, ip1_sorted ) {
        if (progress_callback.abort_requested())
          vw_throw( Aborted() << "Aborted by ProgressCallback" );
        progress_callback.report_incremental_progress(inc_amt);


#if VW_HAVE_PKG_FLANN
        kd.knn_search( ip.descriptor, indices, distances, 2 );
        if ( distances[0] < m_threshold * distances[1] ) {
          matched_ip1.push_back(ip);
          matched_ip2.push_back(ip2_sorted[indices[0]]);
        }
#else
        std::vector<InterestPoint> nearest_records(2);
        int num_records = kd.m_nearest_neighbors(ip, nearest_records, 2);
        if (num_records != 2)
          continue; // Ignore if there are no matches
        double dist0 = m_metric(nearest_records[0], ip);
        double dist1 = m_metric(nearest_records[1], ip);

        if (dist0 < m_threshold * dist1) {
          matched_ip1.push_back(ip);
          matched_ip2.push_back(nearest_records[0]);
        }
#endif
      }

      progress_callback.report_finished();
      delete total;
    }
  };

  // A even more basic interest point matcher that doesn't rely on
  // KDTree as it sometimes produces incorrect results
  template < class MetricT, class ConstraintT >
  class InterestPointMatcherSimple {
    ConstraintT m_constraint;
    MetricT m_distance_metric;
    double m_threshold;

  public:

  InterestPointMatcherSimple(double threshold = 0.5, MetricT metric = MetricT(), ConstraintT constraint = ConstraintT())
    : m_constraint(constraint), m_distance_metric(metric), m_threshold(threshold) { }

    /// Given two lists of interest points, this routine returns the two lists
    /// of matching interest points based on the Metric and Constraints
    /// provided by the user.
    template <class ListT, class MatchListT>
    void operator()( ListT const& ip1, ListT const& ip2,
                     MatchListT& matched_ip1, MatchListT& matched_ip2,
                     bool /*bidirectional*/ = false,
                     const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) const {
      typedef typename ListT::const_iterator IterT;

      Timer *total = new Timer("Total elapsed time", DebugMessage, "interest_point");

      std::vector<ip::InterestPoint> ip1_sorted, ip2_sorted;
      sort_interest_points(ip1, ip2, ip1_sorted, ip2_sorted);

      matched_ip1.clear(); matched_ip2.clear();
      if (!ip1_sorted.size() || !ip2_sorted.size()) {
        vw_out(InfoMessage,"interest_point") << "No points to match, exiting\n";
        progress_callback.report_finished();
        return;
      }

      float inc_amt = 1.0f / float(ip1_sorted.size());
      progress_callback.report_progress(0);
      std::vector<int> match_index( ip1_sorted.size() );
      for (unsigned i = 0; i < ip1_sorted.size(); i++ ) {
        if (progress_callback.abort_requested())
          vw_throw( Aborted() << "Aborted by ProgressCallback" );
        progress_callback.report_incremental_progress(inc_amt);

        double first_pick = 1e100, second_pick = 1e100;
        match_index[i] = -1;
        // Comparing ip1_sorted's feature against all of ip2_sorted's
        for (unsigned j = 0; j < ip2_sorted.size(); j++ ) {

          if ( ip1_sorted[i].polarity != ip2_sorted[j].polarity )
            continue;

          double distance = m_distance_metric( ip1_sorted[i], ip2_sorted[j] );

          if ( distance < first_pick ) {
            match_index[i] = j;
            second_pick = first_pick;
            first_pick = distance;
          } else if ( distance < second_pick ) {
            second_pick = distance;
          }

        }

        // Checking to see if the match is strong enough
        if ( first_pick > m_threshold * second_pick )
          match_index[i] = -1;
      }

      progress_callback.report_finished();

      // Building matched_ip1 & matched_ip2
      for (unsigned i = 0; i < ip1_sorted.size(); i++ ) {
        if ( match_index[i] != -1 ) {
          matched_ip1.push_back( ip1_sorted[i] );
          matched_ip2.push_back( ip2_sorted[match_index[i]] );
        }
      }

      delete total;
    }
  };


  //////////////////////////////////////////////////////////////////////////////
  // Convenience Typedefs
  typedef InterestPointMatcher< L2NormMetric, NullConstraint > DefaultMatcher;
  typedef InterestPointMatcher< L2NormMetric, ScaleOrientationConstraint > ConstraintedMatcher;

  // Matching doesn't constraint a point to being matched to only one
  // other point. Here's a way to remove duplicates and have only
  // pairwise points.
  void remove_duplicates(std::vector<InterestPoint>& ip1,
                         std::vector<InterestPoint>& ip2);

}} // namespace vw::ip

#endif // _INTEREST_POINT_MATCHER_H_
