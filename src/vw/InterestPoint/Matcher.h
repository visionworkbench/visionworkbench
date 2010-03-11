// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Matcher.h
///
/// Classes and functions for matching image interest points.
///
#ifndef _INTERESTPOINT_MATCHER_H_
#define _INTERESTPOINT_MATCHER_H_

#include <algorithm>

#include <vw/InterestPoint/Descriptor.h>
#include <vw/Math/KDTree.h>
#include <vw/Core/Log.h>
#include <vector>
#include <cmath>
#include <cfloat>

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
    double operator() (InterestPoint const& ip1, InterestPoint const& ip2, float maxdist = DBL_MAX) const {
      double dist = 0.0;
      for (unsigned int i = 0; i < ip1.descriptor.size(); i++) {
        dist += (ip1.descriptor[i] - ip2.descriptor[i])*(ip1.descriptor[i] - ip2.descriptor[i]);
        if (dist > maxdist) break;  // abort calculation if distance exceeds upper bound
      }
      return dist;
    }
  };

  /// KL distance (relative entropy) between interest point
  /// descriptors. Optional argument "maxdist" to provide early
  /// termination of computation if result will exceed maxdist.
  struct RelativeEntropyMetric {
    float operator() (InterestPoint const& ip1, InterestPoint const& ip2, float maxdist = DBL_MAX) const {
      float dist = 0.0;
      for (unsigned int i = 0; i < ip1.descriptor.size(); i++) {
        dist += ip1.descriptor[i] * logf(ip1.descriptor[i]/(ip2.descriptor[i]+1e-16)+1e-16)/logf(2.) ;
        if (dist > maxdist) break;  // abort calculation if distance exceeds upper bound
      }
      return dist;
    }
  };

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
                               double odmin = -0.1, double odmax = 0.1) {
      scale_ratio_min = srmin;
      scale_ratio_max = srmax;
      ori_diff_min = odmin;
      ori_diff_max = odmax;
    }

    bool operator()(InterestPoint const& baseline_ip, InterestPoint const& test_ip) const {

      double sr = test_ip.scale / baseline_ip.scale;
      double od = test_ip.orientation - baseline_ip.orientation;
      // Bring orientation delta (od) into range -M_PI to M_PI
      if (od < -M_PI) od += M_PI*2;
      else if (od > M_PI) od -= M_PI*2;

      if (sr >= scale_ratio_min && sr <= scale_ratio_max &&
          od >= ori_diff_min && od <= ori_diff_max) {
        return true;
      }

      // Otherwise...
      return false;
    }

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
                       double _min_y = -10.0, double _max_y = 10.0) {
      min_x = _min_x; max_x = _max_x; min_y = _min_y; max_y = _max_y;
    }

    bool operator()(InterestPoint const& baseline_ip, InterestPoint const& test_ip) const {

      double dx = test_ip.x - baseline_ip.x;
      double dy = test_ip.y - baseline_ip.y;
      if (dx >= min_x && dx <= max_x && dy >= min_y && dy <= max_y) {
        return true;
      }

      // Otherwise...
      return false;
    }

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

      matched_ip1.clear(); matched_ip2.clear();
      if (!ip1.size() || !ip2.size()) {
        vw_out(InfoMessage,"interest_point") << "KD-Tree: no points to match, exiting\n";
        progress_callback.report_finished();
        return;
      }


      math::KDTree<ListT> kd(ip2.begin()->descriptor.size(), ip2);
      vw_out(InfoMessage,"interest_point") << "KD-Tree created with " << kd.size() << " nodes and depth ranging from " << kd.min_depth() << " to " << kd.max_depth() << ".  Searching...\n";

      progress_callback.report_progress(0);

      int size = ip1.size();
      int n = 0;
      for (IterT iter = ip1.begin(); iter != ip1.end(); ++iter, ++n) {
        if (progress_callback.abort_requested())
          vw_throw( Aborted() << "Aborted by ProgressCallback" );
        progress_callback.report_progress(float(n)/float(size));

        std::vector<InterestPoint> nearest_records;
        int num_records = kd.m_nearest_neighbors(*iter, nearest_records, 2);
        if (num_records != 2)
          continue; // Ignore if there are no matches

        bool constraint_satisfied = false;
        if (bidirectional) {
          if (m_constraint(nearest_records[0], *iter) &&
              m_constraint(*iter, nearest_records[0]))
            constraint_satisfied = true;
        } else {
          if (m_constraint(nearest_records[0], *iter))
            constraint_satisfied = true;
        }

        if (constraint_satisfied) {
          double dist0 = m_distance_metric(nearest_records[0], *iter);
          double dist1 = m_distance_metric(nearest_records[1], *iter);

          if (dist0 < m_threshold * dist1) {
            matched_ip1.push_back(*iter);
            matched_ip2.push_back(nearest_records[0]);
          }
        }
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

      matched_ip1.clear(); matched_ip2.clear();
      if (!ip1.size() || !ip2.size()) {
        vw_out(InfoMessage,"interest_point") << "No points to match, exiting\n";
        progress_callback.report_finished();
        return;
      }

      progress_callback.report_progress(0);
      std::vector<int> match_index( ip1.size() );
      for (unsigned i = 0; i < ip1.size(); i++ ) {
        if (progress_callback.abort_requested())
          vw_throw( Aborted() << "Aborted by ProgressCallback" );
        progress_callback.report_progress(float(i)/float(ip1.size()));

        double first_pick = 1e100, second_pick = 1e100;
        match_index[i] = -1;
        // Comparing ip1's feature against all of ip2's
        for (unsigned j = 0; j < ip2.size(); j++ ) {

          if ( ip1[i].polarity != ip2[j].polarity )
            continue;

          double distance = m_distance_metric( ip1[i], ip2[j] );

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

      // Building matched_ip1 & matched ip 2
      for (unsigned i = 0; i < ip1.size(); i++ ) {
        if ( match_index[i] != -1 ) {
          matched_ip1.push_back( ip1[i] );
          matched_ip2.push_back( ip2[match_index[i]] );
        }
      }

      delete total;
    }
  };


//////////////////////////////////////////////////////////////////////////////
// Convenience Typedefs
typedef InterestPointMatcher< L2NormMetric, NullConstraint > DefaultMatcher;
typedef InterestPointMatcher< L2NormMetric, ScaleOrientationConstraint > ConstraintedMatcher;

}} // namespace vw::ip

#endif // _INTEREST_POINT_MATCHER_H_
