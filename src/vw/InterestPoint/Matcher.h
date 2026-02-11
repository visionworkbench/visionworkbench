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

#include <algorithm>

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

  /// KL distance (relative entropy) between interest point
  /// descriptors. Optional argument "maxdist" to provide early
  /// termination of computation if result will exceed maxdist.
  struct RelativeEntropyMetric {
    float operator() (InterestPoint const& ip1, InterestPoint const& ip2,
                      float maxdist = std::numeric_limits<float>::max()) const;
    static const math::FLANN_DistType flann_type = math::FLANN_DistType_Unsupported;
  };


  //======================================================================

  template <class ListT>
  inline void sort_interest_points(ListT const& ip1, ListT const& ip2,
                                   std::vector<ip::InterestPoint> & ip1_sorted,
                                   std::vector<ip::InterestPoint> & ip2_sorted);

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
    std::string m_flann_method; 
    ConstraintT m_constraint;
    MetricT     m_distance_metric;
    double      m_threshold;
    bool        m_bidirectional;

    // Helper function to help reduce conditionals in the event of
    // NullConstraint. (Which is common).
    template <class InConstraintT>
    typename boost::disable_if<boost::is_same<InConstraintT,NullConstraint>, bool>::type
    check_constraint( InterestPoint const& ip1, InterestPoint const& ip2 ) const;

    template <class InConstraintT>
    typename boost::enable_if<boost::is_same<InConstraintT,NullConstraint>, bool>::type
    check_constraint( InterestPoint const& ip1, InterestPoint const& ip2 ) const {
      return true;
    }

  public:

    InterestPointMatcher(std::string const& flann_method, 
                         double threshold = 0.5, 
                         MetricT metric = MetricT(), ConstraintT constraint = ConstraintT(), bool bidirectional = false): 
    m_flann_method(flann_method), m_constraint(constraint), m_distance_metric(metric),
    m_threshold(threshold), m_bidirectional(bidirectional) {}

    /// Given two lists of interest points, this write to index_list
    /// the corresponding matching index in ip2. index_list is the
    /// same length as ip1. index_list will be filled with max value
    /// of size_t in the event that a match was not found.
    template <class ListT, class IndexListT >
    void operator()(ListT const& ip1, ListT const& ip2,
                    IndexListT& index_list,
                    const ProgressCallback &progress_callback 
                      = ProgressCallback::dummy_instance(), 
                    bool quiet = false) const;

    /// Given two lists of interest points, this routine returns the two lists
    /// of matching interest points based on the Metric and Constraints provided by the user.
    template <class ListT, class MatchListT>
    void operator()( ListT const& ip1, ListT const& ip2,
                     MatchListT& matched_ip1, MatchListT& matched_ip2,
                     const ProgressCallback &progress_callback = ProgressCallback::dummy_instance(), 
                     bool quiet = false) const;
  };


  /// A even more basic interest point matcher that doesn't rely on
  /// KDTree as it sometimes produces incorrect results
  template < class MetricT, class ConstraintT >
  class InterestPointMatcherSimple {
    ConstraintT m_constraint;
    MetricT     m_distance_metric;
    double      m_threshold;

  public:

  InterestPointMatcherSimple(double threshold = 0.5, MetricT metric = MetricT(), ConstraintT constraint = ConstraintT())
    : m_constraint(constraint), m_distance_metric(metric), m_threshold(threshold) { }

    /// Given two lists of interest points, this routine returns the two lists
    /// of matching interest points based on the Metric and Constraints provided by the user.
    template <class ListT, class MatchListT>
    void operator()( ListT const& ip1, ListT const& ip2,
                     MatchListT& matched_ip1, MatchListT& matched_ip2,
                     const ProgressCallback &progress_callback = ProgressCallback::dummy_instance(), 
                     bool quiet = false ) const;
  };


  //////////////////////////////////////////////////////////////////////////////
  // Convenience typedef
  typedef InterestPointMatcher<L2NormMetric, NullConstraint> DefaultMatcher;

  /// Matching doesn't constraint a point to being matched to only one
  /// other point. Here's a way to remove duplicates and have only pairwise points.
  void remove_duplicates(std::vector<InterestPoint>& ip1,
                         std::vector<InterestPoint>& ip2);

//==========================================================================
// Fuction definitions

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

//---------------------------------------------------------------------------
// InterestPointMatcher

// Helper function to help reduce conditionals in the event of
// NullConstraint. (Which is common).
template <class MetricT, class ConstraintT>
template <class InConstraintT>
typename boost::disable_if<boost::is_same<InConstraintT,NullConstraint>, bool>::type
InterestPointMatcher<MetricT, ConstraintT>::
check_constraint( InterestPoint const& ip1, InterestPoint const& ip2 ) const {
  if (m_bidirectional) {
    if (m_constraint(ip2, ip1) &&
        m_constraint(ip1, ip2))
      return true;
  } else {
    if (m_constraint(ip1, ip2))
      return true;
  }
  return false;
}

// Given two lists of interest points, this write to index_list
// the corresponding matching index in ip2. index_list is the
// same length as ip1. index_list will be filled with max value
// of size_t in the event that a match was not found.
template <class MetricT, class ConstraintT>
template <class ListT, class IndexListT >
void InterestPointMatcher<MetricT, ConstraintT>::operator()
    (ListT const& ip1, ListT const& ip2,
    IndexListT& index_list,
    const ProgressCallback &progress_callback,
    bool quiet) const {

  Timer total_time("Total elapsed time", DebugMessage, "interest_point");
  size_t ip1_size = ip1.size(), ip2_size = ip2.size();

  index_list.clear();
  if (!ip1_size || !ip2_size) {
    if (!quiet) {
      vw_out(InfoMessage,"interest_point") << "KD-Tree: no points to match, exiting\n";
      progress_callback.report_finished();
    }
    return;
  }

  float inc_amt = 1.0f/float(ip1_size);

  // Set up FLANNTree objects of all the different types we may need.
  math::FLANNTree<float>         kd_float(m_flann_method);
  math::FLANNTree<unsigned char> kd_uchar(m_flann_method);

  Matrix<float>         ip2_matrix_float;
  Matrix<unsigned char> ip2_matrix_uchar;

  // Pack the IP descriptors into a matrix and feed it to the chosen FLANNTree object
  const bool use_uchar_FLANN = (MetricT::flann_type == math::FLANN_DistType_Hamming);
  if (use_uchar_FLANN) {
    ip_list_to_matrix(ip2, ip2_matrix_uchar);
    kd_uchar.load_match_data(ip2_matrix_uchar, MetricT::flann_type);
  }else {
    ip_list_to_matrix(ip2, ip2_matrix_float);
    kd_float.load_match_data(ip2_matrix_float,  MetricT::flann_type);
  }

  if (!quiet)
    vw_out(InfoMessage,"interest_point") << "FLANN-Tree created. Searching...\n";

  const size_t KNN = 2; // Find this many matches
  Vector<int>    indices(KNN);
  Vector<double> distances(KNN);

  if (!quiet)
    progress_callback.report_progress(0);

  BOOST_FOREACH( InterestPoint ip, ip1 ) {
    if (progress_callback.abort_requested())
      vw_throw( Aborted() << "Aborted by ProgressCallback");

    if (!quiet)
      progress_callback.report_incremental_progress(inc_amt);

    size_t num_matches_found = 0;
    if (use_uchar_FLANN) {
      // Convert the descriptor to unsigned chars, then call FLANN
      vw::Vector<unsigned char> uchar_descriptor(ip.descriptor.size());
      for (size_t i=0; i<ip.descriptor.size(); i++)
        uchar_descriptor[i] = static_cast<unsigned char>(ip.descriptor[i]);
      num_matches_found = kd_uchar.knn_search( uchar_descriptor, indices, distances, KNN );
    }
    else {
      // Use float
      num_matches_found = kd_float.knn_search( ip.descriptor, indices, distances, KNN );
    }

    //vw_out() << "KNN matches "<< num_matches_found <<": indices = " << indices << ",    distances = " << distances << std::endl;

    if (num_matches_found < KNN) {
      // If we did not get two nearest neighbors, return no match for this point.
      //vw_out() << "Bad descriptor = " << ip.descriptor << std::endl;
      index_list.push_back( (size_t)(-1) ); // Last value of size_t
    } else {

      // Copy the two nearest matches
      std::vector<InterestPoint> nearest_records(KNN);
      typename ListT::const_iterator iterator = ip2.begin();
      std::advance( iterator, indices[0] );
      nearest_records[0] = *iterator;
      iterator = ip2.begin();
      std::advance( iterator, indices[1] );
      nearest_records[1] = *iterator;

      // Check the user constraint on the record
      if ( check_constraint<ConstraintT>( nearest_records[0], ip ) ) {
        double dist0 = m_distance_metric(nearest_records[0], ip);
        double dist1 = m_distance_metric(nearest_records[1], ip);

        // As a final check, make sure the nearest record is significantly closer than the next one.
        if (dist0 < m_threshold * dist1) {
          index_list.push_back( indices[0] );
        } else {
          index_list.push_back( (size_t)(-1) ); // Last value of size_t
        }
      } // End check constraint
    } // End both valid case
  }

  if (!quiet)
    progress_callback.report_finished();

} // End InterestPointMatcher::operator()

// Given two lists of interest points, this routine returns the two lists
// of matching interest points based on the Metric and Constraints
// provided by the user.
template <class MetricT, class ConstraintT>
template <class ListT, class MatchListT>
void InterestPointMatcher<MetricT, ConstraintT>::operator()
    (ListT const& ip1, ListT const& ip2,
     MatchListT& matched_ip1, MatchListT& matched_ip2,
     const ProgressCallback &progress_callback, 
     bool quiet) const {

  // Clear output lists
  matched_ip1.clear();
  matched_ip2.clear();

  // Redirect to the other version of this function, getting the results in an index list.
  std::list<size_t> index_list;
  this->operator()(ip1, ip2, index_list, progress_callback, quiet);

  // Now convert from the index output to the pairs output

  // Loop through ip1 and index_list
  std::list<size_t>::const_iterator index_list_iter = index_list.begin();
  BOOST_FOREACH( InterestPoint ip, ip1 ) {

    // Get and check the match index
    typename ListT::const_iterator ip2_iter = ip2.begin();
    size_t list_position = *index_list_iter;
    if (list_position < ip2.size()) {
      // Skip points without a match

      // Get the ip2 that corresponds to the current ip1 and store the point pair
      std::advance( ip2_iter, list_position);
      matched_ip1.push_back(ip);
      matched_ip2.push_back(*ip2_iter);
    }
    ++index_list_iter;
  } // End loop through ip1

}

//-----------------------------------------------------------
// InterestPointMatcherSimple

// Given two lists of interest points, this routine returns the two lists
// of matching interest points based on the Metric and Constraints
// provided by the user.
template <class MetricT, class ConstraintT>
template <class ListT, class MatchListT>
void InterestPointMatcherSimple<MetricT, ConstraintT>::operator()
  (ListT const& ip1, ListT const& ip2,
   MatchListT& matched_ip1, MatchListT& matched_ip2,
   const ProgressCallback &progress_callback, 
   bool quiet) const {

  Timer total_time("Total elapsed time", DebugMessage, "interest_point");

  matched_ip1.clear(); matched_ip2.clear();
  if (!ip1.size() || !ip2.size()) {
    if (!quiet) {
      vw_out(InfoMessage,"interest_point") << "No points to match, exiting\n";
      progress_callback.report_finished();
    }
    return;
  }

  std::vector<int> match_index( ip1.size() );

  float inc_amt = 1.0f / float(ip1.size());
  if (!quiet)
    progress_callback.report_progress(0);

  // Loop through one vector of IPs
  for (size_t i = 0; i < ip1.size(); i++ ) {
    if (progress_callback.abort_requested())
      vw_throw( Aborted() << "Aborted by ProgressCallback" );

    if (!quiet)
      progress_callback.report_incremental_progress(inc_amt);

    double first_pick = 1e100, second_pick = 1e100;
    match_index[i] = -1;
    // Comparing ip1_sorted's feature against all of ip2_sorted's
    for (size_t j = 0; j < ip2.size(); j++ ) {

      if ( ip1[i].polarity != ip2[j].polarity )
        continue;

      // Use the selected distance metric to measure the distance between the points
      double distance = m_distance_metric( ip1[i], ip2[j] );

      if ( distance < first_pick ) {
        match_index[i] = j;
        second_pick = first_pick;
        first_pick = distance;
      } else if ( distance < second_pick ) {
        second_pick = distance;
      }

    }

    //vw_out() << "Best 2 distances: " << first_pick <<", " << second_pick <<"\n";

    // Checking to see if the match is strong enough
    if ( first_pick > m_threshold * second_pick )
      match_index[i] = -1;
  } // End double loop through IPs

  if (!quiet)
    progress_callback.report_finished();

  // Building matched_ip1 & matched_ip2
  for (size_t i = 0; i < ip1.size(); i++ ) {
    if ( match_index[i] != -1 ) {
      matched_ip1.push_back( ip1[i] );
      matched_ip2.push_back( ip2[match_index[i]] );
    }
  }

  return;
} // end function operator()

}} // namespace vw::ip

#endif // _INTEREST_POINT_MATCHER_H_
