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

/// \file Matcher.cc
///
/// Classes and functions for matching image interest points.
///
#include <vw/InterestPoint/Matcher.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Core/Exception.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>

#include <set>

namespace vw { namespace ip {

//==================================================================================
// IP descriptor distance metrics

float
L2NormMetric::operator()( InterestPoint const& ip1, InterestPoint const& ip2,
                          float maxdist ) const {
                            
  float dist = 0.0;
  for (size_t i = 0; i < ip1.descriptor.size(); i++) {
    dist += (ip1.descriptor[i] - ip2.descriptor[i])*(ip1.descriptor[i] - ip2.descriptor[i]);
    if (dist > maxdist) break;  // abort calculation if distance exceeds upper bound
  }
  return dist;
}

// This could be simplified a lot if the description length was constant!
float HammingMetric::operator()( InterestPoint const& ip1, 
                                  InterestPoint const& ip2,
                                  float maxdist ) const {
  float dist = 0.0;
  
  // To use optimized code, first break the descriptor into units of 8, then 4, then 1.
  const int desc_len = static_cast<int>(ip1.descriptor.size());
  int d = desc_len;
  const int uint64_len = 8;
  const int uint32_len = 4;
  const int num_64 = d / uint64_len;
  d -= num_64*uint64_len;
  const int num_32 = d / uint32_len;
  
  // Pack the date into vectors as raw bytes (instead of integer values casted as floats)
  std::vector<uint8> packed1(desc_len), packed2(desc_len);
  
  for (int k=0; k<desc_len; ++k) {
    packed1[k] = static_cast<uint8>(ip1.descriptor[k]);
    packed2[k] = static_cast<uint8>(ip2.descriptor[k]);
  }
  
  int i=0;  // Never reset this.
  for (int i64=0; i64<num_64; ++i64) {
    // Compute the hamming distance between the next 8 bytes
    uint64 desc1    = *reinterpret_cast<uint64*>(&packed1[i]);
    uint64 desc2    = *reinterpret_cast<uint64*>(&packed2[i]);
    size_t dist_int = hamming_distance(desc1, desc2);
    i += uint64_len;

    // Accumulate the floating point distance
    dist += static_cast<float>(dist_int);
    if (dist > maxdist) return dist;  // abort calculation if distance exceeds upper bound
  }

  for (int i32=0; i32<num_32; ++i32) {
    // Compute the hamming distance between the next 4 bytes
    uint32 desc1    = *reinterpret_cast<uint32*>(&packed1[i]);
    uint32 desc2    = *reinterpret_cast<uint32*>(&packed2[i]);
    size_t dist_int = hamming_distance(desc1, desc2);
    i += uint32_len;

    // Accumulate the floating point distance
    dist += static_cast<float>(dist_int);
    if (dist > maxdist) return dist;  // abort calculation if distance exceeds upper bound
  }
  
  for (i=i; i<desc_len; ++i) {
    // Compute the hamming distance between any remaining bytes
    size_t dist_int = hamming_distance(packed1[i], packed2[i]);

    // Accumulate the floating point distance
    dist += static_cast<float>(dist_int);
  }
  return dist;
}

// Remove duplicates in interest point matches. Use a set for speed. When
// duplicates are encountered, keep the last repeated value. The code is this
// way to get the same results as a prior version of the code.
void remove_duplicates(std::vector<InterestPoint>& ip1,
                       std::vector<InterestPoint>& ip2) {
  
  // Sanity check
  VW_ASSERT(ip1.size() == ip2.size(),
            ArgumentErr() << "Input vectors are not the same size.");

  // Create a set of pairs of floats to be able to remove duplicates fast
  std::set<std::pair<float, float>> ip1_set, ip2_set;

  std::vector<InterestPoint> ip1_fltr, ip2_fltr;
  ip1_fltr.reserve(ip1.size());
  ip2_fltr.reserve(ip2.size());

  for (size_t i = 0; i < ip1.size(); i++) {
    
    // Make a pair for ip1 and ip2
    int j = ip1.size() - i - 1; // note the reverse order to keep last repeated one
    std::pair<float, float> ip1_pair(ip1[j].x, ip1[j].y);
    std::pair<float, float> ip2_pair(ip2[j].x, ip2[j].y);
    
    // If in the set, continue
    if (ip1_set.find(ip1_pair) != ip1_set.end() ||
        ip2_set.find(ip2_pair) != ip2_set.end())
      continue;
      
    // Add to the sets
    ip1_set.insert(ip1_pair);
    ip2_set.insert(ip2_pair);
      
    // New unique points, add to the filtered vectors  
    ip1_fltr.push_back(ip1[j]);
    ip2_fltr.push_back(ip2[j]);
  }
  
  // Undo the reverse order
  ip1.clear();
  ip2.clear();
  for (size_t i = 0; i < ip1_fltr.size(); i++)  {
    int j = ip1_fltr.size() - i - 1;
    ip1.push_back(ip1_fltr[j]);
    ip2.push_back(ip2_fltr[j]);
  }
}

//---------------------------------------------------------------------------
// InterestPointMatcher

// Write to index_list the matching index in ip2 for each point in ip1.
// index_list is the same length as ip1, with the max value of size_t where no
// match was found.
template <class MetricT>
void InterestPointMatcher<MetricT>::operator()
    (std::vector<InterestPoint> const& ip1, std::vector<InterestPoint> const& ip2,
     std::vector<size_t>& index_list,
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

    if (num_matches_found < KNN) {
      // If we did not get two nearest neighbors, return no match for this point.
      index_list.push_back( (size_t)(-1) ); // Last value of size_t
    } else {

      // Copy the two nearest matches
      std::vector<InterestPoint> nearest_records(KNN);
      auto iterator = ip2.begin();
      std::advance( iterator, indices[0] );
      nearest_records[0] = *iterator;
      iterator = ip2.begin();
      std::advance( iterator, indices[1] );
      nearest_records[1] = *iterator;

      double dist0 = m_distance_metric(nearest_records[0], ip);
      double dist1 = m_distance_metric(nearest_records[1], ip);

      // As a final check, make sure the nearest record is significantly closer than the next one.
      if (dist0 < m_threshold * dist1) {
        index_list.push_back( indices[0] );
      } else {
        index_list.push_back( (size_t)(-1) ); // Last value of size_t
      }
    } // End both valid case
  }

  if (!quiet)
    progress_callback.report_finished();

} // End InterestPointMatcher::operator()

// Return the two lists of matching interest points based on the metric.
template <class MetricT>
void InterestPointMatcher<MetricT>::operator()
    (std::vector<InterestPoint> const& ip1, std::vector<InterestPoint> const& ip2,
     std::vector<InterestPoint>& matched_ip1, std::vector<InterestPoint>& matched_ip2,
     const ProgressCallback &progress_callback,
     bool quiet) const {

  // Clear output lists
  matched_ip1.clear();
  matched_ip2.clear();

  // Redirect to the index-list version, then convert the indices to pairs.
  std::vector<size_t> index_list;
  this->operator()(ip1, ip2, index_list, progress_callback, quiet);

  // Loop through ip1 and index_list
  auto index_list_iter = index_list.begin();
  BOOST_FOREACH( InterestPoint ip, ip1 ) {

    // Get and check the match index
    auto ip2_iter = ip2.begin();
    size_t list_position = *index_list_iter;
    if (list_position < ip2.size()) {
      // Skip points without a match

      // Get the ip2 that corresponds to the current ip1 and store the point pair
      std::advance( ip2_iter, list_position);
      matched_ip1.push_back(ip);
      matched_ip2.push_back(*ip2_iter);
    }
    index_list_iter++;
  } // End loop through ip1

}

//-----------------------------------------------------------
// InterestPointMatcherSimple

// Return the two lists of matching interest points based on the metric.
template <class MetricT>
void InterestPointMatcherSimple<MetricT>::operator()
  (std::vector<InterestPoint> const& ip1, std::vector<InterestPoint> const& ip2,
   std::vector<InterestPoint>& matched_ip1, std::vector<InterestPoint>& matched_ip2,
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
    // Compare ip1's feature against all of ip2's
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

// Explicit instantiations for the two metrics actually used.
template class InterestPointMatcher<L2NormMetric>;
template class InterestPointMatcher<HammingMetric>;
template class InterestPointMatcherSimple<L2NormMetric>;
template class InterestPointMatcherSimple<HammingMetric>;

}} // namespace vw::ip
