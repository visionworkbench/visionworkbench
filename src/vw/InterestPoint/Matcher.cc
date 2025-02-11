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
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

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
  uint8 packed1[desc_len];
  uint8 packed2[desc_len];
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

float
RelativeEntropyMetric::operator()( InterestPoint const& ip1,
                                    InterestPoint const& ip2,
                                    float maxdist ) const {
  float dist = 0.0;
  for (size_t i = 0; i < ip1.descriptor.size(); i++) {
    dist += ip1.descriptor[i] * logf(ip1.descriptor[i]/(ip2.descriptor[i]+1e-16)+1e-16)/logf(2.) ;
    if (dist > maxdist) break;  // abort calculation if distance exceeds upper bound
  }
  return dist;
}

//==================================================================================
// Constraints

bool ScaleOrientationConstraint::operator()( InterestPoint const& baseline_ip,
                                              InterestPoint const& test_ip ) const {
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

bool PositionConstraint::operator()( InterestPoint const& baseline_ip,
                                      InterestPoint const& test_ip ) const {

  double dx = test_ip.x - baseline_ip.x;
  double dy = test_ip.y - baseline_ip.y;
  if (dx >= min_x && dx <= max_x && dy >= min_y && dy <= max_y) {
    return true;
  }

  // Otherwise...
  return false;
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

}} // namespace vw::ip
