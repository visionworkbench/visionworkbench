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
#include <vw/InterestPoint/Matcher.h>
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;

namespace vw {
namespace ip {

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


  /// Simple, unoptimized code for computing the hamming distance of two bytes.
  size_t hamming_helper(unsigned char a, unsigned char b) {
      int dist = 0;
      unsigned char val = a ^ b; // XOR

      // Count the number of bits set
      while (val != 0) {
          // A bit is set, so increment the count and clear the bit
          dist++;
          val &= val - 1;
      }
      return dist; // Return the number of differing bits
  }

  // This is a very slow Hamming distance implementation which can be replaced if needed.
  float HammingMetric::operator()( InterestPoint const& ip1, 
                                   InterestPoint const& ip2,
                                   float maxdist ) const {
    float dist = 0.0;
    for (size_t i = 0; i < ip1.descriptor.size(); i++) {
      // Cast the two elements to bytes which they should have originally been
      unsigned char desc1 = static_cast<unsigned char>(ip1.descriptor[i]);
      unsigned char desc2 = static_cast<unsigned char>(ip2.descriptor[i]);

      // Compute the hamming distance between just these two bytes
      size_t dist_int = hamming_helper(desc1, desc2);

      // Accumulate the floating point distance
      dist += static_cast<float>(dist_int);

      if (dist > maxdist) break;  // abort calculation if distance exceeds upper bound
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


//==================================================================================

  void remove_duplicates(std::vector<InterestPoint>& ip1,
                         std::vector<InterestPoint>& ip2) {
    VW_ASSERT( ip1.size() == ip2.size(),
               ArgumentErr() << "Input vectors are not the same size.");
    std::vector<InterestPoint> ip1_fltr, ip2_fltr;
    ip1_fltr.reserve( ip1.size() );
    ip2_fltr.reserve( ip2.size() );

    for ( size_t i = 0; i < ip1.size(); ++i ) {
      bool bad_entry = false;
      for ( size_t j = i + 1; j < ip1.size(); ++j ) {
        if ( (ip1[i].x == ip1[j].x && ip1[i].y == ip1[j].y) ||
             (ip2[i].x == ip2[j].x && ip2[i].y == ip2[j].y) ) {
          bad_entry = true;
          break;
        }
      }
      if (!bad_entry) {
        ip1_fltr.push_back( ip1[i] );
        ip2_fltr.push_back( ip2[i] );
      }
    }
    ip1 = ip1_fltr;
    ip2 = ip2_fltr;
  }

  std::string strip_path(std::string out_prefix, std::string filename){

    // If filename starts with out_prefix followed by dash, strip both.
    // Also strip filename extension.

    std::string ss = out_prefix + "-";
    size_t found = filename.find(ss);

    if (found != std::string::npos)
      filename.erase(found, ss.length());

    filename = fs::path(filename).stem().string();

    return filename;
  }


  std::string match_filename(std::string const& out_prefix,
                             std::string const& input_file1,
                             std::string const& input_file2){

    // filenames longer than this must be chopped, as too long names
    // cause problems later with boost.
    int max_len = 40;
    std::string name1 = strip_path(out_prefix, input_file1).substr(0, max_len);
    std::string name2 = strip_path(out_prefix, input_file2).substr(0, max_len);

    return out_prefix + "-" + name1 + "__" + name2 + ".match";
  }

  void ip_filenames(std::string const& out_prefix,
                    std::string const& input_file1,
                    std::string const& input_file2,
                    std::string & output_ip1,
                    std::string & output_ip2
                    ){
    output_ip1 = out_prefix + "-" + strip_path(out_prefix, input_file1) + ".vwip";
    output_ip2 = out_prefix + "-" + strip_path(out_prefix, input_file2) + ".vwip";
  }

}} // namespace vw::ip
