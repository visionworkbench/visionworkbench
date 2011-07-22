// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

/// \file Matcher.h
///
/// Classes and functions for matching image interest points.
///
#include <vw/InterestPoint/Matcher.h>

namespace vw {
namespace ip {

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

}} // namespace vw::ip
