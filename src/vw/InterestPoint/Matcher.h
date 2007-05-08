// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

/// \file Matcher.h
/// 
/// Classes and functions for matching image interest points.
/// 
#ifndef _INTERESTPOINT_MATCHER_H_
#define _INTERESTPOINT_MATCHER_H_

#include <vw/InterestPoint/Descriptor.h>
#include <vector>
#include <math.h>
#include <float.h>

namespace vw { 
namespace ip {

//////////////////////////////////////////////////////////////////////////////////////////////////
// Interest point metrics
//  Must be of form: 
//      float operator() (const InterestPoint& ip1, const InterestPoint& ip2, float maxdist = DBL_MAX) 

/// L2 Norm: returns Euclidian distance squared between pair of keypoint descriptors.  Optional 
/// argument "maxdist" to provide early termination of computation if result will exceed maxdist.
class L2NormMetric {
public:
  float operator() (const InterestPoint& ip1, const InterestPoint& ip2, float maxdist = DBL_MAX) const {
    float dist = 0.0;
    for (unsigned int i = 0; i < ip1.descriptor.size(); i++) {
      dist += pow((ip1.descriptor[i] - ip2.descriptor[i]),2);
      if (dist > maxdist) break;  // abort calculation if distance exceeds upper bound
    }
    return dist;
  }
};

/// KL distance (relative entropy) between interest point descriptors. Optional
/// argument "maxdist" to provide early termination of computation if result will exceed maxdist.
class RelativeEntropyMetric {
public:
  float operator() (const InterestPoint& ip1, const InterestPoint& ip2, float maxdist = DBL_MAX) const {
    float dist = 0.0;
    for (unsigned int i = 0; i < ip1.descriptor.size(); i++) {
      dist += ip1.descriptor[i] * log2f(ip1.descriptor[i]/ip2.descriptor[i]) ;
      if (dist > maxdist) break;  // abort calculation if distance exceeds upper bound
    }
    return dist;
  }
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Interest Point Match contraints functors to return a list of
// allowed match candidates to an interest point.
// 
// Must have be of form:
//     std::vector<int> operator()( const InterestPoint& ip, const std::vector<InterestPoint>& candidates)
// returning list of points in candidates that are sufficiently close to ip.
//
// To use bimatch, must also include the inverse function:
//     std::vector<int> inverse ()( const InterestPoint& ip, const std::vector<InterestPoint>& candidates)
//  returning list of points in candidates that are close to ip, under the inverse constraint mapping.
//
//  0) Default: no constraints.  Return a list with -1 as first element.
//  i) The ratio of scales between two interest points and difference in orientations
//  ii) Geometrical constraints on point locations, such as epipolar constraints, homography's etc.
// be set.

/// default matching constraint functor providing NO constraints on putative matches.  Returns 
/// vector<int> = {-1}, interpreted by IP_Best_Match as indicating there are no constraints (as
/// opposed to an empty vector, which indicates overly tight constraints.
class NullConstraint {
public:
  std::vector<int> operator()( const InterestPoint& ip, const std::vector<InterestPoint>& candidates) const {
    return std::vector<int>(1,-1);  // signals no constraint.
  }
  std::vector<int> inverse( const InterestPoint& ip, const std::vector<InterestPoint>& candidates) const {
    return std::vector<int>(1,-1);  // signals no constraint.
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
  
  std::vector<int> operator()(const InterestPoint& ip, const std::vector<InterestPoint>& candidates ) const {

    std::vector<int> result;
    double sr, od;
    for (unsigned int i = 0; i < candidates.size(); i++) {
      sr = candidates[i].scale/ip.scale;
      od = candidates[i].orientation - ip.orientation;
      
      if (sr >= scale_ratio_min && sr < scale_ratio_max &&
        od >= ori_diff_min && od < ori_diff_max) {
        result.push_back(i);
      }
    }
    return result;
  }

  std::vector<int> inverse(const InterestPoint& ip, const std::vector<InterestPoint>& candidates ) const {
    std::vector<int> result;
    double sr, od;
    for (unsigned int i = 0; i < candidates.size(); i++) {
      sr = ip.scale/candidates[i].scale;
      od = ip.orientation - candidates[i].orientation;      
      
      if (sr >= scale_ratio_min && sr < scale_ratio_max &&
        od >= ori_diff_min && od < ori_diff_max) {
        result.push_back(i);
      }
    }
    return result;
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
                     double _min_y = -10.0, double _max_y = 10.0)
  {
    min_x = _min_x; max_x = _max_x; min_y = _min_y; max_y = _max_y;
  }
  
  std::vector<int> operator()(const InterestPoint& ip, const std::vector<InterestPoint>& candidates ) const {
    std::vector<int> result;
    double dx, dy;
    for (unsigned int i = 0; i < candidates.size(); i++) {
      dx = candidates[i].x - ip.x;
      dx = candidates[i].y - ip.y;
      if (dx >= min_x && dx < max_x && dy >= min_y && dy < max_y) {
        result.push_back(i);
      }
    }
    return result;
  }
  
  std::vector<int> inverse(const InterestPoint& ip, const std::vector<InterestPoint>& candidates ) const {
    std::vector<int> result;
    double dx, dy;
    for (unsigned int i = 0; i < candidates.size(); i++) {
      dx = ip.x - candidates[i].x;
      dx = ip.y - candidates[i].y;
      if (dx >= min_x && dx < max_x && dy >= min_y && dy < max_y) {
        result.push_back(i);
      }
    }
    return result;
  }
  
};

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Interest Point Matcher class
template < class MetricT, class ConstraintT >
class InterestPointMatcher {
public:
  
  ConstraintT constraint;
  MetricT distance_metric;
  double threshold;
  
  InterestPointMatcher(MetricT metric = MetricT(), ConstraintT constraint = ConstraintT(), double threshold = 0.4) 
    : distance_metric(metric), constraint(constraint), threshold(threshold) { }
  
  /// Given two lists of keypoints, this routine returns the two lists
  /// of matching keypoints based on the Metric and Constraints
  /// provided by the user.
  void match( std::vector<InterestPoint> const& ip1,
              std::vector<InterestPoint> const& ip2,
              std::vector<InterestPoint>& matched_ip1,
              std::vector<InterestPoint>& matched_ip2,
              bool bidirectional = false) const {

    std::vector<int> index;
    if (bidirectional) 
      index = bimatch_index(ip1,ip2);
    else 
      index = match_index(ip1,ip2);
    
    matched_ip1.resize(0); matched_ip2.resize(0);
    for (unsigned int i = 0; i < index.size(); i++) {
      if (index[i] >= 0) {
        matched_ip1.push_back(ip1[i]);
        matched_ip2.push_back(ip2[index[i]]);
      }
    }
  }
    
private:
  /// Return vector of indices of matches from each point in ip1 to best match (if any) in ip2.
  /// Usage: indx = matchindx(ip1,ip2) means that 
  ///     ip1[i] has no match in ip2  indx v[i] == -1
  ///     otherwise, ip2[indx[i]] is best match to ip[i]. 
  std::vector<int> match_index(std::vector<InterestPoint> const& ip1, std::vector<InterestPoint> const& ip2) const
  {
    std::vector<int> indx(ip1.size(),-1);
    int progress_inc = ip1.size() / 10;
    vw_out(InfoMessage) << "\tFinding keypoint matches" << std::flush;
    for (unsigned int i=0; i<ip1.size(); i++) {
      if (i % progress_inc == 0)
        vw_out(InfoMessage) << "." << std::flush;
      std::vector<int> candidates = constraint(ip1[i], ip2);
      indx[i] = best_match(ip1[i], ip2, candidates);
    }
    vw_out(InfoMessage) << " done.\n";
    
    return indx;
  }
  
  /// Return vector of indices of matches from each point in ip1 to best match (if any) in ip2, using bidirectional
  /// search (ie matches are one to one).
  /// Usage: indx = matchindx(ip1,ip2) means that 
  ///     ip1[i] has no match in ip2  indx v[i] == -1
  ///     otherwise, ip2[indx[i]] is best match to ip[i] AND ip[i] is best match to ip2[indx[i]].
  std::vector<int> bimatch_index(std::vector<InterestPoint> const& ip1, std::vector<InterestPoint> const& ip2) const
  {
    std::vector<int> indx12 = match_index(ip1,ip2);   // FIXME: use pthreads to create seperate thread for this, to run in
                                                    // parallel with computation of indx21 (both readonly access to ip1, ip2).
    
    //std::vector<int> indx21 = match_index(ip2, ip1);
    std::vector<int> indx21(ip2.size(),-1);       // Cannot repeat call to matchindx, since need to use inverse constraints.
    for (unsigned int i=0; i<ip2.size(); i++) 
    {
      std::vector<int> candidates = constraint.inverse(ip2[i], ip1);
      indx21[i] = best_match(ip2[i], ip1, candidates);
    }
    
    for (unsigned int i=0; i<indx12.size(); i++) 
      if ( indx12[i] < 0 || indx21[indx12[i]] < 0 || indx21[indx12[i]] != i)
        indx12[i] = -1; 
    
    return indx12;
  };

  // Returns the index of the best match of 'query' in 'interest_points'.  
  // If candidate_indices is provided, only the subset of candidates will 
  // be compared.  Returns the index of the best match or -1 if no suitable 
  // match can be found.
  int best_match(InterestPoint const& query, std::vector<InterestPoint> const& interest_points, 
                 std::vector<int> const& candidate_indices = std::vector<int>()) const {
    
    float dist = 0.0, dist0 = DBL_MAX, dist1 = DBL_MAX;
    int best_match = 0;
    
    // An empty candidates vector indicates that we should search all of interest_points.
    if (candidate_indices.empty()) return -1;
    if (candidate_indices[0] == -1) {

      // find closest in FULL SET of interest_points
      for (unsigned int i = 0; i < interest_points.size(); i++) {
        dist = distance_metric(query, interest_points[i], dist1);
        if (dist < dist0) {
          dist1 = dist0;
          dist0 = dist;
          best_match = i;
        }
        else if (dist < dist1) dist1 = dist;
      }     

    } else {

      // find closest in entire SUBSET of ip2
      for (unsigned int i = 0; i < candidate_indices.size(); i++) {
        dist = distance_metric(query, interest_points[candidate_indices[i]], dist1);
        if (dist < dist0) {
          dist1 = dist0;
          dist0 = dist;
          best_match = candidate_indices[i];
        }
        else if (dist < dist1) dist1 = dist;
      }     
    }
    
    // best match acceptance test
    if (dist0 < threshold * dist1)
      return best_match;
    else
      return -1;
  }

};
  
//////////////////////////////////////////////////////////////////////////////
// Convencience Typedefs
typedef InterestPointMatcher< L2NormMetric, NullConstraint > DefaultMatcher;
typedef InterestPointMatcher< L2NormMetric, ScaleOrientationConstraint > ConstraintedMatcher;

}} // namespace vw::ip

#endif // _INTEREST_POINT_MATCHER_H_
