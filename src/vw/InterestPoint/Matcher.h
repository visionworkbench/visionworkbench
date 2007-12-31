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

#include <algorithm>

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

/// L2 Norm: returns Euclidian distance squared between pair of interest point descriptors.  Optional 
/// argument "maxdist" to provide early termination of computation if result will exceed maxdist.
class L2NormMetric {
public:
  float operator() (InterestPoint const& ip1, InterestPoint const& ip2, float maxdist = DBL_MAX) const {
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
  float operator() (InterestPoint const& ip1, InterestPoint const& ip2, float maxdist = DBL_MAX) const {
    float dist = 0.0;
    for (unsigned int i = 0; i < ip1.descriptor.size(); i++) {
      //dist += ip1.descriptor[i] * log2f(ip1.descriptor[i]/ip2.descriptor[i]) ;
      // log2(x) = log(x)/log(2);  MSVC 8 does not have log2
      dist += ip1.descriptor[i] * logf(ip1.descriptor[i]/ip2.descriptor[i])/logf(2.) ;
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
  template <class ListT>
  std::vector<typename ListT::const_iterator> operator()(InterestPoint const& /*ip*/, ListT const& candidates) const {
    return std::vector<typename ListT::const_iterator>(1, candidates.end());  // signals no constraint.
  }

  template <class ListT>
  std::vector<typename ListT::const_iterator> inverse(InterestPoint const& /*ip*/, ListT const& candidates) const {
    return std::vector<typename ListT::const_iterator>(1, candidates.end());  // signals no constraint.
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
  
  template <class ListT>
  std::vector<typename ListT::const_iterator> operator()(InterestPoint const& ip, ListT const& candidates ) const {
    typedef typename ListT::const_iterator IterT;

    std::vector<IterT> result;
    double sr, od;
    for (IterT i = candidates.begin(); i != candidates.end(); ++i) {
      sr = (*i).scale / ip.scale;
      od = (*i).orientation - ip.orientation;
      // Bring orientation delta (od) into range -M_PI to M_PI
      if (od < -M_PI) od += M_PI*2;
      else if (od > M_PI) od -= M_PI*2;
      
      if (sr >= scale_ratio_min && sr <= scale_ratio_max &&
          od >= ori_diff_min && od <= ori_diff_max) {
        result.push_back(i);
      }
    }
    return result;
  }

  template <class ListT>
  std::vector<typename ListT::const_iterator> inverse(InterestPoint const& ip, ListT const& candidates ) const {
    typedef typename ListT::const_iterator IterT;

    std::vector<IterT> result;
    double sr, od;
    for (IterT i = candidates.begin(); i != candidates.end(); ++i) {
      sr = ip.scale / (*i).scale;
      od = ip.orientation - (*i).orientation;      
      
      if (sr >= scale_ratio_min && sr <= scale_ratio_max &&
          od >= ori_diff_min && od <= ori_diff_max) {
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
  
  template <class ListT>
  std::vector<typename ListT::const_iterator> operator()(InterestPoint const& ip, ListT const& candidates ) const {
    typedef typename ListT::const_iterator IterT;

    std::vector<IterT> result;
    double dx, dy;
    for (IterT i = candidates.begin(); i != candidates.end(); ++i) {
      dx = (*i).x - ip.x;
      dy = (*i).y - ip.y;
      if (dx >= min_x && dx < max_x && dy >= min_y && dy < max_y) {
        result.push_back(i);
      }
    }
    return result;
  }

  template <class ListT>
  std::vector<typename ListT::const_iterator> inverse(InterestPoint const& ip, ListT const& candidates ) const {
    typedef typename ListT::const_iterator IterT;

    std::vector<IterT> result;
    double dx, dy;
    for (IterT i = candidates.begin(); i != candidates.end(); ++i) {
      dx = ip.x - (*i).x;
      dy = ip.y - (*i).y;
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
  
  InterestPointMatcher(double threshold = 0.5, MetricT metric = MetricT(), ConstraintT constraint = ConstraintT()) 
    : constraint(constraint), distance_metric(metric), threshold(threshold) { }
  
  /// Given two lists of interest points, this routine returns the two lists
  /// of matching interest points based on the Metric and Constraints
  /// provided by the user.
  template <class ListT, class MatchListT>
  void match( ListT const& ip1, ListT const& ip2,
              MatchListT& matched_ip1, MatchListT& matched_ip2,
              bool bidirectional = false) const {
    typedef typename ListT::const_iterator IterT;

    std::vector<IterT> index;
    if (bidirectional) 
      index = bimatch_index(ip1,ip2);
    else 
      index = match_index(ip1,ip2);
    
    matched_ip1.clear(); matched_ip2.clear();
    IterT ip1_iter = ip1.begin();
    for (unsigned int i = 0; i < index.size(); ++i, ++ip1_iter) {
      if ( index[i] != ip2.end() ) {
        matched_ip1.push_back(*ip1_iter);
        matched_ip2.push_back(*(index[i]));
      }
    }
  }
    
private:
  /// Return vector of indices of matches from each point in ip1 to best match (if any) in ip2.
  /// Usage: indx = matchindx(ip1,ip2) means that 
  ///     ip1[i] has no match in ip2  indx v[i] == -1
  ///     otherwise, ip2[indx[i]] is best match to ip[i]. 
  template <class ListT>
  std::vector<typename ListT::const_iterator> match_index(ListT const& ip1, ListT const& ip2) const {
    typedef typename ListT::const_iterator IterT;

    std::vector<IterT> indx(ip1.size(), ip2.end());
    int progress_inc = std::max((int)(ip1.size() / 10), 1);
    vw_out(InfoMessage, "interest_point") << "\tFinding interest point matches" << std::flush;
    int n = 0;
    for (IterT ip1_i = ip1.begin(); ip1_i != ip1.end(); ++ip1_i, ++n) {
      if (n % progress_inc == 0)
        vw_out(InfoMessage, "interest_point") << "." << std::flush;
      std::vector<IterT> candidates = constraint(*ip1_i, ip2);
      indx[n] = best_match(*ip1_i, ip2, candidates);
    }
    vw_out(InfoMessage, "interest_point") << " done.\n";
    
    return indx;
  }
  
  /// Return vector of indices of matches from each point in ip1 to best match (if any) in ip2, using bidirectional
  /// search (ie matches are one to one).
  /// Usage: indx = matchindx(ip1,ip2) means that 
  ///     ip1[i] has no match in ip2  indx v[i] == -1
  ///     otherwise, ip2[indx[i]] is best match to ip[i] AND ip[i] is best match to ip2[indx[i]].
  template <class ListT>
  std::vector<typename ListT::const_iterator> bimatch_index(ListT const& ip1, ListT const& ip2) const {
    typedef typename ListT::const_iterator IterT;

    std::vector<IterT> indx12 = match_index(ip1,ip2);
    // FIXME: use pthreads to create separate thread for this, to run in
    // parallel with computation of indx21 (both readonly access to ip1, ip2).
    
    //std::vector<int> indx21 = match_index(ip2, ip1);
    std::vector<IterT> indx21(ip2.size(), ip1.end());       // Cannot repeat call to matchindx, since need to use inverse constraints.
    int n = 0;
    for (IterT ip2_i = ip2.begin(); ip2_i != ip2.end(); ++ip2_i, ++n) {
      std::vector<IterT> candidates = constraint.inverse(*ip2_i, ip1);
      indx21[n] = best_match(*ip2_i, ip1, candidates);
    }

    // TODO: Well, I broke this. Save n's in 2nd loop above? (PM)
    /*
    for (unsigned int i = 0; i < indx12.size(); ++i)
      if ( indx12[i] == ip2.end() || indx21[indx12[i]] == ip1.end() || indx21[indx12[i]] != i )
        indx12[i] = ip2.end();
    */

    return indx12;
  };

  // Returns the index of the best match of 'query' in 'interest_points'.  
  // If candidate_indices is provided, only the subset of candidates will 
  // be compared.  Returns the index of the best match or -1 if no suitable 
  // match can be found.
  template <class ListT>
  typename ListT::const_iterator
  best_match(InterestPoint const& query, ListT const& interest_points, 
             std::vector<typename ListT::const_iterator> const&
             candidate_indices = std::vector<typename ListT::const_iterator>()) const {
    typedef typename ListT::const_iterator IterT;
    float dist = 0.0, dist0 = FLT_MAX, dist1 = FLT_MAX;
    IterT best_match = interest_points.end();
    
    // An empty candidates vector indicates that we should search all of interest_points.
    if (candidate_indices.empty()) return interest_points.end();
    if (candidate_indices[0] == interest_points.end()) {

      // find closest in FULL SET of interest_points
      for (IterT i = interest_points.begin(); i != interest_points.end(); ++i) {
        dist = distance_metric(query, *i, dist1);
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
        dist = distance_metric(query, *(candidate_indices[i]), dist1);
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
      return interest_points.end();
  }

};
  
//////////////////////////////////////////////////////////////////////////////
// Convenience Typedefs
typedef InterestPointMatcher< L2NormMetric, NullConstraint > DefaultMatcher;
typedef InterestPointMatcher< L2NormMetric, ScaleOrientationConstraint > ConstraintedMatcher;

}} // namespace vw::ip

#endif // _INTEREST_POINT_MATCHER_H_
