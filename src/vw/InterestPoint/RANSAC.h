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

/// \file RANSAC.h
/// 
/// Robust outlier rejection using the Random Sample Consensus
/// (RANSAC) algorithm.  
/// 
/// This technique was introduced in: 
///
/// Fischler, Martin A. and Bolles, Robert C. "Random Sample
/// Consensus: a Paradigm for Model Fitting with Applications to
/// Image Analysis and Automated Cartography" (1981)
///  
/// This generic implementation of RANSAC requires that you supply two
/// functors:
///
/// 1. A "Fitting" functor that, given two vectors of putative
///    matches, returns an object that is the best fit between the two
///    data sets.  For example, if you are attempting to robustly
///    compute the best homography that best describes the
///    transformation between two sets of N-vectors, the result of the
///    fitting function would be an NxN matrix.
/// 
///    The fitting function must also provide a result_type typedef
///    and a method called min_elements_needed_for_fit(example), which
///    returns the minimum number of matching data points required to
///    compute a fit.  Again, as an example, to fit a homography
///    between two sets of 2-D points, you need at least 4 point pairs.
///
/// 2. A "Error" functor that, given a pair of data p1 and p2, and
///    the computed fit function H, returns the distance between p2
///    and H(p1).  You can think of this as the error for a given p1,
///    p2 and H.  In the case where you are fitting a homography to a
///    set of points, this routine could compute the 2-norm of the
///    error: || p2 - H * p1 ||
///

#ifndef __VW_INTERESTPOINT_RANSAC_H__
#define __VW_INTERESTPOINT_RANSAC_H__

#include <vw/Math/Vector.h>

namespace vw { 
namespace ip {

  VW_DEFINE_EXCEPTION(RANSACErr, Exception);

  /// This is a basic error metric can be used when the mathematical
  /// multiplication and subtraction operators are defined for p1, p2,
  /// and H.
  struct L2NormErrorMetric {
    template <class RelationT, class ContainerT>
    double operator() (RelationT const& H,
                       ContainerT const& p1, 
                       ContainerT const& p2) const {
      return vw::math::norm_2( p2 - H * p1 );
    }
  };

  /// This metric can be used to measure the error between a interest
  /// point p2 and a second interest point p1 that is transformed by a
  /// 3x3 matrix H.  This is predominately used when matching interest
  /// points using RANSAC.
  struct InterestPointErrorMetric {
    template <class RelationT, class ContainerT>
    double operator() (RelationT const& H,
                       ContainerT const& p1, 
                       ContainerT const& p2) const {
      return vw::math::norm_2( Vector3(p2.x,p2.y,1) - H * Vector3(p1.x,p1.y,1));
    }
  };

  /// RANSAC Driver class
  template <class FittingFuncT, class ErrorFuncT>
  class RandomSampleConsensus {
    FittingFuncT m_fitting_func;
    ErrorFuncT m_error_func;
    double m_inlier_threshold;

    // Returns the number of inliers for a given threshold.
    template <class ContainerT>
    unsigned num_inliers(typename FittingFuncT::result_type const& H,
                    std::vector<ContainerT> const& p1, std::vector<ContainerT> const& p2) const {
        
      unsigned result = 0;
      for (unsigned i=0; i<p1.size(); i++) {
        if (m_error_func(H,p1[i],p2[i]) < m_inlier_threshold) 
          ++result;
      }
      return result;
    }
    
    /// \cond INTERNAL  
    // Utility Function: Pick N UNIQUE, random integers in the range [0, size] 
    inline void _vw_get_n_unique_integers(unsigned int size, unsigned n, int* samples) {
      VW_ASSERT(size >= n, ArgumentErr() << "Not enough samples (" << n << " / " << size << ")\n");
      
      for (unsigned i=0; i<n; ++i) {
        bool done = false;
        while (!done) {
          samples[i] = (int)(((double)random() / (double)RAND_MAX) * size);
          done = true;
          for (unsigned j = 0; j < i; j++) 
            if (samples[i] == samples[j]) 
              done = false;
        }
      }
    }
    /// \endcond
    
  public:

    // Returns the list of inlier indices.
    template <class ContainerT>
    void inliers(typename FittingFuncT::result_type const& H,
                 std::vector<ContainerT> const& p1,std::vector<ContainerT> const& p2,
                 std::vector<ContainerT> &inliers1, std::vector<ContainerT> &inliers2) const {
        
      inliers1.clear();
      inliers2.clear();
      
      for (unsigned int i=0; i<p1.size(); i++) {
        if (m_error_func(H,p1[i],p2[i]) < m_inlier_threshold) {
          inliers1.push_back(p1[i]);
          inliers2.push_back(p2[i]);
        }
      }
    }

    RandomSampleConsensus(FittingFuncT const& fitting_func, ErrorFuncT const& error_func, double inlier_threshold)
      : m_fitting_func(fitting_func), m_error_func(error_func), m_inlier_threshold(inlier_threshold) {}
    
    template <class ContainerT>
    typename FittingFuncT::result_type operator()(std::vector<ContainerT> const& p1, 
                                                  std::vector<ContainerT> const& p2,
                                                  int ransac_iterations = 0) {
      // check consistency
      VW_ASSERT( p1.size() == p2.size(), 
                 RANSACErr() << "RANSAC Error.  data vectors are not the same size." );
      VW_ASSERT( p1.size() != 0,  
                 RANSACErr() << "RANSAC Error.  Insufficient data.\n");
      VW_ASSERT( p1.size() >= m_fitting_func.min_elements_needed_for_fit(p1[0]),  
                 RANSACErr() << "RANSAC Error.  Not enough potential matches for this fitting funtor. ("<<p1.size() << "/" << m_fitting_func.min_elements_needed_for_fit(p1[0]) << ")\n");

      unsigned inliers_max = 0;
      typename FittingFuncT::result_type H;
      typename FittingFuncT::result_type H_max;

      /////////////////////////////////////////
      // First part: 
      //   1. choose N points at random
      //   2. find a fit for those N points
      //   3. check for consensus
      //   4. keep fit with best consensus so far
      /////////////////////////////////////////

      // Seed random number generator
      srandom((unsigned int) clock());

      // This is a rough value, but it seems to produce reasonably good results.
      if (ransac_iterations == 0) 
        ransac_iterations = p1.size() * 2;

      int n = m_fitting_func.min_elements_needed_for_fit(p1[0]);
      std::vector<ContainerT> try1(n), try2(n);
      int random_indices[n];
      for (int iteration=0; iteration < ransac_iterations; ++iteration) {
        // Get four points at random, taking care not 
        // to select the same point twice.
        _vw_get_n_unique_integers(p1.size(), n, random_indices);

        for (int i=0; i < n; ++i) {
          try1[i] = p1[random_indices[i]];
          try2[i] = p2[random_indices[i]];
        }
        
        // Compute the fit using these samples
        H = m_fitting_func(try1, try2);
        
        // Compute consensuss
        unsigned n_inliers = num_inliers(H, p1, p2);

        // Keep best consensus
        if (n_inliers > inliers_max) {
          inliers_max = n_inliers;
          H_max = H;
        }
      }
      
      if (inliers_max < m_fitting_func.min_elements_needed_for_fit(p1[0])) {
        vw_throw( RANSACErr() << "RANSAC was unable to find a fit that matched the supplied data." );
      }

      // For debugging
      vw_out(InfoMessage, "interest_point") << "\nBest overall:" << std::endl;
      vw_out(InfoMessage, "interest_point") << "\tFit = " << H_max << std::endl;
      vw_out(InfoMessage, "interest_point") << "\tInliers / Total  = " << inliers_max << " / " << p1.size() << "\n\n";

      ////////////////////////////////////
      // Second part:
      //    1. find all inliers the best fit
      //    2. re-estimate the fit using all inliers
      //    3. repeat until # of inliers stabilizes
      ///////////////////////////////////
      unsigned int num_old = 0;
      inliers(H_max, p1, p2, try1, try2 );
      while( try1.size() > num_old ){
        num_old = try1.size();
        H = m_fitting_func(try1, try2);
        inliers( H, p1, p2, try1, try2 );
      }
      return H;
    }
    
  };

  // Free function
  template <class ContainerT, class FittingFuncT, class ErrorFuncT>
  typename FittingFuncT::result_type ransac(std::vector<ContainerT> const& p1, 
                                            std::vector<ContainerT> const& p2,
                                            FittingFuncT const& fitting_func, 
                                            ErrorFuncT const& error_func,
                                            double inlier_threshold = 50) {
    RandomSampleConsensus<FittingFuncT, ErrorFuncT> ransac_instance(fitting_func, 
                                                                    error_func, 
                                                                    inlier_threshold);
    return ransac_instance(p1,p2);
  }

}} // namespace vw::ip

#endif // __INTERESTPOINT_RANSAC_H__
