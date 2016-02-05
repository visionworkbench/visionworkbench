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

#ifndef __VW_MATH_RANSAC_H__
#define __VW_MATH_RANSAC_H__

#include <vw/Math/Vector.h>
#include <vw/Core/Log.h>

namespace vw {
namespace math {

  VW_DEFINE_EXCEPTION(RANSACErr, Exception);

  /// This is a basic error metric can be used when the mathematical
  /// multiplication and subtraction operators are defined for p1, p2, and H.
  struct L2NormErrorMetric {
    template <class RelationT, class ContainerT>
    double operator() (RelationT  const& H,
                       ContainerT const& p1,
                       ContainerT const& p2) const {
      return vw::math::norm_2( p2 - H * p1 );
    }
  };

  template <size_t dim>
  struct HomogeneousL2NormErrorMetric {
    template <class RelationT, class ContainerT>
    double operator() (RelationT  const& H,
                       ContainerT const& p1,
                       ContainerT const& p2) const {
      // Copy the input data into homogenous vector objects
      Vector<double, dim+1> p1_h, p2_h;
      for (unsigned i = 0; i < dim; i++) {
        p1_h[i] = p1[i];
        p2_h[i] = p2[i];
      }
      p1_h[dim] = 1;
      p2_h[dim] = 1;

      Vector<double, dim+1> inter_result = H*p1_h;
      // Re-normalizing. This conditional should only throw if H is
      // an homography matrix
      if ( inter_result[dim] != 1 )
        inter_result /= inter_result[dim];

      return vw::math::norm_2(p2_h - inter_result);
    }
  };

  /// This metric can be used to measure the error between an interest
  /// point p2 and a second interest point p1 that is transformed by a
  /// 3x3 matrix H.  This is predominately used when matching interest
  /// points using RANSAC.
  typedef HomogeneousL2NormErrorMetric<2> InterestPointErrorMetric;

  /// RANSAC Driver class
  template <class FittingFuncT, class ErrorFuncT>
  class RandomSampleConsensus {
    const FittingFuncT& m_fitting_func;
    const ErrorFuncT  & m_error_func;
          int           m_num_iterations;
          double        m_inlier_threshold;
          int           m_min_num_output_inliers;
          bool          m_reduce_min_num_output_inliers_if_no_fit;
    
    /// \cond INTERNAL
    // Utility Function: Pick N UNIQUE, random integers in the range [0, size]
    inline void get_n_unique_integers(int size, std::vector<int> & samples) const {

      // Note: We do not modify the initial random seed. As such, if
      // a program uses RANSAC, repeatedly running this program will
      // always return the same results. However, if that program
      // calls RANSAC twice while within the same instance of the
      // program, the second time the result of RANSAC will be
      // different, since we keep on pulling new random numbers.
        
      int n = samples.size();
      VW_ASSERT(size >= n, ArgumentErr() << "Not enough samples (" << n << " / " << size << ")\n");

      const double divisor = static_cast<double>(RAND_MAX) + 1.0;
      for (int i = 0; i < n; ++i) {
        bool done = false;
        while (!done) {
          samples[i] = static_cast<int>( (static_cast<double>(std::rand()) / divisor) * size );
          done = true;
          for (int j = 0; j < i; j++)
            if (samples[i] == samples[j])
              done = false;
        }
      }
    }
    /// \endcond

  public:

    // Returns the list of inliers.
    template <class ContainerT1, class ContainerT2>
    void inliers(typename FittingFuncT::result_type const& H,
                          std::vector<ContainerT1>  const& p1, 
                          std::vector<ContainerT2>  const& p2,
                          std::vector<ContainerT1>       & inliers1, 
                          std::vector<ContainerT2>       & inliers2) const {

      inliers1.clear();
      inliers2.clear();

      for (size_t i=0; i<p1.size(); i++) {
        if (m_error_func(H,p1[i],p2[i]) < m_inlier_threshold) {
          inliers1.push_back(p1[i]);
          inliers2.push_back(p2[i]);
        }
      }
    }

    // Returns the list of inlier indices.
    template <class ContainerT1, class ContainerT2>
    std::vector<size_t> inlier_indices(typename FittingFuncT::result_type const& H,
                                                std::vector<ContainerT1>  const& p1,
                                                std::vector<ContainerT2>  const& p2) const {
      std::vector<size_t> result;
      for (size_t i=0; i<p1.size(); i++)
        if (m_error_func(H,p1[i],p2[i]) < m_inlier_threshold)
          result.push_back(i);
      return result;
    }

    void reduce_min_num_output_inliers(){
      m_min_num_output_inliers = int(m_min_num_output_inliers/1.5);
    }
      
    /// Constructor - Stores all the inputs in member variables
    RandomSampleConsensus(FittingFuncT const& fitting_func, 
                          ErrorFuncT   const& error_func,
                          int    num_iterations,
                          double inlier_threshold,
                          int    min_num_output_inliers,
                          bool   reduce_min_num_output_inliers_if_no_fit = false
                          ):
      m_fitting_func(fitting_func), m_error_func(error_func),
      m_num_iterations(num_iterations), 
      m_inlier_threshold(inlier_threshold),
      m_min_num_output_inliers(min_num_output_inliers),
      m_reduce_min_num_output_inliers_if_no_fit(reduce_min_num_output_inliers_if_no_fit){}

    /// As attempt_ransac but keep trying with smaller numbers of required inliers.
    template <class ContainerT1, class ContainerT2>
    typename FittingFuncT::result_type operator()(std::vector<ContainerT1> const& p1,
                                                  std::vector<ContainerT2> const& p2) {

      // Try to fit using RANSAC. Perform repeated fits with smaller
      // m_min_num_output_inliers if the fit fails and
      // m_reduce_min_num_output_inliers_if_no_fit is true.

      typename FittingFuncT::result_type H;
      bool success = false;
      
      for (int attempt = 0; attempt < 10; attempt++){
        try{
          H = attempt_ransac(p1, p2);
          success = true;
          break; 
        } catch ( const std::exception& e ) { 
          vw_out() << e.what() << "\n";
          if (!m_reduce_min_num_output_inliers_if_no_fit) 
            break;
          reduce_min_num_output_inliers();
          if (m_min_num_output_inliers < 2) // Can't possibly compute a transform with 1 or 0 samples!
            break;
          vw_out() << "Attempting RANSAC with " << m_min_num_output_inliers << " number of output inliers.\n";
          
        }
      }

      if (!success) 
        vw_throw( RANSACErr() << "RANSAC was unable to find a fit that matched the supplied data." );

      return H;
    }

    /// Run RANSAC on two input data lists using the current parameters.
    template <class ContainerT1, class ContainerT2>
    typename FittingFuncT::result_type attempt_ransac(std::vector<ContainerT1> const& p1,
                                                      std::vector<ContainerT2> const& p2) const {

      VW_ASSERT( !p1.empty(),
                 RANSACErr() << "RANSAC Error.  Insufficient data.\n");
      VW_ASSERT( p1.size() == p2.size(),
                 RANSACErr() << "RANSAC Error.  Data vectors are not the same size." );

      int min_elems_for_fit = m_fitting_func.min_elements_needed_for_fit(p1[0]);

      VW_ASSERT( (int)p1.size() >= min_elems_for_fit,
                 RANSACErr() << "RANSAC Error.  Not enough potential matches for this fitting functor. (" << p1.size() << "/" << min_elems_for_fit << ")\n");

      VW_ASSERT( m_min_num_output_inliers >= min_elems_for_fit,
                 RANSACErr() << "RANSAC Error.  Number of requested inliers is less than min number of elements needed for fit. (" << m_min_num_output_inliers << "/" << min_elems_for_fit << ")\n");

      typename FittingFuncT::result_type best_H;

      std::vector<ContainerT1> try1;
      std::vector<ContainerT2> try2;
      std::vector<int> random_indices(min_elems_for_fit);

      int num_inliers = 0;
      double min_err = std::numeric_limits<double>::max();
      for (int iteration = 0; iteration < m_num_iterations; ++iteration) {

        // 0. Get min_elems_for_fit points at random, taking care not
        //    to select the same point twice.
        get_n_unique_integers(p1.size(), random_indices);
        // Resizing below is essential, as by now their size may have changed
        try1.resize(min_elems_for_fit);
        try2.resize(min_elems_for_fit);
        for (int i = 0; i < min_elems_for_fit; ++i) {
          try1[i] = p1[random_indices[i]];
          try2[i] = p2[random_indices[i]];
        }

        // 1. Compute the fit using these samples.
        typename FittingFuncT::result_type H = m_fitting_func(try1, try2);

        // 2. Find all the inliers for this fit.
        inliers(H, p1, p2, try1, try2);

        // 3. Skip this model if too few inliers.
        if ((int)try1.size() < m_min_num_output_inliers) 
          continue;

        // 4. Re-estimate the model using the inliers.
        H = m_fitting_func(try1, try2, H);
        
        // 5. Find the mean error for the inliers.
        double err_val = 0.0;
        for (size_t i = 0; i < try1.size(); i++) 
          err_val += m_error_func(H, try1[i], try2[i]);
        err_val /= try1.size();

        // 6. Save this model if its error is lowest so far.
        if (err_val < min_err){
          min_err     = err_val;
          best_H      = H;
          num_inliers = try1.size();
        }

      }

      if (num_inliers < m_min_num_output_inliers) {
        vw_throw( RANSACErr() << "RANSAC was unable to find a fit that matched the supplied data." );
      }

      // For debugging
      VW_OUT(InfoMessage, "interest_point") << "\nRANSAC Summary:"     << std::endl;
      VW_OUT(InfoMessage, "interest_point") << "\tFit = "              << best_H      << std::endl;
      VW_OUT(InfoMessage, "interest_point") << "\tInliers / Total  = " << num_inliers << " / " << p1.size() << "\n\n";
      
      return best_H;
    }

  }; // End of RandomSampleConsensus class definition

  // Helper function to instantiate a RANSAC class object and immediately call it
  template <class ContainerT1, class ContainerT2, class FittingFuncT, class ErrorFuncT>
  typename FittingFuncT::result_type ransac(std::vector<ContainerT1> const& p1,
                                            std::vector<ContainerT2> const& p2,
                                            FittingFuncT             const& fitting_func,
                                            ErrorFuncT               const& error_func,
                                            int     num_iterations,
                                            double  inlier_threshold,
                                            int     min_num_output_inliers,
                                            bool    reduce_min_num_output_inliers_if_no_fit = false
                                            ) {
    RandomSampleConsensus<FittingFuncT, ErrorFuncT> ransac_instance(fitting_func,
                                                                    error_func,
                                                                    num_iterations,
                                                                    inlier_threshold,
                                                                    min_num_output_inliers,
                                                                    reduce_min_num_output_inliers_if_no_fit
                                                                    );
    return ransac_instance(p1,p2);
  }

}} // namespace vw::math

#endif // __MATH_RANSAC_H__
