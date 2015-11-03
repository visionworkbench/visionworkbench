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


#ifndef __VW_MATH_FLANNTREE_H__
#define __VW_MATH_FLANNTREE_H__

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Core/Log.h>

#include <stddef.h>

#include <boost/noncopyable.hpp>

namespace vw {
namespace math {

  enum FLANN_DistType {FLANN_DistType_Unsupported, 
                       FLANN_DistType_L2, 
                       FLANN_DistType_Hamming};


  /// FLANN = FLANN is a library for performing fast approximate nearest
  ///         neighbor searches in high dimensional spaces.
  /// - Currently supports T = float, double, or unsigned char.
  /// - Make sure that the input features match the requested distance type.
  /// - Currently the real types support L2 and unsigned char supports Hamming.
  template <class T>
  class FLANNTree : boost::noncopyable {

  private:

    size_t m_num_features_loaded;
    void* m_index_ptr;
    FLANN_DistType m_dist_type;
    Matrix<T> m_features_cast; // The index makes pointers to this object. So we copy it.

    /// Returns the number of results found (usually knn)
    size_t knn_search_help( void* data_ptr, size_t rows, size_t cols, // Values we are looking for
                            Vector<int   >& indices,  // Index of each result
                            Vector<double>& dists,    // Distance of each result
                            size_t knn );             // Number of results to return

    /// Make a FLANN index wrapping a matrix of feature data
    void construct_index( void* data_ptr, size_t rows, size_t cols );

  public: // Functions

    /// Simple constructor.  Call load_match_data() before calling knn_search()!
    FLANNTree() 
      : m_num_features_loaded(0), m_index_ptr(NULL), m_dist_type(FLANN_DistType_Unsupported) {}

    /// Destructor
    ~FLANNTree();

    /// Load the matrix of feature data to match to and set the match distance type.
    /// - It is up to the user to make sure that the features matrix matches the distance type being used!
    /// - Hamming distance can take any data type but just compares the bits of the input data.
    template <class MatrixT>
    void load_match_data( MatrixBase<MatrixT> const& features,  FLANN_DistType dist_type) {
      if (features.impl().rows() == 0)
        vw_throw( ArgumentErr() << "Cannot create a FLANN tree with no input data!" );
      m_dist_type           = dist_type;
      m_features_cast       = features;
      m_num_features_loaded = m_features_cast.rows();
      construct_index( (void*)&m_features_cast(0,0), m_features_cast.rows(), m_features_cast.cols() );

      //vw_out() << "Loading " << m_features_cast.rows() << " FLANN targets of size " 
      //         << m_features_cast.cols() << "\n";
    }

    /// Multiple query access via VW's Matrix
    template <class MatrixT>
    size_t knn_search( MatrixBase<MatrixT> const& query,  // Values we are looking for
                       Vector<int   >& indices,           // Index of each result
                       Vector<double>& dists,             // Distance of each result
                       size_t knn ) {                     // Number of results to return

      // Convert query to our type
      Matrix<T> query_cast = query.impl();
      int flann_found = knn_search_help( (void*)&query_cast(0,0), query_cast.rows(), query_cast.cols(),
                                         indices, dists, knn );
      // Double check the number of valid indices found
      size_t num_found = 0;
      for (int i=0; i<flann_found; ++i)
        if ( (indices[i] >= 0) && (indices[i] < static_cast<int>(m_num_features_loaded)) )
          ++num_found;
      return num_found;
    }

    /// Single query access via ASP's Matrix
    template <class VectorT>
    size_t knn_search( VectorBase<VectorT> const& query, // Values we are looking for
                       Vector<int   >& indices,          // Index of each result
                       Vector<double>& dists,            // Distance of each result
                       size_t knn ) {                    // Number of results to return

      // Convert query to our type
      Vector<T> query_cast = query.impl();
      int flann_found = knn_search_help( (void*)&query_cast[0], 1, query_cast.size(), indices, dists, knn );

      // Count the number of valid indices obtained
      size_t num_found = 0;
      for (int i=0; i<flann_found; ++i)
        if ( (indices[i] >= 0) && (indices[i] < static_cast<int>(m_num_features_loaded)) )
          ++num_found;
      return num_found;
    }

    size_t size1() const;
    size_t size2() const;

  };

}} // end namespace vw::math

#endif//__VW_MATH_FLANNTREE_H__
