// __BEGIN_LICENSE__
//  Copyright (c) 2006-2012, United States Government as represented by the
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
#include <boost/utility.hpp>

namespace vw {
namespace math {

  // This provides only access to L2 distance metric in the FLANN
  // tree. We would need to make new objects for each distance metric
  // if we want to avoid having the FLANN header show up
  // everywhere. (Their header has lots of warnings).
  template <class FloatT>
  class FLANNTree : boost::noncopyable {
    void* m_index_ptr;
    Matrix<FloatT> m_features_cast; // The index makes pointers to this object. So we copy it.

    void knn_search_help( void* data_ptr, size_t rows, size_t cols,
                          Vector<int>& indices,
                          Vector<FloatT>& dists,
                          size_t knn );

    void construct_index( void* data_ptr, size_t rows, size_t cols );

  public:
    template <class MatrixT>
    FLANNTree( MatrixBase<MatrixT> const& features ) : m_index_ptr(NULL), m_features_cast( features ) {
      construct_index( (void*)&m_features_cast(0,0), m_features_cast.rows(), m_features_cast.cols() );
    }

    ~FLANNTree();

    // Multiple query access via VW's Matrix
    template <class MatrixT>
    void knn_search( MatrixBase<MatrixT> const& query,
                     Vector<int>& indices,
                     Vector<FloatT>& dists,
                     size_t knn ) {

      // Convert query to our type
      Matrix<FloatT> query_cast = query.impl();
      knn_search_help( (void*)&query_cast(0,0), query_cast.rows(), query_cast.cols(),
                  indices, dists, knn );
    }

    // Single query access via ASP's Matrix
    template <class VectorT>
    void knn_search( VectorBase<VectorT> const& query,
                     Vector<int>& indices,
                     Vector<FloatT>& dists,
                     size_t knn ) {

      // Convert query to our type
      Vector<FloatT> query_cast = query.impl();
      knn_search_help( (void*)&query_cast[0], 1, query_cast.size(),
                  indices, dists, knn );
    }

    size_t size1() const;
    size_t size2() const;
  };

}} // end namespace vw::math

#endif//__VW_MATH_FLANNTREE_H__
