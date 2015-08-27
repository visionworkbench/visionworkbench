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


#include <vw/Math/FLANNTree.h>
#include <flann/flann.hpp>

namespace vw {
namespace math {

  template <>
  size_t FLANNTree<float>::knn_search_help( void* data_ptr, size_t rows, size_t cols,
                                            Vector<int  >& indices,
                                            Vector<float>& dists,
                                            size_t knn ) {
    // Constrain the number of results that we can return to the number of loaded objects
    size_t maxNumReturns = m_features_cast.rows();
    if (knn > maxNumReturns)
      knn = maxNumReturns;

    if ( indices.size() != knn )
      indices.set_size( knn );
    if ( dists.size() != knn )
      dists.set_size( knn );

    flann::Matrix<float> query_mat ( (float*)data_ptr, rows, cols );
    flann::Matrix<int  > indice_mat( &indices[0], 1, knn );
    flann::Matrix<float> dists_mat ( &dists[0], 1, dists.size() );
    ((flann::Index<flann::L2<float> >*)(m_index_ptr))->knnSearch( query_mat, indice_mat, dists_mat, knn, flann::SearchParams(128) );
    return knn;
  }

  template <>
  size_t FLANNTree<double>::knn_search_help( void* data_ptr, size_t rows, size_t cols,
                                             Vector<int   >& indices,
                                             Vector<double>& dists,
                                             size_t knn ) {
    // Constrain the number of results that we can return to the number of loaded objects
    size_t maxNumReturns = m_features_cast.rows();
    if (knn > maxNumReturns)
      knn = maxNumReturns;

    if ( indices.size() != knn )
      indices.set_size( knn );
    if ( dists.size() != knn )
      dists.set_size( knn );

    flann::Matrix<double> query_mat ( (double*)data_ptr, rows, cols );
    flann::Matrix<int   > indices_mat( &indices[0], 1, knn );
    flann::Matrix<double> dists_mat ( &dists[0], 1, dists.size() );
    ((flann::Index<flann::L2<double> >*)(m_index_ptr))->knnSearch( query_mat, indices_mat, dists_mat, knn, flann::SearchParams(128) );
    return knn;
  }

  template <>
  void FLANNTree<float>::construct_index( void* data_ptr, size_t rows, size_t cols ) {
    if ( m_index_ptr != NULL )
      vw_throw( IOErr() << "FLANNTree: Void ptr is not null, this is unexpected." );
    m_index_ptr = new flann::Index<flann::L2<float> >( flann::Matrix<float>( (float*)data_ptr, rows, cols ),
                                                       flann::KDTreeIndexParams(4),
                                                       flann::L2<float>() );
    ((flann::Index<flann::L2<float> >*)m_index_ptr)->buildIndex();
  }
  template <>
  void FLANNTree<double>::construct_index( void* data_ptr, size_t rows, size_t cols ) {
    if ( m_index_ptr != NULL )
      vw_throw( IOErr() << "FLANNTree: Void ptr is not null, this is unexpected." );
    m_index_ptr = new flann::Index<flann::L2<double> >( flann::Matrix<double>( (double*)data_ptr, rows, cols ),
                                                       flann::KDTreeIndexParams(4),
                                                       flann::L2<double>() );
    ((flann::Index<flann::L2<double> >*)m_index_ptr)->buildIndex();
  }


  template <>
  FLANNTree<float>::~FLANNTree() {
    delete (flann::Index<flann::L2<float> >*)(m_index_ptr);
  }
  template <>
  FLANNTree<double>::~FLANNTree() {
    delete (flann::Index<flann::L2<double> >*)(m_index_ptr);
  }

  template <>
  size_t FLANNTree<float >::size1() const { return ((flann::Index<flann::L2<float > >*)m_index_ptr)->size(); }
  template <>
  size_t FLANNTree<double>::size1() const { return ((flann::Index<flann::L2<double> >*)m_index_ptr)->size(); }
  template <>
  size_t FLANNTree<float >::size2() const { return ((flann::Index<flann::L2<float > >*)m_index_ptr)->veclen(); }
  template <>
  size_t FLANNTree<double>::size2() const { return ((flann::Index<flann::L2<double> >*)m_index_ptr)->veclen(); }

}}
