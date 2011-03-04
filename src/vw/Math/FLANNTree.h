#ifndef __VW_MATH_FLANNTREE_H__
#define __VW_MATH_FLANNTREE_H__

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <flann/flann.hpp>

namespace vw {
namespace math {

  template <typename DistanceT>
  class FLANNTree {
    flann::Index<DistanceT> m_index;

  public:
    typedef typename DistanceT::ElementType element_type;
    typedef typename DistanceT::ResultType distance_type;

    template <class MatrixT>
    FLANNTree( MatrixBase<MatrixT> const& features,
               flann::IndexParams const& params = flann::KDTreeIndexParams(4),
               DistanceT d = DistanceT() ) :
      m_index( flann::Matrix<typename MatrixT::value_type>( const_cast<typename MatrixT::value_type*>(features.impl().data()), features.impl().rows(), features.impl().cols() ), params, d ) {
      m_index.buildIndex();
    }

    // Multiple query access
    template <class MatrixT>
    void knn_search( MatrixBase<MatrixT> const& query,
                     Vector<int>& indices,
                     Vector<distance_type>& dists,
                     size_t knn,
                     flann::SearchParams const& params = flann::SearchParams(128) ) {
      typedef typename MatrixT::value_type element_type;

      if ( indices.size() != knn )
        indices.set_size( knn );
      if ( dists.size() != knn )
        dists.set_size( knn );

      flann::Matrix<element_type> query_mat( const_cast<element_type*>(query.impl().data()),
                                             query.impl().rows(), query.impl().cols() );
      flann::Matrix<int> indice_mat( &indices[0], 1, indices.size() );
      flann::Matrix<distance_type> dists_mat( &dists[0], 1, dists.size() );
      m_index.knnSearch( query_mat, indice_mat, dists_mat, knn, params );
    }

    // Single query access
    template <class VectorT>
    void knn_search( VectorBase<VectorT> const& query,
                     Vector<int>& indices,
                     Vector<distance_type>& dists,
                     size_t knn,
                     flann::SearchParams const& params = flann::SearchParams(128) ) {
      typedef typename VectorT::value_type element_type;

      if ( indices.size() != knn )
        indices.set_size( knn );
      if ( dists.size() != knn )
        dists.set_size( knn );

      flann::Matrix<element_type> query_mat( const_cast<element_type*>(&query.impl()[0]),
                                             1, query.impl().size() );
      flann::Matrix<int> indice_mat( &indices[0], 1, indices.size() );
      flann::Matrix<distance_type> dists_mat( &dists[0], 1, dists.size() );
      m_index.knnSearch( query_mat, indice_mat, dists_mat, knn, params );
    }

    size_t size1() const { return m_index.size(); }
    size_t size2() const { return m_index.veclen(); }
  };

}} // end namespace vw::math

#endif//__VW_MATH_FLANNTREE_H__
