// __BEGIN_LICENSE__
// __END_LICENSE__

/// \file CameraGeometry.h
///
/// Contains assorted functions for solving for camera matrices.
/// This is an expansion on Math/Geometry.h but requires that
/// container p1 is Vector3 and container p2 is Vector4 since
/// cameras map 3D points to 2D locations.
///
#ifndef __CAMERA_CAMERA_GEOMETRY_H__
#define __CAMERA_CAMERA_GEOMETRY_H__

#include <boost/foreach.hpp>

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/LevenbergMarquardt.h>

namespace vw {
namespace camera {

  // Transform Functors
  template <class VectorT>
  Vector3 Convert2HomogenousVec2( VectorT const& v ) {
    if ( v.size() == 2 )
      return Vector3( v[0], v[1], 1 );
    return Vector3( v[0], v[1], v[2] );
  }

  // Measurement Prep Functors
  struct SimilarityNormalizingFunctor {
    // Matrix apply functor
    template <class ContainerT>
    struct MatrixApplyFunc {
      Matrix<double> m_matrix;
      template <class MatrixT>
      MatrixApplyFunc( MatrixBase<MatrixT> const& m ) : m_matrix(m.impl()) {}

      // In-place modifier (use with for_each)
      void operator()( ContainerT & v ) const { v = m_matrix*v; }

      // Copy modifier (use with transform)
      ContainerT operator()( ContainerT const& v ) const { return m_matrix*v; }
    };

    // Solve for normalizing similarity matrix
    template <class ContainerT>
    vw::Matrix<double>
    NormSimilarity( std::vector<ContainerT> const& pts ) const {
      unsigned num_points = pts.size();
      unsigned dimension = pts[0].size();

      Vector<double> translation;
      translation.set_size(dimension-1);
      BOOST_FOREACH( const ContainerT& p, pts ) {
        translation += subvector(p,0,dimension-1);
      }
      translation /= num_points;

      double scale = 0;
      BOOST_FOREACH( const ContainerT& p, pts ) {
        Vector<double> delta = subvector(p,0,dimension-1) - translation;
        scale += norm_2(delta);
      }
      scale = num_points*sqrt(2.)/scale;

      Matrix<double> s;
      s.set_size(dimension,dimension);
      s.set_identity();
      for ( unsigned i = 0; i < dimension-1; i++ ) {
        s(i,i) *= scale;
        s(i,dimension-1) = -scale*translation[i];
      }
      return s;
    }
  };

  // Camera Matrix Functor
  //
  // Solves for 3x4 matrix which describes the mapping of a pinhole
  // camera from 3D points in the world to 2D points in an image.

  struct CameraMatrixErrorMetric {
    template <class RelationT, class ContainerT>
    double operator()( RelationT const& H,
                       ContainerT const& p1,
                       ContainerT const& p2 ) const {
      Vector<double> projection = H*p1;
      projection /= projection[2];
      return norm_2(p2-projection);
    }
  };

  // This fitting functor fits a 3x4 camera matrix P
  struct CameraMatrixFittingFunctor : public SimilarityNormalizingFunctor {
    typedef vw::Matrix<double> result_type;

    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& /*example*/) const { return 6; }

    // Simple linear solution
    vw::Matrix<double> BasicDLT( std::vector<Vector<double> > const& input,
                                 std::vector<Vector<double> > const& output ) const {
      VW_ASSERT( input[0].size() == 4 && output[0].size() == 3,
                 vw::ArgumentErr() << "Camera Matrix requires Vector4 inputs and Vector3 outputs." );

      vw::Matrix<double,12,12> A;
      for ( unsigned i = 0; i < 6; i++ ) { // Measure iterator
        for ( unsigned j = 0; j < 4; j++ ) { // Measure's X elem iterator
          // Filling w*Xt
          A(2*i+1,j) = output[i][2]*input[i][j];
          // Filling -w*Xt
          A(2*i,j+4) = -output[i][2]*input[i][j];
          // Filling y*Xt
          A(2*i,j+8) = output[i][1]*input[i][j];
          // Filling -x*Xt
          A(2*i+1,j+8) = -output[i][0]*input[i][j];
        }
      }

      Matrix<double,3,4> p;
      // SVD for smallest singular value
      Matrix<double> U,VT;
      Vector<double> S;
      svd(A,U,S,VT);
      submatrix(p,0,0,1,4) = submatrix(VT,VT.rows()-1,0,1,4);
      submatrix(p,1,0,1,4) = submatrix(VT,VT.rows()-1,4,1,4);
      submatrix(p,2,0,1,4) = submatrix(VT,VT.rows()-1,8,1,4);
      return p;
    }

    // Non-linear solution to an over-constrained problem
    class CameraMatrixModelLMA : public math::LeastSquaresModelBase<CameraMatrixModelLMA> {
      std::vector<Vector<double> > m_world_input;
      std::vector<Vector<double> > m_image_output;

    public:
      // What is returned by evaluating functor. It's the reprojection error.
      typedef Vector<double> result_type;
      // Define the search space. This is the camera matrix flattened out.
      typedef Vector<double> domain_type;
      // Jacobian form
      typedef Matrix<double> jacobian_type;

      // Constructor
      inline CameraMatrixModelLMA( std::vector<Vector<double> > input,
                                   std::vector<Vector<double> > output ) :
        m_world_input(input), m_image_output(output) {}

      // Evaluator
      inline result_type operator()( domain_type const& x ) const {
        result_type output;
        output.set_size( m_world_input.size() );
        Matrix<double,3,4> P = unflatten(x);

        for ( uint32 i = 0; i < m_world_input.size(); i++ ) {
          Vector3 reproj = P*m_world_input[i];
          reproj /= reproj[2];
          output[i] = norm_2( subvector(m_image_output[i],0,2) - subvector(reproj,0,2) );
        }
        return output;
      }

      // Help functions
      Vector<double> flatten( Matrix<double> const& input ) const {
        Vector<double,12> output;
        for ( uint8 i = 0; i < 3; i++ ) {
          for ( uint8 j = 0; j < 4; j++ ) {
            output(4*i+j) = input(i,j);
          }
        }
        return output;
      }

      Matrix<double> unflatten( Vector<double> const& input ) const {
        Matrix<double,3,4> output;
        for ( uint8 i = 0; i < 3; i++ ) {
          for ( uint8 j = 0; j < 4; j++ ) {
            output(i,j) = input(4*i+j);
          }
        }
        return output;
      }
    };

    // Interface for solving for Camera Matrix (switches between DLT
    // and Gold Standard)
    template <class ContainerT>
    result_type
    operator()( std::vector<ContainerT> const& p1,
                std::vector<ContainerT> const& p2,
                vw::Matrix<double> const& seed_input = vw::Matrix<double>() ) const {
      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "Camera Matrix requires equal number of input and output measures." );
      VW_ASSERT( p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute camera matrix. Insufficient measurements." );
      VW_ASSERT( p1[0].size() == 4 && p2[0].size() == 3,
                 vw::ArgumentErr() << "Incorrect input format: p1 needs to Vector4 and p2 need to be Vector3." );

      // Converting to internal format
      std::vector<Vector<double> > input, output;
      BOOST_FOREACH( const ContainerT& p, p1 ) {
        input.push_back( Vector4( p[0], p[1], p[2], p[3] ) );
      }
      BOOST_FOREACH( const ContainerT& p, p2 ) {
        output.push_back( Vector3( p[0], p[1], p[2] ) );
      }

      // Normalizing
      Matrix<double> S_in = NormSimilarity( input );
      Matrix<double> S_out = NormSimilarity( output );
      std::for_each( input.begin(), input.end(),
                     MatrixApplyFunc<Vector<double> >( S_in ) );
      std::for_each( output.begin(), output.end(),
                     MatrixApplyFunc<Vector<double> >( S_out ) );

      Matrix<double> p;
      if ( p1.size() == min_elements_needed_for_fit(p1[0] ) ) {
        // Use DLT
        Matrix<double> P_prime = BasicDLT( input, output );

        p = inverse(S_out)*P_prime*S_in;
      } else {
        Matrix<double> seed_matrix = seed_input;

        // Perform DLT if needed
        if ( seed_matrix.cols() == 0 )
          seed_matrix = BasicDLT( input, output );

        // Iterative solution be here
        CameraMatrixModelLMA model( input, output );
        Vector<double> seed = model.flatten( seed_matrix );
        int status = 0;
        Vector<double> objective;
        objective.set_size( input.size() );
        Vector<double> result = levenberg_marquardt( model, seed,
                                                     objective, status );
        seed_matrix = model.unflatten( result );

        // Denormalization
        p = inverse(S_out)*seed_matrix*S_in;
      }

      p /= p(2,3);
      return p;
    }
  };

  // Fundamental Matrix Functor
  //
  // Solves for 3x3 matrix which relates stereo images. Fundamental
  // matrix ties points in first image to lines in the second image. The
  // lines are defined be the epipole (which is the location of the first
  // image in the perspective of the second camera (and vice versa)).

  struct FundamentalMatrixErrorMetric {
    template <class RelationT, class ContainerT>
    double operator()( RelationT const& F,
                       ContainerT const& p1,
                       ContainerT const& p2 ) const {
      return fabs(transpose(p2)*F*p1) + fabs(transpose(p1)*transpose(F)*p2);
    }
  };

  // Fundamental Matrix solver using normalized 8 point algorithm. (pg
  // 282 or Algorithm 11.1 in Multiple View Geometry )
  struct FundamentalMatrix8PFittingFunctor : public SimilarityNormalizingFunctor {
    typedef vw::Matrix<double> result_type;

    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& /*example*/) const { return 8; }

    template <class ContainerT>
    result_type
    operator()( std::vector<ContainerT> const& p1,
                std::vector<ContainerT> const& p2,
                vw::Matrix<double> const& /*seed_input*/ = vw::Matrix<double>() ) const {

      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "Fundamental Matrix requires equal number of input and output measures." );

      // Converting to internal format
      std::vector<Vector<double> > input(p1.size()), output(p2.size());
      std::transform( p1.begin(), p1.end(), input.begin(),
                      Convert2HomogenousVec2<ContainerT> );
      std::transform( p2.begin(), p2.end(), output.begin(),
                      Convert2HomogenousVec2<ContainerT> );

      // Normalizing
      Matrix<double> S_in = NormSimilarity( input );
      Matrix<double> S_out = NormSimilarity( output );
      std::for_each( input.begin(), input.end(),
                     MatrixApplyFunc<Vector<double> >( S_in ) );
      std::for_each( output.begin(), output.end(),
                     MatrixApplyFunc<Vector<double> >( S_out ) );

      // Constructing A
      Matrix<double> A(p1.size(),9);
      for ( size_t i = 0; i < p1.size(); i++ ) {
        A(i,0) = output[i][0]*input[i][0];
        A(i,1) = output[i][0]*input[i][1];
        A(i,2) = output[i][0];
        A(i,3) = output[i][1]*input[i][0];
        A(i,4) = output[i][1]*input[i][1];
        A(i,5) = output[i][1];
        A(i,6) = input[i][0];
        A(i,7) = input[i][1];
        A(i,8) = 1;
      }

      VW_ASSERT( rank(A) >= 8,
                 vw::MathErr() << "Measurements produce rank deficient A." );

      // Pulling singular vector of smallest singular value of A
      Matrix<double> n = nullspace(A);
      Matrix3x3 F;
      int i = 0;
      for ( Matrix<double,3,3>::iterator it = F.begin();
            it != F.end(); it++ ) {
        (*it) = n(i,0);
        i++;
      }

      // Constraint Enforcement
      Matrix<double> U,VT;
      Vector<double> S;
      svd(F,U,S,VT);
      S[2] = 0;
      F = U*diagonal_matrix(S)*VT;

      // Denormalizing
      return transpose(S_out)*F*S_in;
    }
  };

}} // end vw::camera

#endif//__CAMERA_CAMERA_GEOMETRY_H__
