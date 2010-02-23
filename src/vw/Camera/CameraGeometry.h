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

namespace vw {
namespace camera {

  // This fitting functor fits a 3x4 camera matrix P
  struct CameraMatrixFittingFunctor {

    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const { return 6; }

    // Apply a transform to a list of vectors
    // Applies a transform matrix to a list of points;
    std::vector<Vector<double> > apply_matrix( Matrix<double> const& m,
                                               std::vector<Vector<double> > const& pts ) const {
      std::vector<Vector<double> > out;
      for ( unsigned i = 0; i < pts.size(); i++ ) {
        out.push_back( m*pts[i] );
      }
      return out;
    }

    // Solve for normalizing similarity matrix
    template <class ContainerT>
    vw::Matrix<double> NormSimilarity( std::vector<ContainerT> const& pts ) const {
      unsigned num_points = pts.size();
      unsigned dimension = pts[0].size();

      Vector<double> translation;
      translation.set_size(dimension-1);
      BOOST_FOREACH( ContainerT p, pts ) {
        translation += subvector(p,0,dimension-1);
      }
      translation /= num_points;

      double scale = 0;
      BOOST_FOREACH( ContainerT p, pts ) {
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

    // Simple linear solution
    vw::Matrix<double> BasicDLT( std::vector<Vector<double> > const& input,
                                 std::vector<Vector<double> > const& output ) const {
      VW_ASSERT( input.size() == 6 && output.size() == 6,
                 vw::ArgumentErr() << "Camera Matrix DLT requires exactly 6 measures." );
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

      Matrix<double> nullsp = nullspace(A);
      nullsp /= nullsp(11,0);
      Matrix<double,3,4> p;
      for ( unsigned i = 0; i < 3; i++ )
        for ( unsigned j = 0; j < 4; j++ )
          p(i,j) = nullsp(i*4+j,0);
      return p;
    }

    // Interface for solving for Camera Matrix (switches between DLT and Gold Standard)
    template <class ContainerT>
    vw::Matrix<double> operator()( std::vector<ContainerT> const& p1,
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
      BOOST_FOREACH( ContainerT p, p1 ) {
        input.push_back( Vector4( p[0], p[1], p[2], p[3] ) );
      }
      BOOST_FOREACH( ContainerT p, p2 ) {
        output.push_back( Vector3( p[0], p[1], p[2] ) );
      }

      if ( p1.size() == min_elements_needed_for_fit(p1[0] ) ) {
        // Use DLT
        Matrix<double> S_in = NormSimilarity( input );
        Matrix<double> S_out = NormSimilarity( output );
        std::vector<Vector<double> > input_prime = apply_matrix( S_in, input );
        std::vector<Vector<double> > output_prime = apply_matrix( S_out, output );
        Matrix<double> P_prime = BasicDLT( input_prime, output_prime );
        Matrix<double> P = inverse(S_out)*P_prime*S_in;
        P /= P(2,3);
        return P;
      } else {
        vw_throw( NoImplErr() << "Camera Matrix iterative solution is unfinished." );
      }
    }
  };

}} // end vw::camera

#endif//__CAMERA_CAMERA_GEOMETRY_H__
