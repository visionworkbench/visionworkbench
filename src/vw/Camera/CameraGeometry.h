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

#include <vw/Core/Exception.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/LevenbergMarquardt.h>

namespace vw {

// Forward Declare
namespace ip { struct InterestPoint; }

namespace camera {

  // Transform Functors
  template <class VectorT>
  Vector3 Convert2HomogenousVec2( VectorT const& v ) {
    if ( v.size() == 2 )
      return Vector3( v[0], v[1], 1 );
    return Vector3( v[0], v[1], v[2] );
  }
  template <>
  Vector3 Convert2HomogenousVec2( ip::InterestPoint const& ip );

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
    Matrix<double> BasicDLT( std::vector<Vector<double> > const& input,
                             std::vector<Vector<double> > const& output ) const;

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
      CameraMatrixModelLMA( std::vector<Vector<double> > const& input,
                            std::vector<Vector<double> > const& output );

      // Evaluator
      result_type operator()( domain_type const& x ) const;

      // Helper functions
      Vector<double> flatten( Matrix<double> const& input ) const;
      Matrix<double> unflatten( Vector<double> const& input ) const;
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

  struct FundamentalMatrixSampsonErrorMetric {
    template <class RelationT, class ContainerT>
    double operator()( RelationT const& F,
                       ContainerT const& p1,
                       ContainerT const& p2 ) const {
      return fabs(transpose(p2)*F*p1) + fabs(transpose(p1)*transpose(F)*p2);
    }
  };

  struct FundamentalMatrixDistanceErrorMetric {
    template <class RelationT, class ContainerT>
    double operator()( RelationT const& F,
                       ContainerT const& p1,
                       ContainerT const& p2 ) const {
      Vector3 line = F*p1;
      return fabs(dot_prod(line,p2))/norm_2(subvector(line,0,2));
    }
  };

  /// Fundamental Matrix solver using normalized 8 point algorithm.
  /// - Page 282 or Algorithm 11.1 in Multiple View Geometry
  struct FundamentalMatrix8PFittingFunctor : public SimilarityNormalizingFunctor {
    typedef vw::Matrix<double> result_type;

    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& /*example*/) const { return 8; }

    /// Solve for the fundamental matrix F using two sets of observation of a set
    ///  of points from the same camera at different positions.
    template <class ContainerT>
    result_type
    operator()( std::vector<ContainerT> const& p1,
                std::vector<ContainerT> const& p2,
                vw::Matrix<double> const& /*seed_input*/ = vw::Matrix<double>() ) const {

      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "Fundamental Matrix requires equal number of input and output measures." );

      // Converting to internal format
      std::vector<Vector<double> > input(p1.size()), output(p2.size());
      std::transform( p1.begin(), p1.end(), input.begin(),  Convert2HomogenousVec2<ContainerT> );
      std::transform( p2.begin(), p2.end(), output.begin(), Convert2HomogenousVec2<ContainerT> );

      // Normalizing
      Matrix<double> S_in  = NormSimilarity( input  );
      Matrix<double> S_out = NormSimilarity( output );
      std::for_each( input.begin(),  input.end(),  MatrixApplyFunc<Vector<double> >( S_in  ) );
      std::for_each( output.begin(), output.end(), MatrixApplyFunc<Vector<double> >( S_out ) );

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

      VW_ASSERT( math::rank(A) >= 8, MathErr() << "Measurements produce rank deficient A." );

      // Pulling singular vector of smallest singular value of A
      Matrix<double> U,VT;
      Vector<double> S;
      svd(A,U,S,VT);
      Matrix3x3 F;
      int i = 0;
      for ( Matrix<double,3,3>::iterator it = F.begin(); it != F.end(); it++ ) {
        (*it) = VT(VT.rows()-1,i);
        i++;
      }

      // Constraint Enforcement
      svd(F,U,S,VT);
      S[2] = 0;
      F = U*diagonal_matrix(S)*VT;

      // Denormalizing
      return transpose(S_out)*F*S_in;
    }
  };

  // Fundamental Matrix solver using the Maximum Likelihood method
  // Page 285 or Algorithm 11.3 in Multiple View Geometry
  struct FundamentalMatrixMLFittingFunctor {
    typedef vw::Matrix<double> result_type;

    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& /*example*/) const { return 8; }

    // Skew symmetric matrix or []x operation
    template <class VectorT>
    Matrix<double> skew_symm( VectorBase<VectorT> const& b ) const {
      VW_ASSERT( b.impl().size() == 3,
                 vw::ArgumentErr() << "Skew Symmetric is only elemented for 3 element vectors." );
      VectorT const& v = b.impl();
      Matrix3x3 skew;
      skew(0,1) = -v[2];
      skew(0,2) = v[1];
      skew(1,0) = v[2];
      skew(1,2) = -v[0];
      skew(2,0) = -v[1];
      skew(2,1) = v[0];
      return skew;
    }

    template <class VectorT>
    Vector4
    OneSideTriangulation( Matrix<double> const& P2,
                          VectorBase<VectorT> const& meas1,
                          VectorBase<VectorT> const& meas2) const {
      VectorT const& v1 = meas1.impl();
      VectorT const& v2 = meas2.impl();
      Matrix<double> A(4,4);
      select_row(A,0) = Vector4(-1,0,v1[0],0);
      select_row(A,1) = Vector4(0,-1,v1[1],0);
      select_row(A,2) = v2[0]*select_row(P2,2)-select_row(P2,0);
      select_row(A,3) = v2[1]*select_row(P2,2)-select_row(P2,1);

      Matrix<double> U, VT;
      Vector<double> S;
      svd(A,U,S,VT);

      Vector<double> solution = select_row(VT,3);
      solution /= solution[3];
      return solution;
    }

    // Non-linear solution to an over-constrained problem
    class ProjectiveModelLMA : public math::LeastSquaresModelBase<ProjectiveModelLMA> {

    public:
      typedef Vector<double> result_type; // reprojective error
      typedef Vector<double> domain_type; // searchspace 3n+12 variables
      typedef Matrix<double> jacobian_type;

      // Evaluator
      result_type operator()( domain_type const& x )  const;
    };

    /// Solve for the fundamental matrix F using two sets of observation of a set
    ///  of points from the same camera at different positions.
    template <class ContainerT>
    result_type
    operator()( std::vector<ContainerT> const& p1,
                std::vector<ContainerT> const& p2,
                vw::Matrix<double> const& seed_input = vw::Matrix<double>() ) const {

      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "Fundamental Matrix requires equal number of input and output measures." );

      // Determine if to use 8P
      if ( p1.size() == 8 )
        return FundamentalMatrix8PFittingFunctor()(p1,p2);

      // Converting to internal format
      std::vector<Vector<double> > input(p1.size()), output(p2.size());
      std::transform( p1.begin(), p1.end(), input.begin(),
                      Convert2HomogenousVec2<ContainerT> );
      std::transform( p2.begin(), p2.end(), output.begin(),
                      Convert2HomogenousVec2<ContainerT> );

      // Getting epipole
      Matrix<double> U, VT;
      Vector<double> S;
      svd(seed_input,U,S,VT);
      Vector<double> epipole_prime = select_col(U,2);

      // Seed initial camera guesses
      Matrix<double> P2(3,4);
      submatrix(P2,0,0,3,3) = skew_symm(epipole_prime)*seed_input;
      select_col(P2,3) = epipole_prime;

      // Solve for projective
      ProjectiveModelLMA model;
      Vector<double> seed(12+3*input.size());
      subvector(seed,0,4) = select_row(P2,0);
      subvector(seed,4,4) = select_row(P2,1);
      subvector(seed,8,4) = select_row(P2,2);

      // Triangulating
      for ( unsigned i = 0, j = 12; i < input.size(); i++, j+=3 ) {
        Vector4 projection = OneSideTriangulation(P2,input[i],output[i]);
        projection /= projection[3];
        subvector(seed,j,3) = subvector(projection,0,3);
      }

      // Setting Objective
      Vector<double> objective(4*input.size());
      for ( unsigned i = 0; i < input.size(); i++ ) {
        subvector(objective,4*i,2) = subvector(input[i],0,2);
        subvector(objective,4*i+2,2) = subvector(output[i],0,2);
      }

      int status = 0;
      Vector<double> result =
        levenberg_marquardt( model, seed,
                             objective, status );

      select_row(P2,0) = subvector(result,0,4);
      select_row(P2,1) = subvector(result,4,4);
      select_row(P2,2) = subvector(result,8,4);

      Matrix<double> F =
        skew_symm(select_col(P2,3))*submatrix(P2,0,0,3,3);

      return F;
    }
  };

}} // end vw::camera

#endif//__CAMERA_CAMERA_GEOMETRY_H__
