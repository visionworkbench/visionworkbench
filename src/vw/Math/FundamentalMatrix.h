// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file FundamentalMatrix.h
///
/// Assorted useful geometric functors for calculating the fundamental
/// matrix and working with epipoles.
///
#ifndef __MATH_FUNDAMENTAL_MATRIX_H__
#define __MATH_FUNDAMENTAL_MATRIX_H__

#include <vector>

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>

namespace vw {
namespace math {

  // This fitting functor attempts to find a 7 DoF matrix that charaterizes
  // the relationship of epipolar lines between to images. This implements
  // the algorithm defined on pg 281 of Multiview Geometry (aka Bible)
  struct Fundamental7FittingFunctor {
    typedef vw::Matrix<double> result_type;
    Matrix<double> m_nullspace;
    std::vector<double> m_solutions;

    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const {
      return 7;
    }

    // In this algorithm there can be multiple solutions for the F matrix
    // here we allow the user to see them
    uint num_solutions() const { return m_solutions.size(); }

    Matrix<double> fundamental_matrix( uint which = 0 ) const {
      Matrix<double> output(3,3);
      VW_ASSERT( m_nullspace.rows() == 9 && m_nullspace.cols() == 2,
                 vw::ArgumentErr() << "FundamentalMatrixFittingFunctor7::operator() must have been called once." );
      uint current_index = 0;
      double a = m_solutions[which];
      double ia = 1 - a;
      for ( uint i = 0; i < 3; i++ ) {
        for ( uint j = 0; j < 3; j++ ) {
          output(i,j) = a*m_nullspace(current_index,0)+ia*m_nullspace(current_index,1);
          current_index++;
        }
      }
      return output;
    }

    // Interface to solve for F matrix. Will throw error if more than
    // 7 elements
    template <class ContainerT>
    vw::Matrix<double> operator()( std::vector<ContainerT> const& p1,
                                   std::vector<ContainerT> const& p2 ) {
      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "p1 and p2's size not equal" );
      VW_ASSERT( p1.size() == 7,
                 vw::ArgumentErr() << "Only seven elements are used in the 7 point Fundamental Matrix algorithm" );
      VW_ASSERT( p1[0].size() == 3 && p1[0][2] == 1,
                 vw::ArgumentErr() << "p1 does not appert to be normalized homogeneous 2D vectors." );

      // Building A-Matrix
      Matrix<double> A(7,9);
      for ( uint i = 0; i < 7; i++ ) {
        A(i,0) = p1[i].x()*p2[i].x();
        A(i,1) = p2[i].x()*p1[i].y();
        A(i,2) = p2[i].x();
        A(i,3) = p2[i].y()*p1[i].x();
        A(i,4) = p2[i].y()*p1[i].y();
        A(i,5) = p2[i].y();
        A(i,6) = p1[i].x();
        A(i,7) = p1[i].y();
        A(i,8) = 1;
      }

      // Matrix9x2
      Matrix<double> n = nullspace(A);
      VW_ASSERT( n.cols() == 2,
                 vw::MathErr() << "A Matrix has incorrect nullity due insufficient rank or threshold error in math ops. Did you input the data correctly? A = " << A << " Rank = " << rank(A) << " Null = " << nullity(A) << "\n" );
      m_nullspace = n;

      // Solving for alpha cubic
      Vector<double,4> acubic;
      // The following equations are expansions of 0=det(a*F1+(1-a)*F2)
      acubic[0] = n(0,1)*n(4,1)*n(8,1) - n(0,1)*n(5,1)*n(7,1) - n(1,1)*n(3,1)*n(8,1) + n(1,1)*n(5,1)*n(6,1) + n(2,1)*n(3,1)*n(7,1) - n(2,1)*n(4,1)*n(6,1);
      acubic[1] = n(0,0)*n(4,1)*n(8,1) - n(0,0)*n(5,1)*n(7,1) + n(0,1)*n(4,0)*n(8,1) + n(0,1)*n(4,1)*n(8,0) - n(0,1)*n(5,0)*n(7,1) - n(0,1)*n(5,1)*n(7,0) - n(1,0)*n(3,1)*n(8,1) + n(1,0)*n(5,1)*n(6,1) - n(1,1)*n(3,0)*n(8,1) - n(1,1)*n(3,1)*n(8,0) + n(1,1)*n(5,0)*n(6,1) + n(1,1)*n(5,1)*n(6,0) + n(2,0)*n(3,1)*n(7,1) - n(2,0)*n(4,1)*n(6,1) + n(2,1)*n(3,0)*n(7,1) + n(2,1)*n(3,1)*n(7,0) - n(2,1)*n(4,0)*n(6,1) - n(2,1)*n(4,1)*n(6,0) - 3*n(0,1)*n(4,1)*n(8,1) + 3*n(0,1)*n(5,1)*n(7,1) + 3*n(1,1)*n(3,1)*n(8,1) - 3*n(1,1)*n(5,1)*n(6,1) - 3*n(2,1)*n(3,1)*n(7,1) + 3*n(2,1)*n(4,1)*n(6,1);
      acubic[2] = n(0,0)*n(4,0)*n(8,1) + n(0,0)*n(4,1)*n(8,0) - n(0,0)*n(5,0)*n(7,1) - n(0,0)*n(5,1)*n(7,0) + n(0,1)*n(4,0)*n(8,0) - n(0,1)*n(5,0)*n(7,0) - n(1,0)*n(3,0)*n(8,1) - n(1,0)*n(3,1)*n(8,0) + n(1,0)*n(5,0)*n(6,1) + n(1,0)*n(5,1)*n(6,0) - n(1,1)*n(3,0)*n(8,0) + n(1,1)*n(5,0)*n(6,0) + n(2,0)*n(3,0)*n(7,1) + n(2,0)*n(3,1)*n(7,0) - n(2,0)*n(4,0)*n(6,1) - n(2,0)*n(4,1)*n(6,0) + n(2,1)*n(3,0)*n(7,0) - n(2,1)*n(4,0)*n(6,0)  - 2*n(0,0)*n(4,1)*n(8,1) + 2*n(0,0)*n(5,1)*n(7,1) - 2*n(0,1)*n(4,0)*n(8,1) - 2*n(0,1)*n(4,1)*n(8,0) + 2*n(0,1)*n(5,0)*n(7,1) + 2*n(0,1)*n(5,1)*n(7,0) + 2*n(1,0)*n(3,1)*n(8,1) - 2*n(1,0)*n(5,1)*n(6,1) + 2*n(1,1)*n(3,0)*n(8,1) + 2*n(1,1)*n(3,1)*n(8,0) - 2*n(1,1)*n(5,0)*n(6,1) - 2*n(1,1)*n(5,1)*n(6,0) - 2*n(2,0)*n(3,1)*n(7,1) + 2*n(2,0)*n(4,1)*n(6,1) - 2*n(2,1)*n(3,0)*n(7,1) - 2*n(2,1)*n(3,1)*n(7,0) + 2*n(2,1)*n(4,0)*n(6,1) + 2*n(2,1)*n(4,1)*n(6,0) + 3*n(0,1)*n(4,1)*n(8,1) - 3*n(0,1)*n(5,1)*n(7,1) - 3*n(1,1)*n(3,1)*n(8,1) + 3*n(1,1)*n(5,1)*n(6,1) + 3*n(2,1)*n(3,1)*n(7,1) - 3*n(2,1)*n(4,1)*n(6,1);
      acubic[3] = n(0,0)*n(4,0)*n(8,0) - n(0,0)*n(5,0)*n(7,0) - n(1,0)*n(3,0)*n(8,0) + n(1,0)*n(5,0)*n(6,0) + n(2,0)*n(3,0)*n(7,0) - n(2,0)*n(4,0)*n(6,0) - n(0,0)*n(4,0)*n(8,1) - n(0,0)*n(4,1)*n(8,0) + n(0,0)*n(5,0)*n(7,1) + n(0,0)*n(5,1)*n(7,0) - n(0,1)*n(4,0)*n(8,0) + n(0,1)*n(5,0)*n(7,0) + n(1,0)*n(3,0)*n(8,1) + n(1,0)*n(3,1)*n(8,0) - n(1,0)*n(5,0)*n(6,1) - n(1,0)*n(5,1)*n(6,0) + n(1,1)*n(3,0)*n(8,0) - n(1,1)*n(5,0)*n(6,0) - n(2,0)*n(3,0)*n(7,1) - n(2,0)*n(3,1)*n(7,0) + n(2,0)*n(4,0)*n(6,1) + n(2,0)*n(4,1)*n(6,0) - n(2,1)*n(3,0)*n(7,0) + n(2,1)*n(4,0)*n(6,0) + n(0,0)*n(4,1)*n(8,1) - n(0,0)*n(5,1)*n(7,1) + n(0,1)*n(4,0)*n(8,1) + n(0,1)*n(4,1)*n(8,0) - n(0,1)*n(5,0)*n(7,1) - n(0,1)*n(5,1)*n(7,0) - n(1,0)*n(3,1)*n(8,1) + n(1,0)*n(5,1)*n(6,1) - n(1,1)*n(3,0)*n(8,1) - n(1,1)*n(3,1)*n(8,0) + n(1,1)*n(5,0)*n(6,1) + n(1,1)*n(5,1)*n(6,0) + n(2,0)*n(3,1)*n(7,1) - n(2,0)*n(4,1)*n(6,1) + n(2,1)*n(3,0)*n(7,1) + n(2,1)*n(3,1)*n(7,0) - n(2,1)*n(4,0)*n(6,1) - n(2,1)*n(4,1)*n(6,0) - n(0,1)*n(4,1)*n(8,1) + n(0,1)*n(5,1)*n(7,1) + n(1,1)*n(3,1)*n(8,1) - n(1,1)*n(5,1)*n(6,1) - n(2,1)*n(3,1)*n(7,1) + n(2,1)*n(4,1)*n(6,1);
      acubic /= acubic[3];

      //Finding for solutions of cubic function
      Matrix3x3 companion;
      companion(0,2) = -acubic[0];
      companion(1,0) = 1;
      companion(1,2) = -acubic[1];
      companion(2,1) = 1;
      companion(2,2) = -acubic[2];

      Vector<std::complex<double> > roots;
      eigen(companion,roots);

      m_solutions.clear();
      for ( uint i = 0; i < roots.size(); i++ ) {
        if ( roots[i].imag() < 1e-10 && roots[i].imag() > -1e-10 ) {
          m_solutions.push_back(roots[i].real());
        }
      }
      VW_ASSERT( roots.size() > 0,
                 vw::MathErr() << "FundamentalMatrixFittingFunctor7 didn't find a solution.\n" );

      return this->fundamental_matrix();
    }
  };

  // The normalized 8 point and greater algorithm for F matrix. On
  // page 282. This algorithm is not gold standard as it doesn't
  // enforce output F has norm =1 to 1. It only enforces that the det is
  // zero.
  struct Fundamental8FittingFunctor {
    typedef vw::Matrix<double> result_type;

    /// This requires 8 and greater
    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const { return 8; }

    /// Solve for Normalization Similarity Matrix used for noise rej.
    template <class ContainerT>
    vw::Matrix<double> NormSimilarity( std::vector<ContainerT> const& pts ) const {
      unsigned num_points = pts.size();
      unsigned dimension = pts[0].size();

      Matrix<double> translation;
      translation.set_identity(dimension);

      Vector<double> sum;
      sum.set_size(dimension-1);
      for ( unsigned i = 0; i < num_points; i++ )
        sum+=subvector(pts[i],0,dimension-1);
      sum /= num_points;
      for ( unsigned i = 0; i < dimension-1; i++ )
        translation(i,dimension-1) = -sum(i);

      std::vector<Vector<double> > pts_int;
      for ( uint i = 0; i < pts.size(); i++ )
        pts_int.push_back( translation*pts[i] );

      Matrix<double> scalar;
      scalar.set_identity(dimension);
      double scale = 0;
      for ( unsigned i = 0; i < num_points; i++ )
        scale += norm_2( subvector(pts_int[i],0,dimension-1) );
      scale = num_points*sqrt(2.)/scale;
      scalar *= scale;
      scalar(dimension-1,dimension-1) = 1;
      return scalar*translation;
    }

    // Interface to solve for F matrix
    template <class ContainerT>
    vw::Matrix<double> operator()( std::vector<ContainerT> const& p1,
                                   std::vector<ContainerT> const& p2 ) const {
      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "Cannot compute Fundamental M. p1 and p2 not the same size." );
      VW_ASSERT( p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute Fundamental M. Insufficient data.");
      VW_ASSERT( p1[0].size() == 3 && p1[0][2] == 1,
                 vw::ArgumentErr() << "Cannot compute Fundamental M. Vectors don't seem to be 2D homogeneous vectors." );

      // Normalizing
      Matrix<double> S_in = NormSimilarity(p1);
      Matrix<double> S_out = NormSimilarity(p2);
      std::vector<Vector<double> > input_p, output_p;
      for ( uint i = 0; i < p1.size(); i++ ) {
        input_p.push_back( S_in*p1[i] );
        output_p.push_back( S_out*p2[i] );
      }

      // Linear solution
      Matrix<double> A(p1.size(),9);
      for ( uint i = 0; i < p1.size(); i++ ) {
        A(i,0) = input_p[i][0]*output_p[i][0];
        A(i,1) = output_p[i][0]*input_p[i][1];
        A(i,2) = output_p[i][0];
        A(i,3) = output_p[i][1]*input_p[i][0];
        A(i,4) = output_p[i][1]*input_p[i][1];
        A(i,5) = output_p[i][1];
        A(i,6) = input_p[i][0];
        A(i,7) = input_p[i][1];
        A(i,8) = 1;
      }
      Matrix<double> n = nullspace(A);

      // Rearranging
      Matrix<double,3,3> rn;
      int i = 0;
      for ( Matrix<double,3,3>::iterator it = rn.begin(); it != rn.end(); it++ ) {
        (*it) = n(i,0);
        i++;
      }

      // Constraint enforcement
      Matrix<double> U, V;
      Vector<double> S;
      complete_svd(rn, U, S, V);
      S[2] = 0;
      rn = U*diagonal_matrix(S)*V;

      // Denormalization
      rn = inverse(S_out)*rn*S_in;

      return rn;
    }

  };


}}

#endif//__MATH_FUNDAMENTAL_MATRIX_H__
