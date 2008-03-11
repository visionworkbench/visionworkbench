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

/// \file Geometry.h
/// 
/// Assorted useful geometric routines and functors.
///  
#ifndef __MATH_GEOMETRY_H__
#define __MATH_GEOMETRY_H__

#include <vector>
#include <boost/shared_ptr.hpp>

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/Statistics.h>
#include <vw/Math/LevenbergMarquardt.h>

namespace vw { 
namespace math {

  /// This fitting functor attempts to find a homography (8 degrees of
  /// freedom) that transforms point p1 to match points p2.  This fit
  /// is optimal in a least squares sense.
  struct HomographyFittingFunctor {
    typedef vw::Matrix<double> result_type;

    /// A homography requires at least 4 point matches to determine 8 unknowns
    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const {
      return 4;
    }
    
    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1, 
                                   std::vector<ContainerT> const& p2) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(), 
                 vw::ArgumentErr() << "Cannot compute homography.  p1 and p2 are not the same size." );
      VW_ASSERT( p1.size() != 0 && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute homography.  Insufficient data.\n");
      
      unsigned num_points = p1.size();
      unsigned dimensions = p1[0].size();
        
      vw::Matrix<double> A;
      vw::Matrix<double> B;
      
      A.set_size(num_points, dimensions);
      B.set_size(num_points, dimensions);

      for (unsigned r = 0; r < num_points; ++r) {
        for (unsigned c = 0; c < dimensions; ++c) {
          A(r,c) = p1[r][c];
          B(r,c) = p2[r][c];
        }
      }

      // Compute the least squares approximate fit
      vw::Matrix<double> H = transpose( pseudoinverse(A)*B );

      // Renormalize (to achieve 8-DOF)
      H /= H(2,2);
      return H;
    }
  };
 

  /// This fitting functor attempts to find an affine transformation
  /// (rotation, translation, scaling, and skewing -- 6 degrees of
  /// freedom) that transforms point p1 to match points p2.  This fit
  /// is optimal in a least squares sense.
  struct AffineFittingFunctor {
    typedef vw::Matrix<double,3,3> result_type;

    /// A similarity requires 3 pairs of data points to make a fit.
    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const { return 3; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1, 
                                   std::vector<ContainerT> const& p2) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(), 
                 vw::ArgumentErr() << "Cannot compute affine transformation.  p1 and p2 are not the same size." );
      VW_ASSERT( p1.size() != 0 && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute affine transformation.  Insufficient data.\n");
      

      Vector<double> y(p1.size()*2); 
      Vector<double> x(6);
      Matrix<double> A(p1.size()*2, 6);

      // Formulate a linear least squares problem to find the
      // components of the similarity matrix: 
      //       | s00 s01 s02 |
      //  S =  | s10 s11 s12 |
      //       | 0   0   1   |
      //
      for (unsigned i = 0; i < p1.size(); ++i) {
        A(i*2,0) = p1[i][0];
        A(i*2,1) = p1[i][1];
        A(i*2,4) = 1;
        A(i*2+1,2) = p1[i][0];
        A(i*2+1,3) = p1[i][1];
        A(i*2+1,5) = 1;
        
        y(i*2) = p2[i][0];
        y(i*2+1) = p2[i][1];
      }

      x = least_squares(A,y);

      Matrix<double> S(3,3);
      S.set_identity();
      S(0,0) = x(0);
      S(0,1) = x(1);
      S(1,0) = x(2);
      S(1,1) = x(3);
      S(0,2) = x(4);
      S(1,2) = x(5);

      return S;
    }
  };



  /// This fitting functor attempts to find an affine transformation
  /// (rotation, translation, scaling.
  struct SimilarityFittingFunctor {
    typedef vw::Matrix<double,3,3> result_type;

    /// A similarity requires 3 pairs of data points to make a fit.
    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const { return 3; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1, 
                                   std::vector<ContainerT> const& p2) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(), 
                 vw::ArgumentErr() << "Cannot compute affine transformation.  p1 and p2 are not the same size." );
      VW_ASSERT( p1.size() != 0 && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute affine transformation.  Insufficient data.\n");

      unsigned dimensions = p1[0].size()-1;

      // Compute the center of mass of each collection of points.
      MeanFunctor m(true);
      ContainerT mean1 = m(p1);
      ContainerT mean2 = m(p2);

      // Compute the scale factor between the points
      double dist1 = 0, dist2 = 0;
      for (unsigned i = 0; i < p1.size(); ++i) {
        dist1 += norm_2(p1[i]-mean1);
        dist2 += norm_2(p2[i]-mean2);
      }      
      dist1 /= p1.size();
      dist2 /= p2.size();
      double scale_factor = dist2/dist1;
          
      // Compute the rotation
      Matrix<double> H(dimensions, dimensions);
      for (unsigned i = 0; i < p1.size(); ++i) {
        Matrix<double> a(dimensions,1);
        Matrix<double> b(dimensions,1);
        for (unsigned d = 0; d < dimensions; ++d) {
          a(d,0) = p1[i][d]-mean1[d];
          b(d,0) = p2[i][d]-mean2[d];
        }
        H += a * transpose(b);
      }

      Matrix<double> U, VT;
      Vector<double> S;
      svd(H, U, S, VT);

      Matrix<double> R = transpose(VT)*transpose(U);
    
      // Compute the translation
      Vector<double> translation = subvector(mean2,0,2)-scale_factor*R*subvector(mean1,0,2);
  
      Matrix<double> result(3,3);
      submatrix(result,0,0,dimensions,dimensions) = scale_factor*R;
      for (unsigned i = 0; i < result.rows(); ++i) {
        result(i,dimensions) = translation(i);
      }
      result(dimensions,dimensions) = 1;
      return result;
    }
  };

  /// This fitting functor attempts to find an affine transformation
  /// (rotation, translation, scaling.
  struct TranslationRotationFittingFunctor {
    typedef vw::Matrix<double,3,3> result_type;

    /// A similarity requires 3 pairs of data points to make a fit.
    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const { return 3; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1, 
                                   std::vector<ContainerT> const& p2) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(), 
                 vw::ArgumentErr() << "Cannot compute affine transformation.  p1 and p2 are not the same size." );
      VW_ASSERT( p1.size() != 0 && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute affine transformation.  Insufficient data.\n");

      unsigned dimensions = p1[0].size()-1;

      // Compute the center of mass of each collection of points.
      MeanFunctor m(true);
      ContainerT mean1 = m(p1);
      ContainerT mean2 = m(p2);

      // Compute the rotation
      Matrix<double> H(dimensions, dimensions);
      for (unsigned i = 0; i < p1.size(); ++i) {
        Matrix<double> a(dimensions,1);
        Matrix<double> b(dimensions,1);
        for (unsigned d = 0; d < dimensions; ++d) {
          a(d,0) = p1[i][d]-mean1[d];
          b(d,0) = p2[i][d]-mean2[d];
        }
        H += a * transpose(b);
      }

      Matrix<double> U, VT;
      Vector<double> S;
      svd(H, U, S, VT);

      Matrix<double> R = transpose(VT)*transpose(U);
    
      // Compute the translation
      Vector<double> translation = subvector(mean2,0,2)-R*subvector(mean1,0,2);
  
      Matrix<double> result(3,3);
      submatrix(result,0,0,dimensions,dimensions) = R;
      for (unsigned i = 0; i < result.rows(); ++i) {
        result(i,dimensions) = translation(i);
      }
      result(dimensions,dimensions) = 1;
      return result;
    }
  };

}} // namespace vw::math

#endif // __MATH_GEOMETRY_H__
