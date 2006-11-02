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
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>

namespace vw { 
namespace math {

  /// This fitting functor attempts to find a homography (8 degrees of
  /// freedom) that transforms point p1 to match points p2.  This fit
  /// is optimal in a least squares sense.
  class HomographyFittingFunctor {
    bool m_homogeneous;
  public:
    typedef vw::Matrix<double> result_type;

    /// If you pass an example datum to this function, it will return
    /// the minimum number of putative matches needed to compute a fit.
    template <class ContainerT>
    int min_elements_needed_for_fit(ContainerT const& example) const {
      return example.size()*example.size();
    }
    
    /// If this functor is going to be applied to points in a
    /// projective space (i.e. homogeneous coordinates), you should
    /// set this flag to to true.  Otherwise, the functor will augment
    /// the points (make them homegeneous) when computing H.  In
    /// either case, the resulting Matrix will be in projective
    /// coordinates.
    HomographyFittingFunctor(bool homogeneous_points = false) 
      : m_homogeneous(homogeneous_points) {}
    
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
      
      int num_points = p1.size();
      int dimensions = p1[0].size();
        
      vw::Matrix<double> A;
      vw::Matrix<double> B;
      
      if (m_homogeneous) {
        A.set_size(num_points, dimensions);
        B.set_size(num_points, dimensions);
      } else {
        A.set_size(num_points, dimensions+1);
        B.set_size(num_points, dimensions+1);
        for (int i = 0; i < A.rows(); ++i) {
          A(i,dimensions) = 1.0;
          B(i,dimensions) = 1.0;
        }
      }

      for (int r = 0; r < num_points; ++r) {
        for (int c = 0; c < dimensions; ++c) {
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

  /// This fitting functor attempts to find a similarity (rotation,
  /// translation and scaling -- 5 degrees of freedom) that transforms
  /// point p1 to match points p2.  This fit is optimal in a least
  /// squares sense.
  class SimilarityFittingFunctor {
    bool m_homogeneous;
  public:
    typedef vw::Matrix<double> result_type;
    
    /// If you pass an example datum to this function, it will return
    /// the minimum number of putative matches needed to compute a fit.
    template <class ContainerT>
    int min_elements_needed_for_fit(ContainerT const& example) const {
      return example.size()*example.size();
    }
    /// If this functor is going to be applied to points in a
    /// projective space (i.e. homogeneous coordinates), you should
    /// set this flag to to true.  Otherwise, the functor will augment
    /// the points (make them homegeneous) when computing H.  In
    /// either case, the resulting Matrix will be in projective
    /// coordinates.
    SimilarityFittingFunctor(bool homogeneous_points = false) 
      : m_homogeneous(homogeneous_points) {}

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
      
      int num_points = p1.size();
      int dimensions = p1[0].size();
        
      vw::Matrix<double> A;
      vw::Matrix<double> B;
      
      if (m_homogeneous) {
        A.set_size(num_points, dimensions);
        B.set_size(num_points, dimensions-1);
      } else {
        A.set_size(num_points, dimensions+1);
        B.set_size(num_points, dimensions);
        for (int i = 0; i < A.rows(); ++i) {
          A(i,dimensions) = 1.0; 
        }
      }

      for (int r = 0; r < num_points; ++r) {
        for (int c = 0; c < dimensions; ++c) {
          A(r,c) = p1[r][c];
          if (c < B.cols())
            B(r,c) = p2[r][c];
        }
      }
 
      // Compute the least squares approximate fit for this similarity
      if (m_homogeneous) {
        vw::Matrix<double> H(dimensions, dimensions);
        fill(H,0.0);
        H(dimensions-1, dimensions-1) = 1.0;
        submatrix(H, 0, 0, dimensions-1, dimensions) = transpose( pseudoinverse(A)*B );
        return H;
      } else {
        vw::Matrix<double> H(dimensions+1, dimensions+1);
        fill(H,0.0);
        H(dimensions, dimensions) = 1.0;
        submatrix(H, 0, 0, dimensions, dimensions+1) = transpose( pseudoinverse(A)*B );
        return H;
      }
    }
  };

}} // namespace vw::math

#endif // __MATH_GEOMETRY_H__
