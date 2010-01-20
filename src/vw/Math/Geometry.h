// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

  // This fitting functor attempts to find a 7 DoF matrix that charaterizes
  // the relationship of epipolar lines between to images. This implements
  // the algorithm defined on pg 281 of Multiview Geometry (aka Bible)
  struct FundamentalMatrixFittingFunctor7 {
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
      Matrix<double> n = null(A);

      // Nullspace seems to glitch sometimes. It should always by 2D here. -ZMM
      if ( n.cols() > 2 )
        n = submatrix(n,0,n.cols()-2,9,2);
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

  /// This fitting functor attempts to find a homography (8 degrees of
  /// freedom) that transforms point p1 to match points p2.  This fit
  /// is optimal in a least squares sense.
  struct HomographyFittingFunctor {
    typedef vw::Matrix<double> result_type;

    /// A homography requires at least 4 point matches to determine 8 unknowns
    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const { return 4; }

    // Applies a transform matrix to a list of points;
    std::vector<Vector<double> > apply_matrix( Matrix<double> const& m,
                 std::vector<Vector<double> > const& pts ) const {
      std::vector<Vector<double> > out;
      for ( unsigned i = 0; i < pts.size(); i++ ) {
        out.push_back( m*pts[i] );
      }
      return out;
    }

    /// Solve for Normalization Similarity Matrix used for noise
    /// rejection in DLT.
    vw::Matrix<double> NormSimilarity( std::vector<Vector<double> > const& pts ) const {
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

      std::vector<Vector<double> > pts_int = apply_matrix( translation, pts );

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

    vw::Matrix<double> BasicDLT( std::vector<Vector<double> > const& input,
                                 std::vector<Vector<double> > const& output )  const {
      VW_ASSERT( input.size() == 4 && output.size() == 4,
                 vw::ArgumentErr() << "DLT in this implementation expects to have only 4 inputs." );
      VW_ASSERT( input[0][input[0].size()-1] == 1, 
                 vw::ArgumentErr() << "Input data doesn't seem to be normalized.");
      VW_ASSERT( output[0][output[0].size()-1] == 1,
                 vw::ArgumentErr() << "Secondary input data doesn't seem to be normalized.");
      VW_ASSERT( input[0].size() == 3,
                 vw::ArgumentErr() << "Unfortunately at this time, BasicDLT only support homogeneous 2D vectors.");
      
      vw::Matrix<double,8,9> A;
      for ( unsigned i = 0; i < 4; i++ )
        for ( unsigned j = 0; j < 3; j++ ) {
          // Filling in -wi'*xi^T
          A(i,j+3) = -output[i][2]*input[i][j];
          // Filling in yi'*xi^T
          A(i,j+6) = output[i][1]*input[i][j];
          // Filling in wi'*xi^T
          A(i+4,j) = output[i][2]*input[i][j];
          // Filling in -xi'*xi^T
          A(i+4,j+6) = -output[i][0]*input[i][j];
        }

      Matrix<double> nullspace = null(A);
      nullspace /= nullspace(8,0);
      Matrix<double,3,3> H;
      for ( unsigned i = 0; i < 3; i++ )
        for ( unsigned j = 0; j < 3; j++ )
          H(i,j) = nullspace(i*3+j,0);
      return H;
    }

    /// Defining a Levenberg Marquard model that will be used to solve
    /// for a homography matrix.
    class HomographyModelLMA : public LeastSquaresModelBase<HomographyModelLMA> {
      Vector<double> m_measure; // Linear vector all from one side of an image.
      
      // Note that m_measure and result type will be 2*number of measures.
      // The scaling element is not kept and is assumed that it has
      // been normalized to 1.

    public:
      // What is returned by evaluating the functor. This is a vector
      // containing the results of all the measurements
      typedef Vector<double> result_type;
      // Defines the search space. In this case this is a flattened
      // version of the matrix.
      typedef Vector<double> domain_type;
      // The jacobian form. Don't really care
      typedef Matrix<double> jacobian_type;

      // Constructor
      inline HomographyModelLMA( Vector<double> const& measure ) : m_measure(measure) {}

      // Evaluator
      inline result_type operator()( domain_type const& x ) const {
        result_type output;
        output.set_size(m_measure.size());
        Matrix3x3 H;
        for ( unsigned i = 0; i < 3; i++ )
          for ( unsigned j = 0; j < 3; j++ )
            if ( i != 2 || j != 2 )
              H(i,j) = x( 3*i+j );
            else
              H(2,2) = 1;

        for ( unsigned i = 0; i < m_measure.size(); i+=2 ) {
          Vector3 input( m_measure[i], m_measure[i+1], 1 );
          Vector3 output_i = H*input;
          output_i /= output_i(2);
          output(i) = output_i(0);
          output(i+1) = output_i(1);
        }
        return output;
      }
    };

    /// Interface to solve for Homography Matrix. If the number of
    /// samples is 4, than a simple DLT is used. Otherwise it will use
    /// a minimization algorithm.
    template <class ContainerT>
    vw::Matrix<double> operator()( std::vector<ContainerT> const& p1,
           std::vector<ContainerT> const& p2,
           vw::Matrix<double> const& seed_input = vw::Matrix<double>() ) const {
      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "Cannot compute homography. p1 and p2 are not the same size." );
      VW_ASSERT( p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute homography. Insufficient data.");
      VW_ASSERT( p1[0].size() == 3,
                 vw::ArgumentErr() << "Cannot compute homography. Currently only support homogeneous 2D vectors." );
      VW_ASSERT( p1[0][2] == 1,
                 vw::ArgumentErr() << "Cannot compute homography. Vectors have not been normalized.");
      
      // Converting to a container that is used internally.
      std::vector<Vector<double> > input;
      std::vector<Vector<double> > output;
      for ( unsigned i = 0; i < p1.size(); i++ )
        input.push_back( Vector3( p1[i][0], p1[i][1], 1 ) );
      for ( unsigned i = 0; i < p2.size(); i++ ) 
        output.push_back( Vector3( p2[i][0], p2[i][1], 1 ) );

      unsigned num_points = p1.size();
      if ( num_points == min_elements_needed_for_fit(p1[0] ) ) {
        // Use DLT
        Matrix<double> S_in = NormSimilarity(input);
        Matrix<double> S_out = NormSimilarity(output);
        std::vector<Vector<double> > input_prime = apply_matrix( S_in, input );
        std::vector<Vector<double> > output_prime = apply_matrix( S_out, output );
        Matrix<double> H_prime = BasicDLT( input_prime, output_prime );
        Matrix<double> H = inverse(S_out)*H_prime*S_in;
        H /= H(2,2);
        return H;
      } else {
        // Levenberg Marquardt method for solving for homography.
        // - The error metric is x2 - norm(H*x1). 
        // - measure in x1 are fixed.
        // - Unfortunately at this time. Error in x1 are not dealt
        // with and I need more time to read up on other error
        // metrics.

        // Copying seed
        Matrix<double> seed_copy = seed_input;
        if ( seed_copy.cols() != 3 || seed_copy.rows() != 3 ) {
          std::vector<ContainerT> p1_small, p2_small;
          for ( unsigned i = 0; i < 4; i++ ) {
            p1_small.push_back( p1[i] );
            p2_small.push_back( p2[i] );
          }
          seed_copy = HomographyFittingFunctor()(p1_small, p2_small);
        }

        // Flatting input & output;
        Vector<double> input_flat( input.size()*2 );
        Vector<double> output_flat( output.size()*2 );
        for ( unsigned i = 0; i < input.size(); i++ ) {
          subvector(input_flat,i*2,2) = subvector(input[i],0,2);
          subvector(output_flat,i*2,2) = subvector(output[i],0,2);
        }

        HomographyModelLMA model( input_flat );

        // Flatting Homography matrix into 8-vector
        Vector<double> seed(8);
        for ( unsigned i = 0; i < 3; i++ )
          for ( unsigned j = 0; j < 3; j++ )
            if ( i != 2 || j != 2 )
              seed( i*3+j ) = seed_copy(i,j);

        int status = 0;
        Vector<double> result_flat = levenberg_marquardt( model, seed,
                                                          output_flat, status );
        
        // Unflatting result
        Matrix3x3 result;
        for ( unsigned i = 0; i < 3; i++ )
          for ( unsigned j = 0; j < 3; j++ )
            if ( i != 2 || j != 2 )
              result(i,j) = result_flat( 3*i+j );
            else
              result(2,2) = 1;

        return result;
      }
    }
  };

  /// This fitting functor attempts to find an affine transformation in RN
  /// (rotation, translation, scaling, and skewing) This fit
  /// is optimal in a least squares sense.
  template <int dim>
  struct AffineFittingFunctorN {
    typedef vw::Matrix<double,dim+1,dim+1> result_type;

    /// A affine transformation has dim*(dim+1) degrees of freedom so we need
    /// dim*(dim+1)/dim pairs of data points to make a fit.
    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const { return dim+1; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1, 
                                   std::vector<ContainerT> const& p2,
                                   vw::Matrix<double> const& seed_input = vw::Matrix<double>() ) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(), 
                 vw::ArgumentErr() << "Cannot compute affine transformation.  p1 and p2 are not the same size." );
      VW_ASSERT( p1.size() != 0 && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute affine transformation.  Insufficient data.\n");
      

      int dfree = dim * (dim + 1); // Number of parameters that we're solving for

      Vector<double> y(p1.size()*dim); 
      Vector<double> x(dfree);
      Matrix<double> A(p1.size()*dim, dfree);

      // Formulate a linear least squares problem to find the
      // components of the similarity matrix (shown for R2 below):
      //       | x0 x1 x2 |
      //  S =  | x3 x4 x5 |
      //       | x6 x7 x8 |
      //       |  0  0  1 |
      //
      // Least Squares problem:
      //
      // Ax = y
      //
      // x[0:8] represents the components of the similarity matrix, above
      // 
      // A is defined as follows:
      //
      // A(row i)     = | p1[i][0] p1[i][1]        1        0        0        0 |
      // A(row i + 1) = |        0        0        0 p1[i][0] p1[i][1]        1 |
      //
      // y is defined as follows:
      //
      // y(row i)     = | p2[i][0] |
      // y(row i + 1) = | p2[i][1] |
      
      for (unsigned i = 0; i < p1.size(); ++i) {
        for (unsigned j = 0; j < dim; ++j) {
          unsigned row = i*dim+j;
          for (unsigned l = 0; l < dim + 1; ++l) {
            A(row, j*(dim+1)+l) = l == dim ? 1 : p1[i][l];
          }
          y(row) = p2[i][j];
        }
      }
        
      x = least_squares(A,y);

      Matrix<double> S(dim+1,dim+1);
      S.set_identity();

      for (int i = 0; i < dfree; ++i) {
        S(i/(dim+1), i%(dim+1)) = x(i);
      }

      return S;
    }
  };

  /// This fitting functor attempts to find an affine transformation
  /// (rotation, translation, scaling, and skewing -- 6 degrees of
  /// freedom) that transforms point p1 to match points p2.  This fit
  /// is optimal in a least squares sense.
  typedef AffineFittingFunctorN<2> AffineFittingFunctor;

  /// This fitting functor attempts to find a similarity transformation
  /// (rotation, translation, scaling)
  template <int dim>
  struct SimilarityFittingFunctorN {
    typedef vw::Matrix<double,dim+1,dim+1> result_type;

    /// A similarity transformation requires 3 pairs of data points to make a fit.
    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const { return 3; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1, 
                                   std::vector<ContainerT> const& p2,
           vw::Matrix<double> const& seed_input = vw::Matrix<double>() ) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(), 
                 vw::ArgumentErr() << "Cannot compute similarity transformation.  p1 and p2 are not the same size." );
      VW_ASSERT( p1.size() != 0 && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute similarity transformation.  Insufficient data.\n");

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
      Matrix<double> H(dim, dim);
      for (unsigned i = 0; i < p1.size(); ++i) {
        Matrix<double> a(dim,1);
        Matrix<double> b(dim,1);
        for (unsigned d = 0; d < dim; ++d) {
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
      Vector<double> translation = subvector(mean2,0,dim)-scale_factor*R*subvector(mean1,0,dim);
  
      Matrix<double> result(dim+1,dim+1);
      submatrix(result,0,0,dim,dim) = scale_factor*R;
      for (unsigned i = 0; i < result.rows(); ++i) {
        result(i,dim) = translation(i);
      }
      result(dim,dim) = 1;
      return result;
    }
  };

  typedef SimilarityFittingFunctorN<2> SimilarityFittingFunctor;

  /// This fitting functor attempts to find a translation and rotation transformation
  template <int dim>
  struct TranslationRotationFittingFunctorN {
    typedef vw::Matrix<double,dim+1,dim+1> result_type;

    /// A transformation requires 2 pairs of data points to make a fit.
    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const { return 2; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1, 
                                   std::vector<ContainerT> const& p2,
           vw::Matrix<double> const& seed_input = vw::Matrix<double>() ) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(), 
                 vw::ArgumentErr() << "Cannot compute translation rotation transformation.  p1 and p2 are not the same size." );
      VW_ASSERT( p1.size() != 0 && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute translation rotation transformation.  Insufficient data.\n");

      // Compute the center of mass of each collection of points.
      MeanFunctor m(true);
      ContainerT mean1 = m(p1);
      ContainerT mean2 = m(p2);

      // Compute the rotation
      Matrix<double> H(dim, dim);
      for (unsigned i = 0; i < p1.size(); ++i) {
        Matrix<double> a(dim,1);
        Matrix<double> b(dim,1);
        for (unsigned d = 0; d < dim; ++d) {
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
      Vector<double> translation = subvector(mean2,0,dim)-R*subvector(mean1,0,dim);
  
      Matrix<double> result(dim+1,dim+1);
      submatrix(result,0,0,dim,dim) = R;
      for (unsigned i = 0; i < result.rows(); ++i) {
        result(i,dim) = translation(i);
      }
      result(dim,dim) = 1;
      return result;
    }
  };

  typedef TranslationRotationFittingFunctorN<2> TranslationRotationFittingFunctor;
  
  /// This fitting functor attempts to find a translation transformation
  template <int dim>
  struct TranslationFittingFunctorN {
    typedef vw::Matrix<double,dim+1,dim+1> result_type;

    /// A translation transformation needs 1 pair to make a fit
    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& example) const { return 1; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1, 
                                   std::vector<ContainerT> const& p2,
           vw::Matrix<double> const& seed_input = vw::Matrix<double>() ) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(), 
                 vw::ArgumentErr() << "Cannot compute translation transformation.  p1 and p2 are not the same size." );
      VW_ASSERT( p1.size() != 0 && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute translation transformation.  Insufficient data.\n");

      // Compute the center of mass of each collection of points.
      MeanFunctor m(true);
      ContainerT mean1 = m(p1);
      ContainerT mean2 = m(p2);
   
      // Compute the translation
      Vector<double> translation = subvector(mean2,0,dim)-subvector(mean1,0,dim);
  
      Matrix<double> result(dim+1,dim+1);
      result.set_identity();
      for (unsigned i = 0; i < dim; ++i) {
        result(i,dim) = translation(i);
      }
      return result;
    }
  };

  typedef TranslationFittingFunctorN<2> TranslationFittingFunctor;

}} // namespace vw::math

#endif // __MATH_GEOMETRY_H__
