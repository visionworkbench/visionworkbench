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


/// \file Geometry.h
///
/// Assorted useful geometric routines and functors.
///
#ifndef __MATH_GEOMETRY_H__
#define __MATH_GEOMETRY_H__

#include <vw/Core/Exception.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/Statistics.h>
#include <vw/Math/LevenbergMarquardt.h>

#include <stddef.h>
#include <vector>
#include <algorithm>

namespace vw {
namespace math {

  /// Normalizes a degree longitude value into either the -180 to 180 range (default) or the 0-360 range.
  double normalize_longitude(double lon, bool center_on_zero=true);

  /// Computes the absolute distance between to measurements in degrees, accounting for wraparound.
  double degree_diff(double d1, double d2);


  /// This fitting functor attempts to find a homography (8 degrees of
  /// freedom) that transforms point p1 to match points p2.  This fit
  /// is optimal in a least squares sense.
  struct HomographyFittingFunctor {
    typedef vw::Matrix<double> result_type;

    /// A homography requires at least 4 point matches to determine 8 unknowns
    template <class ContainerT>
    size_t min_elements_needed_for_fit(ContainerT const& /*example*/) const { return 4; }

    // Applies a transform matrix to a list of points;
    template <int D>
    struct ApplyMatrixFunctor {
      Matrix<double,D,D> m;
      template <class MatrixT>
      ApplyMatrixFunctor( MatrixBase<MatrixT> const& mat ) : m( mat.impl() ) {}

      Vector<double,D> operator()( Vector<double,D> const& v ) {
        return m * v;
      }
    };

    /// Solve for Normalization Similarity Matrix used for noise rej.
    vw::Matrix3x3 NormSimilarity( std::vector<Vector3> const& pts ) const;
    vw::Matrix3x3 BasicDLT( std::vector<Vector3 > const& input,
                            std::vector<Vector3 > const& output )  const;

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
      HomographyModelLMA( Vector<double> const& measure ) : m_measure(measure) {}

      // Evaluator
      inline result_type operator()( domain_type const& x ) const {
        result_type output;
        output.set_size(m_measure.size());
        Matrix3x3 H;
        H(2,2) = 1;
        VectorProxy<double,8>( H.data() ) = x;

        for ( size_t i = 0; i < m_measure.size(); i+=2 ) {
          Vector3 output_i = H*Vector3( m_measure[i], m_measure[i+1], 1 );
          output_i /= output_i(2);
          subvector(output,i,2) = subvector(output_i,0,2);
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

      // Converting to a container that is used internally.
      std::vector<Vector3 > input(p1.size()), output(p2.size());
      for ( size_t i = 0; i < p1.size(); i++ ) {
        input[i] = Vector3( p1[i][0], p1[i][1], p1[i][2] );
      }
      for ( size_t i = 0; i < p2.size(); i++ ) {
        output[i] = Vector3( p2[i][0], p2[i][1], p2[i][2] );
      }

      if ( p1.size() == min_elements_needed_for_fit(p1[0] ) ) {
        // Use DLT
        Matrix3x3 S_in = NormSimilarity(input);
        Matrix3x3 S_out = NormSimilarity(output);
        std::transform( input.begin(), input.end(), input.begin(),
                        ApplyMatrixFunctor<3>( S_in ) );
        std::transform( output.begin(), output.end(), output.begin(),
                        ApplyMatrixFunctor<3>( S_out ) );
        Matrix3x3 H_prime = BasicDLT( input, output );
        Matrix3x3 H = inverse(S_out)*H_prime*S_in;
        H /= H(2,2);
        return H;
      } else {
        // Levenberg Marquardt method for solving for homography.
        // - The error metric is x2 - norm(H*x1).
        // - measure in x1 are fixed.
        // - Unfortunately at this time. Error in x1 are not dealt
        // with and I need more time to read up on other error
        // metrics.

        // Determine if we have a seed .. if not make one.
        Matrix3x3 seed_copy;
        if ( seed_input == Matrix<double>() ) {
          std::vector<ContainerT> p1_small(4), p2_small(4);
          for ( size_t i = 0; i < 4; i++ ) {
            p1_small[i] = p1[i];
            p2_small[i] = p2[i];
          }
          seed_copy = HomographyFittingFunctor()(p1_small, p2_small);
        } else {
          seed_copy = seed_input;
        }

        // Flatting input & output;
        Vector<double> input_flat( input.size()*2 );
        Vector<double> output_flat( output.size()*2 );
        for ( size_t i = 0; i < input.size(); i++ ) {
          subvector(input_flat,i*2,2) = subvector(input[i],0,2);
          subvector(output_flat,i*2,2) = subvector(output[i],0,2);
        }

        HomographyModelLMA model( input_flat );

        // Flatting Homography matrix into 8-vector
        Vector<double> seed(8);
        std::copy( seed_copy.data(), seed_copy.data()+8,
                   &seed(0) );

        int status = 0;
        VectorProxy<double,8>( seed_copy.data() ) =
          levenberg_marquardt( model, seed,
                               output_flat, status );

        return seed_copy;
      }
    }
  };

  /// This fitting functor attempts to find an affine transformation in RN
  /// (rotation, translation, scaling, and skewing) This fit
  /// is optimal in a least squares sense.
  template <size_t dim>
  struct AffineFittingFunctorN {
    typedef vw::Matrix<double,dim+1,dim+1> result_type;

    /// A affine transformation has dim*(dim+1) degrees of freedom so we need
    /// dim*(dim+1)/dim pairs of data points to make a fit.
    template <class ContainerT>
    size_t min_elements_needed_for_fit(ContainerT const& /*example*/) const { return dim+1; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1,
                                   std::vector<ContainerT> const& p2,
                                   vw::Matrix<double> const& /*seed_input*/ = vw::Matrix<double>() ) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "Cannot compute affine transformation.  p1 and p2 are not the same size." );
      VW_ASSERT( !p1.empty() && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute affine transformation.  Insufficient data.\n");


      size_t dfree = dim * (dim + 1); // Number of parameters that we're solving for

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

      for (size_t i = 0; i < p1.size(); ++i) {
        for (size_t j = 0; j < dim; ++j) {
          size_t row = i*dim+j;
          for (size_t l = 0; l < dim + 1; ++l) {
            A(row, j*(dim+1)+l) = l == dim ? 1 : p1[i][l];
          }
          y(row) = p2[i][j];
        }
      }

      x = least_squares(A,y);

      Matrix<double> S(dim+1,dim+1);
      S.set_identity();

      for (size_t i = 0; i < dfree; ++i) {
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
  template <size_t dim>
  struct SimilarityFittingFunctorN {
    typedef vw::Matrix<double,dim+1,dim+1> result_type;

    /// A similarity transformation requires 3 pairs of data points to make a fit.
    template <class ContainerT>
    size_t min_elements_needed_for_fit(ContainerT const& /*example*/) const { return 3; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1,
                                   std::vector<ContainerT> const& p2,
                                   vw::Matrix<double> const& /*seed_input*/ = vw::Matrix<double>() ) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "Cannot compute similarity transformation.  p1 and p2 are not the same size." );
      VW_ASSERT( !p1.empty() && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute similarity transformation.  Insufficient data.\n");

      // Compute the center of mass of each collection of points.
      MeanFunctor m(true);
      ContainerT mean1 = m(p1);
      ContainerT mean2 = m(p2);

      // Compute the scale factor between the points
      double dist1 = 0, dist2 = 0;
      for (size_t i = 0; i < p1.size(); ++i) {
        dist1 += norm_2(p1[i]-mean1);
        dist2 += norm_2(p2[i]-mean2);
      }
      dist1 /= p1.size();
      dist2 /= p2.size();
      double scale_factor = dist2/dist1;

      // Compute the rotation
      Matrix<double> H(dim, dim);
      for (size_t i = 0; i < p1.size(); ++i) {
        Matrix<double> a(dim,1);
        Matrix<double> b(dim,1);
        for (size_t d = 0; d < dim; ++d) {
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
      for (size_t i = 0; i < translation.size(); ++i) {
        result(i,dim) = translation(i);
      }
      result(dim,dim) = 1;
      return result;
    }
  };

  typedef SimilarityFittingFunctorN<2> SimilarityFittingFunctor;

  /// This fitting functor attempts to find a translation and scale transformation
  /// (each coordinate has its own scale factor).
  template <size_t dim>
  struct TranslationScaleFittingFunctorN {
    typedef vw::Matrix<double,dim+1,dim+1> result_type;

    /// A translation and scale transformation requires 2 pairs of data points to make a fit.
    template <class ContainerT>
    size_t min_elements_needed_for_fit(ContainerT const& /*example*/) const { return 2; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as
    /// well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1,
                                   std::vector<ContainerT> const& p2,
                                   vw::Matrix<double> const& /*seed_input*/ = vw::Matrix<double>() ) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "Cannot compute translation and scale transformation. p1 and p2 are not the same size." );
      VW_ASSERT( !p1.empty() && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute translation and scale transformation. Insufficient data.\n");

      // Compute the center of mass of each collection of points.
      MeanFunctor m(true);
      ContainerT mean1 = m(p1);
      ContainerT mean2 = m(p2);

      // Compute the scale factor between the points
      Vector<double, dim> scale;
      for (size_t d = 0; d < dim; d++) {
        double dist1 = 0.0, dist2 = 0.0;
        for (size_t i = 0; i < p1.size(); ++i) {
          dist1 += std::abs(p1[i][d]-mean1[d]);
          dist2 += std::abs(p2[i][d]-mean2[d]);
        }
        dist1 /= p1.size();
        dist2 /= p2.size();
        if (dist1 != 0)
          scale[d] = dist2/dist1;
        else
          scale[d] = 1.0; // Can't find the scale, assume it is 1.
      }

      // Compute the translation
      Vector<double, dim> translation;
      for (size_t d = 0; d < dim; d++) {
        translation[d] = mean2[d] - scale[d]*mean1[d];
      }

      // Assemble the matrix of the transform
      Matrix<double> result(dim+1,dim+1);
      result.set_identity();
      for (size_t d = 0; d < dim; d++) {
        result(d, d) = scale[d];
      }
      for (size_t d = 0; d < translation.size(); ++d) {
        result(d,dim) = translation(d);
      }
      result(dim,dim) = 1.0;
      return result;
    }
  };

  typedef TranslationScaleFittingFunctorN<2> TranslationScaleFittingFunctor;

  /// This fitting functor attempts to find a translation and rotation transformation
  template <size_t dim>
  struct TranslationRotationFittingFunctorN {
    typedef vw::Matrix<double,dim+1,dim+1> result_type;

    /// A transformation requires 2 pairs of data points to make a fit.
    template <class ContainerT>
    size_t min_elements_needed_for_fit(ContainerT const& /*example*/) const { return 2; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1,
                                   std::vector<ContainerT> const& p2,
      vw::Matrix<double> const& /*seed_input*/ = vw::Matrix<double>() ) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "Cannot compute translation rotation transformation.  p1 and p2 are not the same size." );
      VW_ASSERT( !p1.empty() && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute translation rotation transformation.  Insufficient data.\n");

      // Compute the center of mass of each collection of points.
      MeanFunctor m(true);
      ContainerT mean1 = m(p1);
      ContainerT mean2 = m(p2);

      // Compute the rotation
      Matrix<double> H(dim, dim);
      for (size_t i = 0; i < p1.size(); ++i) {
        Matrix<double> a(dim,1);
        Matrix<double> b(dim,1);
        for (size_t d = 0; d < dim; ++d) {
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
      for (size_t i = 0; i < translation.size(); ++i) {
        result(i,dim) = translation(i);
      }
      result(dim,dim) = 1;
      return result;
    }
  };

  typedef TranslationRotationFittingFunctorN<2> TranslationRotationFittingFunctor;

  /// This fitting functor attempts to find a translation transformation
  template <size_t dim>
  struct TranslationFittingFunctorN {
    typedef vw::Matrix<double,dim+1,dim+1> result_type;

    /// A translation transformation needs 1 pair to make a fit
    template <class ContainerT>
    size_t min_elements_needed_for_fit(ContainerT const& /*example*/) const { return 1; }

    /// This function can match points in any container that supports
    /// the size() and operator[] methods.  The container is usually a
    /// vw::Vector<>, but you could substitute other classes here as well.
    template <class ContainerT>
    vw::Matrix<double> operator() (std::vector<ContainerT> const& p1,
                                   std::vector<ContainerT> const& p2,
      vw::Matrix<double> const& /*seed_input*/ = vw::Matrix<double>() ) const {

      // check consistency
      VW_ASSERT( p1.size() == p2.size(),
                 vw::ArgumentErr() << "Cannot compute translation transformation.  p1 and p2 are not the same size." );
      VW_ASSERT( !p1.empty() && p1.size() >= min_elements_needed_for_fit(p1[0]),
                 vw::ArgumentErr() << "Cannot compute translation transformation.  Insufficient data.\n");

      // Compute the center of mass of each collection of points.
      MeanFunctor m(true);
      ContainerT mean1 = m(p1);
      ContainerT mean2 = m(p2);

      // Compute the translation
      Vector<double> translation = subvector(mean2,0,dim)-subvector(mean1,0,dim);

      Matrix<double> result(dim+1,dim+1);
      result.set_identity();
      for (size_t i = 0; i < translation.size(); ++i) {
        result(i,dim) = translation(i);
      }
      return result;
    }
  };

  typedef TranslationFittingFunctorN<2> TranslationFittingFunctor;

}} // namespace vw::math

#endif // __MATH_GEOMETRY_H__
