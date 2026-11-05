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
///
#ifndef __CAMERA_CAMERA_GEOMETRY_H__
#define __CAMERA_CAMERA_GEOMETRY_H__

#include <vw/Core/Exception.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

#include <vector>
#include <algorithm>

namespace vw {

// Forward Declare
struct InterestPoint;

  // Transform Functors
  template <class VectorT>
  inline Vector3 Convert2HomogenousVec2(VectorT const& v) {
    if (v.size() == 2)
      return Vector3(v[0], v[1], 1);
    return Vector3(v[0], v[1], v[2]);
  }
  template <>
  Vector3 Convert2HomogenousVec2(vw::InterestPoint const& ip);

  // Measurement Prep Functors
  struct SimilarityNormalizingFunctor {
    // Solve for normalizing similarity matrix
    static vw::Matrix<double>
    NormSimilarity(std::vector<Vector<double> > const& pts);
  };

  // Camera Matrix Functor
  //
  // Solves for 3x4 matrix which describes the mapping of a pinhole
  // camera from 3D points in the world to 2D points in an image.

  struct CameraMatrixErrorMetric {
    template <class RelationT, class ContainerT>
    double operator()(RelationT const& H,
                      ContainerT const& p1,
                      ContainerT const& p2) const {
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

    // Fit camera matrix from 4D world coordinates and 3D homogeneous image coordinates.
    // Implementation is compiled into CameraGeometry.cc.
    Matrix<double> fit(std::vector<Vector<double> > const& input,
                       std::vector<Vector<double> > const& output,
                       vw::Matrix<double> const& seed_input = vw::Matrix<double>()) const;

    // Interface for solving for Camera Matrix (switches between DLT
    // and Gold Standard)
    template <class ContainerT>
    result_type
    operator()(std::vector<ContainerT> const& p1,
               std::vector<ContainerT> const& p2,
               vw::Matrix<double> const& seed_input = vw::Matrix<double>()) const {
      VW_ASSERT(p1.size() == p2.size(),
                vw::ArgumentErr() << "Camera Matrix requires equal number of input and output measures.");
      VW_ASSERT(p1.size() >= min_elements_needed_for_fit(p1[0]),
                vw::ArgumentErr() << "Cannot compute camera matrix. Insufficient measurements.");
      VW_ASSERT(p1[0].size() == 4 && p2[0].size() == 3,
                vw::ArgumentErr() << "Incorrect input format: p1 needs to Vector4 and p2 need to be Vector3.");

      // Converting to internal format
      std::vector<Vector<double> > input, output;
      input.reserve(p1.size());
      output.reserve(p2.size());
      for (size_t i = 0; i < p1.size(); i++) {
        input.push_back(Vector4(p1[i][0], p1[i][1], p1[i][2], p1[i][3]));
        output.push_back(Vector3(p2[i][0], p2[i][1], p2[i][2]));
      }

      return fit(input, output, seed_input);
    }
  };

  // Fundamental Matrix Functor
  //
  // Solves for 3x3 matrix which relates stereo images. Fundamental
  // matrix ties points in first image to lines in the second image. The
  // lines are defined by the epipole (which is the location of the first
  // image in the perspective of the second camera (and vice versa)).

  struct FundamentalMatrixSampsonErrorMetric {
    template <class RelationT, class ContainerT>
    double operator()(RelationT const& F,
                      ContainerT const& p1,
                      ContainerT const& p2) const {
      return fabs(transpose(p2)*F*p1) + fabs(transpose(p1)*transpose(F)*p2);
    }
  };

  struct FundamentalMatrixDistanceErrorMetric {
    template <class RelationT, class ContainerT>
    double operator()(RelationT const& F,
                      ContainerT const& p1,
                      ContainerT const& p2) const {
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

    /// Fit fundamental matrix using normalized 8-point algorithm.
    /// Implementation is compiled into CameraGeometry.cc.
    Matrix<double> fit(std::vector<Vector<double> > const& input,
                       std::vector<Vector<double> > const& output) const;

    /// Solve for the fundamental matrix F using two sets of observations of a set
    /// of points from the same camera at different positions.
    template <class ContainerT>
    result_type
    operator()(std::vector<ContainerT> const& p1,
               std::vector<ContainerT> const& p2,
               vw::Matrix<double> const& /*seed_input*/ = vw::Matrix<double>()) const {
      VW_ASSERT(p1.size() == p2.size(),
                vw::ArgumentErr() << "Fundamental Matrix requires equal number of input and output measures.");

      std::vector<Vector<double> > input(p1.size()), output(p2.size());
      std::transform(p1.begin(), p1.end(), input.begin(),  Convert2HomogenousVec2<ContainerT>);
      std::transform(p2.begin(), p2.end(), output.begin(), Convert2HomogenousVec2<ContainerT>);

      return fit(input, output);
    }
  };

  // Fundamental Matrix solver using the Maximum Likelihood method
  // Page 285 or Algorithm 11.3 in Multiple View Geometry
  struct FundamentalMatrixMLFittingFunctor {
    typedef vw::Matrix<double> result_type;

    template <class ContainerT>
    unsigned min_elements_needed_for_fit(ContainerT const& /*example*/) const { return 8; }

    /// Fit fundamental matrix using maximum likelihood optimization.
    /// Implementation is compiled into CameraGeometry.cc.
    Matrix<double> fit(std::vector<Vector<double> > const& input,
                       std::vector<Vector<double> > const& output,
                       vw::Matrix<double> const& seed_input = vw::Matrix<double>()) const;

    /// Solve for the fundamental matrix F using two sets of observation of a set
    /// of points from the same camera at different positions.
    template <class ContainerT>
    result_type
    operator()(std::vector<ContainerT> const& p1,
               std::vector<ContainerT> const& p2,
               vw::Matrix<double> const& seed_input = vw::Matrix<double>()) const {
      VW_ASSERT(p1.size() == p2.size(),
                vw::ArgumentErr() << "Fundamental Matrix requires equal number of input and output measures.");

      if (p1.size() == 8)
        return FundamentalMatrix8PFittingFunctor()(p1, p2);

      std::vector<Vector<double> > input(p1.size()), output(p2.size());
      std::transform(p1.begin(), p1.end(), input.begin(),
                     Convert2HomogenousVec2<ContainerT>);
      std::transform(p2.begin(), p2.end(), output.begin(),
                     Convert2HomogenousVec2<ContainerT>);

      return fit(input, output, seed_input);
    }
  };

} // namespace vw

#endif // __CAMERA_CAMERA_GEOMETRY_H__
