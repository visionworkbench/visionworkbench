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

#include <vw/Camera/CameraGeometry.h>
#include <vw/InterestPoint/InterestPoint.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Math/LevenbergMarquardt.h>

#include <cmath>

namespace vw {

template <>
Vector3 Convert2HomogenousVec2(vw::InterestPoint const& ip) {
  return Vector3(ip.x, ip.y, 1);
}

vw::Matrix<double>
SimilarityNormalizingFunctor::NormSimilarity(std::vector<Vector<double> > const& pts) {
  unsigned num_points = pts.size();
  unsigned dimension = pts[0].size();

  Vector<double> translation;
  translation.set_size(dimension-1);
  for (size_t i = 0; i < num_points; i++) {
    translation += subvector(pts[i], 0, dimension-1);
  }
  translation /= num_points;

  double scale = 0;
  for (size_t i = 0; i < num_points; i++) {
    Vector<double> delta = subvector(pts[i], 0, dimension-1) - translation;
    scale += norm_2(delta);
  }
  scale = num_points * std::sqrt(2.0) / scale;

  Matrix<double> s;
  s.set_size(dimension, dimension);
  s.set_identity();
  for (unsigned i = 0; i < dimension-1; i++) {
    s(i,i) *= scale;
    s(i,dimension-1) = -scale * translation[i];
  }
  return s;
}

namespace {

  // Simple linear solution
  Matrix<double> BasicDLT(std::vector<Vector<double> > const& input,
                          std::vector<Vector<double> > const& output) {
    VW_ASSERT(input[0].size() == 4 && output[0].size() == 3,
              vw::ArgumentErr() << "Camera Matrix requires Vector4 inputs and Vector3 outputs.");

    vw::Matrix<double,12,12> A;
    for (unsigned i = 0; i < 6; i++) { // Measure iterator
      for (unsigned j = 0; j < 4; j++) { // Measure's X elem iterator
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
    Matrix<double> U, VT;
    Vector<double> S;
    svd(A, U, S, VT);
    submatrix(p,0,0,1,4) = submatrix(VT, VT.rows()-1, 0, 1, 4);
    submatrix(p,1,0,1,4) = submatrix(VT, VT.rows()-1, 4, 1, 4);
    submatrix(p,2,0,1,4) = submatrix(VT, VT.rows()-1, 8, 1, 4);
    return p;
  }

  // Non-linear solution to an over-constrained problem
  class CameraMatrixModelLMA : public vw::LeastSquaresModelBase<CameraMatrixModelLMA> {
    std::vector<Vector<double> > m_world_input;
    std::vector<Vector<double> > m_image_output;

  public:
    typedef Vector<double> result_type;
    typedef Vector<double> domain_type;
    typedef Matrix<double> jacobian_type;

    CameraMatrixModelLMA(std::vector<Vector<double> > const& input,
                         std::vector<Vector<double> > const& output) :
      m_world_input(input), m_image_output(output) {}

    result_type operator()(domain_type const& x) const {
      result_type output;
      output.set_size(m_world_input.size());
      Matrix<double,3,4> P = unflatten(x);

      for (uint32 i = 0; i < m_world_input.size(); i++) {
        Vector3 reproj = P*m_world_input[i];
        reproj /= reproj[2];
        output[i] = norm_2(subvector(m_image_output[i],0,2) - subvector(reproj,0,2));
      }
      return output;
    }

    Vector<double> flatten(Matrix<double> const& input) const {
      Vector<double,12> output;
      for (uint8 i = 0; i < 3; i++) {
        for (uint8 j = 0; j < 4; j++) {
          output(4*i+j) = input(i,j);
        }
      }
      return output;
    }

    Matrix<double> unflatten(Vector<double> const& input) const {
      Matrix<double,3,4> output;
      for (uint8 i = 0; i < 3; i++) {
        for (uint8 j = 0; j < 4; j++) {
          output(i,j) = input(4*i+j);
        }
      }
      return output;
    }
  };

  // Skew symmetric matrix
  Matrix<double> skew_symm(Vector3 const& v) {
    Matrix3x3 skew;
    skew(0,1) = -v[2];
    skew(0,2) = v[1];
    skew(1,0) = v[2];
    skew(1,2) = -v[0];
    skew(2,0) = -v[1];
    skew(2,1) = v[0];
    return skew;
  }

  Vector4 OneSideTriangulation(Matrix<double> const& P2,
                               Vector3 const& v1,
                               Vector3 const& v2) {
    Matrix<double> A(4,4);
    select_row(A,0) = Vector4(-1, 0, v1[0], 0);
    select_row(A,1) = Vector4(0, -1, v1[1], 0);
    select_row(A,2) = v2[0]*select_row(P2,2) - select_row(P2,0);
    select_row(A,3) = v2[1]*select_row(P2,2) - select_row(P2,1);

    Matrix<double> U, VT;
    Vector<double> S;
    svd(A, U, S, VT);

    Vector<double> solution = select_row(VT, 3);
    solution /= solution[3];
    return solution;
  }

  class ProjectiveModelLMA : public vw::LeastSquaresModelBase<ProjectiveModelLMA> {
  public:
    typedef Vector<double> result_type;
    typedef Vector<double> domain_type;
    typedef Matrix<double> jacobian_type;

    result_type operator()(domain_type const& x) const {
      unsigned number_of_measures = (x.size()-12)/3;
      Vector<double> output(number_of_measures*4);

      // Create second camera
      Matrix<double> P2(3,4);
      select_row(P2,0) = subvector(x,0,4);
      select_row(P2,1) = subvector(x,4,4);
      select_row(P2,2) = subvector(x,8,4);

      // Calculate reprojection error
      for (unsigned i = 0, j = 12; i < number_of_measures; i++, j+= 3) {
        Vector2 reprojection1(x[j]/x[j+2], x[j+1]/x[j+2]);
        subvector(output, 4*i, 2) = reprojection1;
        Vector3 reprojection2 = P2*Vector4(x[j], x[j+1], x[j+2], 1);
        reprojection2 /= reprojection2[2];
        subvector(output, 4*i+2, 2) = subvector(reprojection2, 0, 2);
      }

      return output;
    }
  };

} // end anonymous namespace

Matrix<double>
CameraMatrixFittingFunctor::fit(std::vector<Vector<double> > const& input_pts,
                                std::vector<Vector<double> > const& output_pts,
                                vw::Matrix<double> const& seed_input) const {

  std::vector<Vector<double> > input = input_pts;
  std::vector<Vector<double> > output = output_pts;

  // Normalizing
  Matrix<double> S_in = NormSimilarity(input);
  Matrix<double> S_out = NormSimilarity(output);
  for (size_t i = 0; i < input.size(); i++)
    input[i] = S_in * input[i];
  for (size_t i = 0; i < output.size(); i++)
    output[i] = S_out * output[i];

  Matrix<double> p;
  if (input.size() == 6) {
    // Use DLT
    Matrix<double> P_prime = BasicDLT(input, output);
    p = inverse(S_out)*P_prime*S_in;
  } else {
    Matrix<double> seed_matrix = seed_input;

    // Perform DLT if needed
    if (seed_matrix.cols() == 0)
      seed_matrix = BasicDLT(input, output);

    // Iterative solution
    CameraMatrixModelLMA model(input, output);
    Vector<double> seed = model.flatten(seed_matrix);
    int status = 0;
    Vector<double> objective;
    objective.set_size(input.size());
    Vector<double> result = levenberg_marquardt(model, seed, objective, status);
    seed_matrix = model.unflatten(result);

    // Denormalization
    p = inverse(S_out)*seed_matrix*S_in;
  }

  p /= p(2,3);
  return p;
}

Matrix<double>
FundamentalMatrix8PFittingFunctor::fit(std::vector<Vector<double> > const& input_pts,
                                       std::vector<Vector<double> > const& output_pts) const {

  std::vector<Vector<double> > input = input_pts;
  std::vector<Vector<double> > output = output_pts;

  // Normalizing
  Matrix<double> S_in  = NormSimilarity(input);
  Matrix<double> S_out = NormSimilarity(output);
  for (size_t i = 0; i < input.size(); i++)
    input[i] = S_in * input[i];
  for (size_t i = 0; i < output.size(); i++)
    output[i] = S_out * output[i];

  // Constructing A
  Matrix<double> A(input.size(), 9);
  for (size_t i = 0; i < input.size(); i++) {
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

  VW_ASSERT(vw::rank(A) >= 8, MathErr() << "Measurements produce rank deficient A.");

  // Pulling singular vector of smallest singular value of A
  Matrix<double> U, VT;
  Vector<double> S;
  svd(A, U, S, VT);
  Matrix3x3 F;
  int idx = 0;
  for (Matrix<double,3,3>::iterator it = F.begin(); it != F.end(); ++it) {
    (*it) = VT(VT.rows()-1, idx);
    idx++;
  }

  // Constraint Enforcement
  svd(F, U, S, VT);
  S[2] = 0;
  F = U*diagonal_matrix(S)*VT;

  // Denormalizing
  return transpose(S_out)*F*S_in;
}

Matrix<double>
FundamentalMatrixMLFittingFunctor::fit(std::vector<Vector<double> > const& input,
                                       std::vector<Vector<double> > const& output,
                                       vw::Matrix<double> const& seed_input) const {

  // Getting epipole
  Matrix<double> U, VT;
  Vector<double> S;
  svd(seed_input, U, S, VT);
  Vector<double> epipole_prime = select_col(U, 2);

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
  for (unsigned i = 0, j = 12; i < input.size(); i++, j+=3) {
    Vector4 projection = OneSideTriangulation(P2, input[i], output[i]);
    projection /= projection[3];
    subvector(seed,j,3) = subvector(projection,0,3);
  }

  // Setting Objective
  Vector<double> objective(4*input.size());
  for (unsigned i = 0; i < input.size(); i++) {
    subvector(objective, 4*i, 2) = subvector(input[i], 0, 2);
    subvector(objective, 4*i+2, 2) = subvector(output[i], 0, 2);
  }

  int status = 0;
  Vector<double> result = levenberg_marquardt(model, seed, objective, status);

  select_row(P2,0) = subvector(result,0,4);
  select_row(P2,1) = subvector(result,4,4);
  select_row(P2,2) = subvector(result,8,4);

  Matrix<double> F = skew_symm(select_col(P2,3))*submatrix(P2,0,0,3,3);
  return F;
}

} // namespace vw
