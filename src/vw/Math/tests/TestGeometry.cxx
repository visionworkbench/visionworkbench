// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>
#include <test/Helpers.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Geometry.h>
#include <vw/Math/EulerAngles.h>

using namespace vw;
using namespace vw::math;

TEST(Geometry, HomographyFittingFunctor) {
  static double A_data[] = {
    0.0153, 0.9318, 1,
    0.6721, 0.6813, 1,
    0.3046, 0.6822, 1,
    0.8600, 0.7468, 1,
    0.5252, 0.8381, 1,
    0.7095, 0.1897, 1,
    0.6979, 0.8537, 1,
    0.4186, 0.2026, 1,
    0.8318, 0.4289, 1,
    0.5417, 0.3784, 1 };

  vw::MatrixProxy<double,10,3> A(A_data);

  static double H_data[] = {
    1.1567,  0.5916,  0.5557,
    0.2814,  1.0851,  0.0225,
    0.7388,  0.9278,  1.0000};

  vw::MatrixProxy<double,3,3> H(H_data);

  // and apply it to some points
  vw::Matrix<double> B = transpose(H*transpose(A));

  // Normalizing B (doesn't work otherwise)
  for (unsigned i = 0; i < B.rows(); ++i)
    submatrix(B,i,0,1,3) /= B(i,2);

  std::vector<Vector3> p1, p2, p1_s, p2_s;
  for (unsigned i = 0; i < A.rows(); ++i) {
    p1.push_back(select_row(A,i));
    p2.push_back(select_row(B,i));
    if ( i < 4 ) {
      p1_s.push_back(select_row(A,i));
      p2_s.push_back(select_row(B,i));
    }
  }

  // DLT version (there's a loss of precision in complete_svd)
  EXPECT_MATRIX_NEAR(H, HomographyFittingFunctor()(p1_s,p2_s), 4.9e-14 );
  // DLT & Levenberg Marquardt Version (it's better with noise and conflict)
  EXPECT_MATRIX_NEAR(H, HomographyFittingFunctor()(p1,p2),  2e-14 );
}

TEST(Geometry, AffinityFittingFunctor) {

  static double A_data[] = {
    0.0153, 0.9318, 1,
    0.6721, 0.6813, 1,
    0.3046, 0.6822, 1,
    0.8600, 0.7468, 1,
    0.5252, 0.8381, 1,
    0.7095, 0.1897, 1,
    0.6979, 0.8537, 1,
    0.4186, 0.2026, 1,
    0.8318, 0.4289, 1,
    0.5417, 0.3784, 1 };
  vw::MatrixProxy<double,10,3> A(A_data);

  static double S_data[] = {
    1.1567,  0.5916,  0.5557,
    0.2814,  1.0851,  0.0225,
    0,  0,  1.0000};

  vw::MatrixProxy<double,3,3> S(S_data);

  vw::Matrix<double> B = transpose(S*transpose(A));
  std::vector<Vector3> p1, p2;
  for (unsigned i = 0; i < A.rows(); ++i) {
    p1.push_back(select_row(A,i));
    p2.push_back(select_row(B,i));
  }

  EXPECT_MATRIX_NEAR( S, AffineFittingFunctor()(p1,p2), 1e-14 );
}

TEST(Geometry, AffinityFittingFunctorN) {

  static double A_data[] = {
    0.0153, 0.9318, 0.2345, 1,
    0.6721, 0.6813, 0.1234, 1,
    0.3046, 0.6822, 0.6543, 1,
    0.8600, 0.7468, 0.1652, 1,
    0.5252, 0.8381, 0.5089, 1,
    0.7095, 0.1897, 0.2258, 1,
    0.6979, 0.8537, 0.0461, 1,
    0.4186, 0.2026, 0.1419, 1,
    0.8318, 0.4289, 0.6408, 1,
    0.3475, 0.3354, 0.2868, 1,
    0.7025, 0.2123, 0.2339, 1,
    0.4944, 0.9446, 0.1233, 1,
    0.7193, 0.6551, 0.3655, 1,
    0.5491, 0.4076, 0.1014, 1 };
  vw::MatrixProxy<double,14,4> A(A_data);

  static double S_data[] = {
    1.1567,  0.5916,  0.5557, -0.7443,
    0.2814,  1.0851,  0.0225, 0.3333,
    0.9074,  0.4688,  0.7468, 0.8321,
         0,       0,       0, 1.0000};

  vw::MatrixProxy<double,4,4> S(S_data);

  vw::Matrix<double> B = transpose(S*transpose(A));
  std::vector<Vector4> p1, p2;
  for (unsigned i = 0; i < A.rows(); ++i) {
    p1.push_back(select_row(A,i));
    p2.push_back(select_row(B,i));
  }

  EXPECT_MATRIX_NEAR( S, AffineFittingFunctorN<3>()(p1,p2), 1e-14 );
}

TEST(Geometry, SimilarityFittingFunctorN) {

  static double A_data[] = {
    0.0153, 0.9318, 0.2345, 1,
    0.6721, 0.6813, 0.1234, 1,
    0.3046, 0.6822, 0.6543, 1,
    0.8600, 0.7468, 0.1652, 1,
    0.5252, 0.8381, 0.5089, 1,
    0.7095, 0.1897, 0.2258, 1,
    0.6979, 0.8537, 0.0461, 1,
    0.4186, 0.2026, 0.1419, 1,
    0.8318, 0.4289, 0.6408, 1,
    0.3475, 0.3354, 0.2868, 1,
    0.7025, 0.2123, 0.2339, 1,
    0.4944, 0.9446, 0.1233, 1,
    0.7193, 0.6551, 0.3655, 1,
    0.5491, 0.4076, 0.1014, 1 };
  vw::MatrixProxy<double,14,4> A(A_data);

  double s = 2.23;
  double xrot = 0.325;
  double yrot = 0.122;
  double zrot = 0.644;
  double dx = 5.6;
  double dy = 8.3;
  double dz = 2.1;

  Matrix4x4 S;

  submatrix(S, 0, 0, 3, 3) = (s * identity_matrix(3)) * euler_to_rotation_matrix(xrot, yrot, zrot, "xyz");
  select_col(S, 3) = Vector4(dx, dy, dz, 1);

  vw::Matrix<double> B = transpose(S*transpose(A));
  std::vector<Vector4> p1, p2;
  for (unsigned i = 0; i < A.rows(); ++i) {
    p1.push_back(select_row(A,i));
    p2.push_back(select_row(B,i));
  }

  EXPECT_MATRIX_NEAR( S, SimilarityFittingFunctorN<3>()(p1,p2), 1e-14 );
}

TEST(Geometry, TranslationRotationFittingFunctorN) {

  static double A_data[] = {
    0.0153, 0.9318, 0.2345, 1,
    0.6721, 0.6813, 0.1234, 1,
    0.3046, 0.6822, 0.6543, 1,
    0.8600, 0.7468, 0.1652, 1,
    0.5252, 0.8381, 0.5089, 1,
    0.7095, 0.1897, 0.2258, 1,
    0.6979, 0.8537, 0.0461, 1,
    0.4186, 0.2026, 0.1419, 1,
    0.8318, 0.4289, 0.6408, 1,
    0.3475, 0.3354, 0.2868, 1,
    0.7025, 0.2123, 0.2339, 1,
    0.4944, 0.9446, 0.1233, 1,
    0.7193, 0.6551, 0.3655, 1,
    0.5491, 0.4076, 0.1014, 1 };
  vw::MatrixProxy<double,14,4> A(A_data);

  double xrot = 0.325;
  double yrot = 0.122;
  double zrot = 0.644;
  double dx = 5.6;
  double dy = 8.3;
  double dz = 2.1;

  Matrix4x4 S;

  submatrix(S, 0, 0, 3, 3) = euler_to_rotation_matrix(xrot, yrot, zrot, "xyz");
  select_col(S, 3) = Vector4(dx, dy, dz, 1);

  vw::Matrix<double> B = transpose(S*transpose(A));
  std::vector<Vector4> p1, p2;
  for (unsigned i = 0; i < A.rows(); ++i) {
    p1.push_back(select_row(A,i));
    p2.push_back(select_row(B,i));
  }

  EXPECT_MATRIX_NEAR( S, SimilarityFittingFunctorN<3>()(p1,p2), 1e-14 );
}

TEST(Geometry, SimilarityFittingFunctor) {
  static double A_data[] = {
    0.0153, 0.9318, 1,
    0.6721, 0.6813, 1,
    0.3046, 0.6822, 1,
    0.8600, 0.7468, 1,
    0.5252, 0.8381, 1,
    0.7095, 0.1897, 1,
    0.6979, 0.8537, 1,
    0.4186, 0.2026, 1,
    0.8318, 0.4289, 1,
    0.5417, 0.3784, 1 };

  vw::MatrixProxy<double,10,3> A(A_data);

  // Create a fictitious similarity matrix
  double s = 2.23;
  double theta = 0.2342;
  double dx = 4.2;
  double dy = 2.9;

  Matrix3x3 S;
  S(0,0) = s*cos(theta);
  S(0,1) = s*sin(theta);
  S(1,0) = s*-sin(theta);
  S(1,1) = s*cos(theta);
  S(0,2) = dx;
  S(1,2) = dy;
  S(2,2) = 1;

  // and apply it to some points
  vw::Matrix<double> B = transpose(S*transpose(A));
  std::vector<Vector3> p1, p2;
  for (unsigned i = 0; i < A.rows(); ++i) {
    p1.push_back(select_row(A,i));
    p2.push_back(select_row(B,i));
  }

  // compute a fit given the points only, and compare it to the true
  // similarity.
  EXPECT_MATRIX_NEAR( S, SimilarityFittingFunctor()(p1,p2), 1.8e-15 );
}

TEST(Geometry, TranslationRotationFittingFunctor) {
  static double A_data[] = {
    0.0153, 0.9318, 1,
    0.6721, 0.6813, 1,
    0.3046, 0.6822, 1,
    0.8600, 0.7468, 1,
    0.5252, 0.8381, 1,
    0.7095, 0.1897, 1,
    0.6979, 0.8537, 1,
    0.4186, 0.2026, 1,
    0.8318, 0.4289, 1,
    0.5417, 0.3784, 1 };

  vw::MatrixProxy<double,10,3> A(A_data);

  // Create a fictitious translation rotation matrix
  double theta = 0.2342;
  double dx = 4.2;
  double dy = 2.9;

  Matrix3x3 S;
  S(0,0) = cos(theta);
  S(0,1) = sin(theta);
  S(1,0) = -sin(theta);
  S(1,1) = cos(theta);
  S(0,2) = dx;
  S(1,2) = dy;
  S(2,2) = 1;

  // and apply it to some points
  vw::Matrix<double> B = transpose(S*transpose(A));
  std::vector<Vector3> p1, p2;
  for (unsigned i = 0; i < A.rows(); ++i) {
    p1.push_back(select_row(A,i));
    p2.push_back(select_row(B,i));
  }

  // compute a fit given the points only, and compare it to the true
  // similarity.
  EXPECT_MATRIX_NEAR( S, TranslationRotationFittingFunctor()(p1,p2), 1.8e-15 );
}

TEST(Geometry, TranslationFittingFunctor) {
  static double A_data[] = {
    0.0153, 0.9318, 1,
    0.6721, 0.6813, 1,
    0.3046, 0.6822, 1,
    0.8600, 0.7468, 1,
    0.5252, 0.8381, 1,
    0.7095, 0.1897, 1,
    0.6979, 0.8537, 1,
    0.4186, 0.2026, 1,
    0.8318, 0.4289, 1,
    0.5417, 0.3784, 1 };

  vw::MatrixProxy<double,10,3> A(A_data);

  // Create a translation matrix
  double dx = 4.2;
  double dy = 2.9;

  Matrix3x3 S;
  S.set_identity();
  S(0,2) = dx;
  S(1,2) = dy;

  // and apply it to some points
  vw::Matrix<double> B = transpose(S*transpose(A));
  std::vector<Vector3> p1, p2;
  for (unsigned i = 0; i < A.rows(); ++i) {
    p1.push_back(select_row(A,i));
    p2.push_back(select_row(B,i));
  }

  // compute a fit given the points only, and compare it to the true
  // similarity.
  EXPECT_MATRIX_NEAR( S, TranslationFittingFunctor()(p1,p2), 1.8e-15 );
}

TEST(Geometry, TranslationFittingFunctorN) {
  static double A_data[] = {
    0.0153, 0.9318, 0.2153, 1,
    0.6721, 0.6813, 0.4356, 1,
    0.3046, 0.6822, 0.7646, 1,
    0.8600, 0.7468, 0.5863, 1,
    0.5252, 0.8381, 0.2958, 1,
    0.7095, 0.1897, 0.9185, 1,
    0.6979, 0.8537, 0.3456, 1,
    0.4186, 0.2026, 0.1124, 1,
    0.8318, 0.4289, 0.7593, 1,
    0.5417, 0.3784, 0.4567, 1 };

  vw::MatrixProxy<double,10,4> A(A_data);

  // Create a translation matrix
  double dx = 4.2;
  double dy = 2.9;
  double dz = 3.1;

  Matrix4x4 S;
  S.set_identity();
  S(0,3) = dx;
  S(1,3) = dy;
  S(2,3) = dz;

  // and apply it to some points
  vw::Matrix<double> B = transpose(S*transpose(A));
  std::vector<Vector4> p1, p2;
  for (unsigned i = 0; i < A.rows(); ++i) {
    p1.push_back(select_row(A,i));
    p2.push_back(select_row(B,i));
  }

  // compute a fit given the points only, and compare it to the true
  // similarity.
  EXPECT_MATRIX_NEAR( S, TranslationFittingFunctorN<3>()(p1,p2), 1.8e-15 );
}


