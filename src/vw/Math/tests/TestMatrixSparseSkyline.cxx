// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <gtest/gtest.h>
#include <vw/Math/MatrixSparseSkyline.h>

using namespace vw;
using namespace vw::math;

static const double DELTA = 1e-5;

Vector<uint> create_test_skyline(unsigned size, int max_offset ) {
  Vector<uint> result(size);
  for ( uint i = 0; i < size; ++i ) {
    int offset = int(i) - rand()%max_offset;
    if ( offset < 0 ) offset = 0;
    result[i] = offset;
  }
  return result;
}

template <class VectorT>
void fill_vector(VectorT& b) {
  for ( unsigned i = 0; i < b.size(); ++i )
    b[i] = lround(double(random())/(pow(2,31)-1)*100)+1;
}

template <class MatrixT>
void fill_symmetric_matrix(MatrixT& A, Vector<uint> const& skyline ) {
  for (unsigned i = 0; i < A.rows(); ++i)
    for (unsigned j = skyline[i]; j < std::min(i+1,A.cols()); ++j) {
      A(i,j) = lround(double(random())/(pow(2,31)-1)*100)+1;
      A(j,i) = A(i,j);
    }
}

template <class SrcMatrixT, class DestMatrixT>
void copy_symmetric_matrix(SrcMatrixT const& src, DestMatrixT& dest, Vector<uint> const& skyline ) {
  for (unsigned i = 0; i < src.rows(); ++i )
    for (unsigned j = skyline[i]; j < std::min(i+1,src.cols()); ++j) {
      dest(i,j) = src(i,j);
      dest(j,i) = dest(i,j);
    }
}

// Unoptimized LDLT decomposition
template <class MatrixT>
void ldl_decomposition(MatrixT& A) {
  VW_ASSERT(A.cols() == A.rows(), ArgumentErr() << "ldl_decomposition: argument must be square and symmetric.\n");
  for (unsigned j = 0; j < A.cols(); ++j) {

    // Compute v(1:j)
    std::vector<double> v(j+1);
    v[j] = A(j,j);
    for (unsigned i = 0; i < j; ++i) {
      v[i] = A(j,i)*A(i,i);
      v[j] -= A(j,i)*v[i];
    }

    // Store d(j) and compute L(j+1:n,j)
    A(j,j) = v[j];
    for (unsigned i = j+1; i < A.cols(); ++i) {
      double row_sum = 0;
      for (unsigned jj = 0; jj < j; ++jj)
        row_sum += A(i,jj)*v[jj];
      A(i,j) = ( A(i,j)-row_sum ) / v[j];
    }
  }
}

TEST(SparseSkyline, LDL_decomp_correctness) {
  uint N = 50;
  uint S = 10;

  srandom((unsigned int) clock());

  Matrix<double> nonsparse_mat(N,N);
  MatrixSparseSkyline<double> sparse_mat(N,N);
  MatrixSparseSkyline<double> original_sparse_mat(N,N);

  Vector<uint> test_skyline = create_test_skyline(N, S);
  //std::cout << "TestSkyline:" << test_skyline << "\n";

  fill_symmetric_matrix(sparse_mat, test_skyline);
  copy_symmetric_matrix(sparse_mat, nonsparse_mat, test_skyline);
  original_sparse_mat = sparse_mat;

  //std::cout << "SparseStructure: " << sparse_mat << "\n";

  sparse_ldl_decomposition(sparse_mat);
  for ( uint i = 0; i < N; i++ )
    for ( uint j = 0; j < N; j++ )
      if ( (sparse_mat(i,j) == 0)^(original_sparse_mat(i,j) == 0) )
        FAIL() << "Sparse structure was not preserved by sparse LDLT decomposition.\nIndex(" << i << "," << j << ") is " << sparse_mat(i,j) << " and " << original_sparse_mat(i,j) << "\n";

  ldl_decomposition(nonsparse_mat);

  for ( uint i = 0; i < N; i++ )
    for ( uint j = 0; j < i; j++ )
      EXPECT_NEAR(sparse_mat(i,j),nonsparse_mat(i,j),DELTA);
}

TEST(SparseSkyline, LDL_decomp_scalability) {
  uint N = 5000;
  uint S = 150;

  srandom((unsigned int) clock());

  MatrixSparseSkyline<double> sparse_mat(N,N);
  MatrixSparseSkyline<double> original_sparse_mat(N,N);
  Vector<uint> test_skyline = create_test_skyline(N,S);

  fill_symmetric_matrix(sparse_mat, test_skyline);
  original_sparse_mat = sparse_mat;

  sparse_ldl_decomposition(sparse_mat);
  for ( uint i = 0; i < N; i++ )
    for ( uint j = 0; j < i; j++ )
      if ( (sparse_mat(i,j) == 0)^(original_sparse_mat(i,j) == 0) )
        FAIL() << "Sparse structure was not preserved by sparse LDLT decomposition.\nIndex(" << i << "," << j << ") is " << sparse_mat(i,j) << " and " << original_sparse_mat(i,j) << "\n";
}

TEST(SparseSkyline, LDL_solve) {
  uint N = 100;
  uint S = 10;

  srandom((unsigned int) clock());

  Matrix<double> A_nonsparse(N,N);
  MatrixSparseSkyline<double> A_sparse(N,N);
  Vector<uint> test_skyline = create_test_skyline(N,S);

  fill_symmetric_matrix(A_sparse, test_skyline);
  copy_symmetric_matrix(A_sparse, A_nonsparse, test_skyline);

  // Create a vector to solve against
  Vector<double> b(N);
  fill_vector(b);

  // Solve using normal
  Vector<double> x_nonsparse = inverse(A_nonsparse)*b;

  // Sparse Version
  Vector<double> x_sparse = sparse_solve(A_sparse, b);

  for ( uint i = 0; i < N; i++ )
    EXPECT_NEAR( x_nonsparse[i], x_sparse[i], DELTA );
}

TEST(SparseSkyline, LDL_solve_scalability) {
  int N = 1000;
  int S = 100;

  srandom((unsigned int) clock());

  MatrixSparseSkyline<double> A_sparse(N,N);
  Vector<uint> test_skyline = create_test_skyline(N,S);
  fill_symmetric_matrix(A_sparse, test_skyline);

  // Create a vector to solve against
  Vector<double> b(A_sparse.cols());
  fill_vector(b);

  // Solving for X
  Vector<double> x_result = sparse_solve(A_sparse, b);

  // Double checking that x_result is correct
  //  Vector<double> b_prime = A_sparse*x_result;

  //for ( uint i = 0; i < b_prime.size(); i++ )
  //  EXPECT_NEAR( b[i], b_prime[i], DELTA );
}
