// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

#include <cxxtest/TestSuite.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>
#include <vw/Camera/BundleAdjust.h>

#include <stdlib.h>
#include <time.h> 

using namespace vw;
using namespace vw::math;
using namespace vw::camera;

using namespace boost::numeric::ublas;

std::vector<uint32> create_test_skyline(unsigned size, int max_offset) {
  std::vector<uint32> result(size);

  for (int i = 0; i < int(size); ++i) {
    int offset = lround(double(random())/((pow(2,31))-1)*max_offset);
    if (i - offset < 0) offset = i;
    result[i] = i - offset;
  }
  return result;
}

template <class VectorT>
void fill_vector(VectorT& b) {
  for (unsigned i = 0; i < b.size(); ++i) {
    b[i] = lround(double(random())/((pow(2,31))-1)*100)+1;
  }
}

template <class MatrixT>
void fill_symmetric_matrix(MatrixT& A, std::vector<uint32> const& skyline) {
  for (unsigned i = 0; i < A.rows(); ++i) {
    for (unsigned j = skyline[i]; j < std::min(i+1,A.cols()); ++j) {
      //      A.push_back(i,j,lround(double(random())/((pow(2,31))-1)*100)+1);
      A(i,j) =  lround(double(random())/((pow(2,31))-1)*100)+1;
      A(j,i) = A(i,j);
    }
  }
}

template <class SrcMatrixT, class DestMatrixT>
void copy_symmetric_matrix(SrcMatrixT const& src, DestMatrixT& dest,std::vector<uint32> const& skyline) {
  for (unsigned i = 0; i < src.rows(); ++i) {
    for (unsigned j = skyline[i]; j < std::min(i+1, src.cols()); ++j) {
      dest(i,j) = src(i,j);
      dest(j,i) = dest(i,j);
    }
  }
}

template <class MatrixT>
void print_matrix(MatrixT const& A) {
  for (int i = 0; i < A.rows(); ++i) {
    for (int j = 0; j < A.cols(); ++j)
      std::cout << A(i,j) << " ";
    std::cout << "\n";
  }
}



// Note: this is a non-skyline optimized factorization that is used to
// test correctness in the unit tests below.  
//
// Perform L*D*L^T decomposition on a symmetric semi-definite
// matrix.  WARNING: The results are stored in place, so this
// operation destroys the previous contents of A.
//
// Once this operation is complete, the diagonal entries of A
// contain the values from D, and the lower left block diagonal of A
// contains L.  (The diagonal entries of L are always 1, so those
// are assumed here...)
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

class TestOptimization : public CxxTest::TestSuite
{
public:

  void test_LDL_decomp_correctness()
  {
    int N = 10;
    int S = 2;

    srandom((unsigned int) clock());

    Matrix<double> nonsparse_mat(N,N);
    SparseSkylineMatrix<double> sparse_mat(N,N);
    SparseSkylineMatrix<double> original_sparse_mat(N,N);

    std::vector<uint32> test_skyline = create_test_skyline(N, S);
//     std::cout << "TEST SKYLINE1\n ";
//     for (int i = 0; i < test_skyline.size(); ++i)
//       std::cout << test_skyline[i] << " ";
//     std::cout << "\n";

    fill_symmetric_matrix(sparse_mat, test_skyline);
    copy_symmetric_matrix(sparse_mat, nonsparse_mat, test_skyline);
    copy_symmetric_matrix(sparse_mat, original_sparse_mat, test_skyline);

//     std::vector<uint32> sline = sparse_mat.skyline();
//     std::cout << "SPARSE_MAT SKYLINE1\n ";
//     for (int i = 0; i < sline.size(); ++i)
//       std::cout << sline[i] << " ";
//     std::cout << "\n";

//     std::cout << "Original Sparse_Mat:\n";
//     print_matrix(sparse_mat);
//     sparse_mat.print_sparse_structure();

//     std::cout << "Original Nonsparse_Mat:\n";
//     print_matrix(nonsparse_mat);


//    std::cout << "Running LDL^T decomposition on sparse skyline matrix..." << std::flush;
    ldl_decomposition(sparse_mat);
    //    std::cout << " done.\n";
    // Check to make sure the decomposod LDL matrix still has the same
    // sparse structure as the original.
    for (unsigned i=0; i < sparse_mat.rows(); ++i) {
      for (unsigned j=0; j < i; ++j) {
        if (sparse_mat.find_element(i,j) && !original_sparse_mat.find_element(i,j))
          TS_FAIL("Sparse structure was not preserved by the LDL^T decomposition!\n");
      }
    }

    //    std::cout << "Running LDL^T decomposition on normal VW matrix..." << std::flush;
    ldl_decomposition(nonsparse_mat);    
    //    std::cout << " done.\n";

//     std::cout << "New Sparse_Mat:\n";
//     print_matrix(sparse_mat);
//     sparse_mat.print_sparse_structure();

//     std::cout << "New Nonsparse_Mat:\n";
//     print_matrix(nonsparse_mat);


    // Check to make sure the skyline-optimized LDL decomp produces
    // the same results as the non-optimized implementation.
    for (unsigned i=0; i < sparse_mat.rows(); ++i) {
      for (unsigned j=0; j < i; ++j) {
        if (fabs(sparse_mat(i,j).ref() - nonsparse_mat(i,j)) > 0.00001)
          std::cout << "Mismatch: " << i << " " << j << "\n";
        TS_ASSERT_DELTA(sparse_mat(i,j).ref(), nonsparse_mat(i,j), 0.0001);
      }
    }

  }

  void test_LDL_decomp_scalability()
  {
    int N = 1000;
    int S = 10;

    srandom((unsigned int) clock());

    SparseSkylineMatrix<double> sparse_mat(N,N);
    SparseSkylineMatrix<double> original_sparse_mat(N,N);

    std::vector<uint32> test_skyline = create_test_skyline(N, S);

    std::cout << "\nBuilding " << N << " x " << N << " sparse matrix.\n";
    fill_symmetric_matrix(sparse_mat, test_skyline);
    copy_symmetric_matrix(sparse_mat, original_sparse_mat, test_skyline);

    std::cout << "Running LDL^T decomposition on sparse skyline matrix..." << std::flush;
    ldl_decomposition(sparse_mat);
    std::cout << " done.\n";
    // Check to make sure the decomposod LDL matrix still has the same
    // sparse structure as the original.
    for (unsigned i=0; i < sparse_mat.rows(); ++i) {
      for (unsigned j=0; j < i; ++j) {
        if (sparse_mat.find_element(i,j) && !original_sparse_mat.find_element(i,j))
          TS_FAIL("Sparse structure was not preserved by the LDL^T decomposition!\n");
      }
    }
  }

  void test_LDL_solve()
  {
    int N = 100;
    int S = 2;

    srandom((unsigned int) clock());

    Matrix<double> A_nonsparse(N,N);
    SparseSkylineMatrix<double> A_sparse(N,N);

    std::vector<uint32> test_skyline = create_test_skyline(N, S);

    fill_symmetric_matrix(A_sparse, test_skyline);
    copy_symmetric_matrix(A_sparse, A_nonsparse, test_skyline);
    
    // Create a 'b' vector with random entries
    vw::Vector<double> b(A_sparse.cols());
    fill_vector(b);

    // First, solve using normal matrix inverse.
    //    std::cout << "Solving non-sparse " << N << "x" << N << " system... " << std::flush;
    Vector<double> x_nonsparse = inverse(A_nonsparse) * b;
    //    std::cout << " done.\n";

    // Next, try our fast sparse skyline solver.
    //    std::cout << "Solving sparse " << N << "x" << N << " system... " << std::flush;
    Vector<double> x_sparse = sparse_solve(A_sparse, b);
    //    std::cout << " done.\n";

    // Check to make sure the skyline-optimized LDL decomp produces
    // the same results as the non-optimized implementation.
    for (unsigned i=0; i < x_sparse.size(); ++i) {
      TS_ASSERT_DELTA(x_sparse[i], x_nonsparse[i], 0.000001);
    }

  }

  void test_LDL_solve_scalability()
  {
    int N =1000;
    int S = 30;

    srandom((unsigned int) clock());

    SparseSkylineMatrix<double> A_sparse(N,N);

    std::vector<uint32> test_skyline = create_test_skyline(N, S);

    fill_symmetric_matrix(A_sparse, test_skyline);
    
    // Create a 'b' vector with random entries
    vw::Vector<double> b(A_sparse.cols());
    fill_vector(b);

    std::cout << "\nSolving sparse " << N << "x" << N << " system... " << std::flush;
    Vector<double> x_sparse = sparse_solve(A_sparse, b);
  }

}; // class TestOptimization
