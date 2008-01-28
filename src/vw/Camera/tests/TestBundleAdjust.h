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


// // Perform L*D*L^T decomposition on a symmetric semi-definite matrix
// // (in place) and then solves an equation of the form Ax=b.
// // 
// // WARNING: Modifies the contents of the matrix A.
// template <class MatrixT, class VectorT>
// vw::Vector<double> sparse_solve(MatrixT& A, VectorT const& b) {
//   VW_ASSERT(A.cols() == A.rows(), ArgumentErr() << "sparse_solve: matrix must be square and symmetric.\n");

//   std::cout << "Running unoptimized LDL^T decomposition.\n";
//   ldl_decomposition(A);

// //   vw::Matrix<double> Ltest(A.cols(), A.rows());
// //   vw::Matrix<double> Dtest(A.cols(), A.rows());
// //   for (unsigned i = 0; i < A.rows(); ++i) {
// //     for (unsigned j = 0; j < i; ++j) {
// //       Ltest(i,j) = A(i,j);
// //     }
// //   }

// //   for (unsigned i = 0; i < A.rows(); ++i) {
// //     Ltest(i,i) = 1;
// //     Dtest(i,i) = A(i,i);
// //   }

// //   std::cout << "A MATRIX:\n";
// //   print_matrix(A);
// //   std::cout << "L MATRIX:\n";
// //   print_matrix(Ltest);
// //   std::cout << "D MATRIX:\n";
// //   print_matrix(Dtest);


// //   std::cout << "Forward substituting.\n";
//   vw::Vector<double> x_prime(A.cols());
//   for (unsigned i = 0; i < x_prime.size(); ++i) {
//     double sum = 0;
//     for (unsigned j = 0; j < i; ++j) 
//       sum += A(i,j)*x_prime(j);
//     x_prime(i) = b(i)-sum;
//   }

// //   std::cout << "Dividing by D..\n";
//   vw::Vector<double> x_doubleprime(A.cols());
//   for (unsigned i = 0; i < x_doubleprime.size(); ++i) {
//     x_doubleprime(i) = x_prime(i)/A(i,i);
//   }

// //   vw::Vector<double> x_doubleprimetest = inverse(Dtest)*x_prime;
// //   std::cout << x_doubleprime << "\n";
// //   std::cout << x_doubleprimetest << "\n";
 
//   vw::Vector<double> x(A.cols());
//   //  std::cout << "Back substituting.\n";
//   for (int i = x.size()-1; i >= 0; --i) {
//     double sum = 0;
//     for (unsigned j = i; j < A.cols(); ++j) 
//       sum += A(j,i)*x(j);
//     x(i) = x_doubleprime(i) - sum;
//   }

// //   vw::Vector<double> x_test = inverse(transpose(Ltest))*x_doubleprime;
// //   std::cout << x << "\n";
// //   std::cout << x_test << "\n";
  
//   return x;
// }

// // Perform L*D*L^T decomposition on a sparse skyline symmetric
// // semi-definite matrix (in place) and then solves an equation of the
// // form Ax=b.
// // 
// // WARNING: Modifies the contents of the matrix A.
// template <class ElemT, class VectorT>
// vw::Vector<double> sparse_solve(SparseSkylineMatrix<ElemT>& A, VectorT const& b) {
//   VW_ASSERT(A.cols() == A.rows(), ArgumentErr() << "sparse_solve: matrix must be square and symmetric.\n");

//   std::cout << "Running skyline optimized LDL^T decomposition.\n";
//   ldl_decomposition(A);

//   //  vw::Matrix<double> Ltest(A.cols(), A.rows());
//   //  vw::Matrix<double> Dtest(A.cols(), A.rows());
// //   print_matrix(A);
// //   A.print_sparse_structure();
// //   for (unsigned i = 0; i < A.rows(); ++i) {
// //     for (unsigned j = 0; j < i; ++j) {
// //       Ltest(i,j) = A(i,j);
// //     }
// //   }
// //   print_matrix(A);
// //   A.print_sparse_structure();

// //   for (unsigned i = 0; i < A.rows(); ++i) {
// //     Ltest(i,i) = 1;
// //     Dtest(i,i) = A(i,i);
// //   }

// //   A.print_sparse_structure();

// //   std::cout << "A MATRIX:\n";
// //   print_matrix(A);
// //   std::cout << "L MATRIX:\n";
// //   print_matrix(Ltest);
// //   std::cout << "D MATRIX:\n";
// //   print_matrix(Dtest);

//   const std::vector<uint32>& skyline = A.skyline();
//   std::vector<uint32> inverse_skyline(skyline.size());

//   // Construct the inverse skyline matrix, which is used to optimize the final
//   // back substitution step below.
//   for (unsigned j = 0; j < inverse_skyline.size(); ++j) {
//     inverse_skyline[j] = 0;
//     for (int i = skyline.size()-1; i>=0; --i) {
//       if (j < skyline[i]) 
//         ++(inverse_skyline[j]);
//       else
//         break; // Break out of the inner loop
//     }
//   }
    

// //   std::cout << "Forward substituting.\n";
//   vw::Vector<double> x_prime(A.cols());
//   for (unsigned i = 0; i < x_prime.size(); ++i) {
//     double sum = 0;
//     for (unsigned j = skyline[i]; j < i; ++j) 
//       sum += A(i,j)*x_prime(j);
//     x_prime(i) = b(i)-sum;
//   }

// //   std::cout << "Dividing by D..\n";
//   vw::Vector<double> x_doubleprime(A.cols());
//   for (unsigned i = 0; i < x_doubleprime.size(); ++i) {
//     x_doubleprime(i) = x_prime(i)/A(i,i);
//   }

// //   vw::Vector<double> x_doubleprimetest = inverse(Dtest)*x_prime;
// //   std::cout << x_doubleprime << "\n";
// //   std::cout << x_doubleprimetest << "\n";
 
//   vw::Vector<double> x(A.cols());
//   //   std::cout << "Back substituting.\n";
//   for (int32 i = x.size()-1; i >= 0; --i) {
//     double sum = 0;
//     for (unsigned j = i+1; j < A.cols()-inverse_skyline[i]; ++j) 
//       sum += A(j,i)*x(j);
//     x(i) = x_doubleprime(i) - sum;
//   }

// //   vw::Vector<double> x_test = inverse(transpose(Ltest))*x_doubleprime;
// //   std::cout << x << "\n";
// //   std::cout << x_test << "\n";
  
//   return x;
// }



  struct BundleAdjustBase {
   
    // Given the 'a' vector (camera model parameters) for the j'th
    // image, and the 'b' vector (3D point location) for the i'th
    // point, return the location of b_i on imager j in pixel
    // coordinates.
    Vector2 operator() ( Vector<double> const& a_j, Vector<double> const& b_i ) const {

      // Apply the "camera model" and return the imaged point
      return Vector2(0,0);

    }

    Vector2 error( Vector<double> const& x_ij, Vector<double> const& a_j, Vector<double> const& b_i ) const {
      return x_ij - this->operator()(a_j,b_i);
    }

    // Approximate the jacobian for small variations in the a_j
    // parameters (camera parameters). 
    Matrix<double> A_jacobian ( Vector<double> const& a_j, Vector<double> const& b_i ) const {
      // Get nominal function value
      Vector2 h0 = this->operator()(a_j, b_i);

      // Jacobian is #outputs x #params
      Matrix<double> J(h0.size(), a_j.size());

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( unsigned i=0; i<a_j.size(); ++i ){
        Vector<double> a_j_prime = a_j;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(a_j(i)*1e-7);
        a_j_prime(i) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        Vector<double> hi = this->operator()(a_j_prime, b_i);
        select_col(J,i) = (hi-h0)/epsilon;
      }
      return J;
    }

    // Approximate the jacobian for small variations in the b_i
    // parameters (3d point locations). 
    Matrix<double> B_jacobian ( Vector<double> const& a_j, Vector<double> const& b_i ) const {
      // Get nominal function value
      Vector2 h0 = this->operator()(a_j, b_i);

      // Jacobian is #outputs x #params
      Matrix<double> J(h0.size(), b_i.size());

      // For each param dimension, add epsilon and re-evaluate h() to
      // get numerical derivative w.r.t. that parameter
      for ( unsigned i=0; i<b_i.size(); ++i ){
        Vector<double> b_i_prime = b_i;

        // Variable step size, depending on parameter value
        double epsilon = 1e-7 + fabs(b_i(i)*1e-7);
        b_i_prime(i) += epsilon;

        // Evaluate function with this step and compute the derivative
        // w.r.t. parameter i
        Vector<double> hi = this->operator()(b_i_prime, b_i);
        select_col(J,i) = (hi-h0)/epsilon;
      }
      return J;
    }

    // Each entry in the outer vector corresponds to a distinct 3D
    // point.  The inner vector contains a list of image IDs and
    // pixel coordinates where that point was imaged.
    typedef std::vector<std::pair<int, Vector2> > bundle_type;
    typedef std::vector< bundle_type > bundles_type;
    
    void update(std::vector<Vector<double> > const& a, std::vector<Vector<double> > const& b, bundles_type const& points, int num_images, double lambda) const {

      mapped_matrix<Matrix<double> > A(points.size(), num_images);
      mapped_matrix<Matrix<double> > B(points.size(), num_images);
      mapped_matrix<Vector2> epsilon(points.size(), num_images);

      mapped_vector< Matrix<double> > U(A.size1());  // size1() == rows
      mapped_vector< Matrix<double> > V(B.size1());  
      mapped_matrix<Matrix<double> > W(points.size(), num_images);
      mapped_vector< Vector2 > epsilon_a(A.size1());
      mapped_vector< Vector2 > epsilon_b(B.size1());
      mapped_matrix<Matrix<double> > Y(points.size(), num_images);
      //      mapped_matrix<Matrix<double> > S(points.size(), num_image);

    
      // Populate the Jacobian, which is broken into two sparse
      // matrices A & B, as well as the error matrix and the W 
      // matrix.
      int i = 0;
      for (bundles_type::const_iterator iter = points.begin(); iter != points.end(); ++iter) {
        for (bundle_type::const_iterator point_iter = (*iter).begin(); point_iter != (*iter).end(); ++point_iter) {
          int j = (*point_iter).first;

          // Store jacobian values
          A(i,j) = A_jacobian(a[j],b[i]);
          B(i,j) = B_jacobian(a[j],b[i]);

          // Compute error vector
          epsilon(i,j) = this->error((*point_iter).second,a[j],b[i]);

          // Store intermediate values
          U(j) += transpose(A(i,j).ref()) * /* inverse(SIGMA)* */ A(i,j).ref();
          V(i) += transpose(B(i,j).ref()) * /* inverse(SIGMA)* */ B(i,j).ref();
          W(i,j) = transpose(A(i,j).ref()) * /* inverse(SIGMA)* */ B(i,j).ref();
          epsilon_a(j) += transpose(A(i,j).ref()) * /* inverse(SIGMA)* */ epsilon(i,j).ref();
          epsilon_b(i) += transpose(B(i,j).ref()) * /* inverse(SIGMA)* */ epsilon(i,j).ref();

          /// MUST ONLY AUGMENT DIAGONAL ENTRIES!!!
          //          Y(i,j) = W(i,j).ref() * inverse( V(i).ref() * (1.0f+lambda) ); 

          ++i;
        }
      }

      // Second pass
      //
      //  NEED TO INCORPORATE THE CORRECT USE OF K!!!
      //

//       i = 0;
//       for (bundles_type::iterator iter = points.begin(); iter != points.end(); ++iter) {
//         for (bundle_type::iterator point_iter = *iter.begin(); point_iter != *iter.end(); ++point_iter) {
//           int j = (*point_iter).first;
//           if (i==j) 
          /// MUST ONLY AUGMENT DIAGONAL ENTRIES!!!
//             S(i,j) -= ( Y(i,j).ref()*transpose(W(i,j).ref()) + (U(j).ref() * (1.0f+lambda)) );
//           else  
//             S(i,j) -= ( Y(i,j).ref()*transpose(W(i,j).ref()) );
          
//           ea(j) += epsilon_a(j).ref() - (Y(i,j).ref() * epsilon_b(i).ref());
//           ++i;
//         }
//       }
      
      //       mapped_vector<double> delta_a = inverse(S) * ea;

      // Third pass
//       i = 0;
//       for (bundles_type::iterator iter = points.begin(); iter != points.end(); ++iter) {
//         for (bundle_type::iterator point_iter = *iter.begin(); point_iter != *iter.end(); ++point_iter) {
//           eb(j) += epsilon_b(i).ref() - (transpose(W(i,j).ref()) * delta_a(j).ref());
//           ++i;
//         }
//       }

      /// MUST ONLY AUGMENT DIAGONAL ENTRIES!!!
      //      mapped_vector<double> delta_b = inverse( V[i] * (1.0f+lambda) ) * eb;
    }

  };

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
