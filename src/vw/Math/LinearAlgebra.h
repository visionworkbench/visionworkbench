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

/// \file LinearAlgebra.h
/// 
/// Provides standard linear algebra functions, such as eigendecomposition, 
/// least-squares, and so on, utilizing LAPACK where possible.
/// Currently provides:
///  - Eigendecomposition via eigen()
///  - Singular value decomposition via svd() and complete_svd()
///  - Pseudoinverse via pseudoinverse()
#ifndef __VW_MATH_LINEAR_ALGEBRA_H__
#define __VW_MATH_LINEAR_ALGEBRA_H__

#include <vector>

#include <vw/config.h>
#include <vw/Core/Exception.h>
#include <vw/Math/Matrix.h>

namespace vw {
/// Numerical linear algebra and computational geometry.
namespace math {

  /// \cond INTERNAL
  extern "C"  int sgeev_(char *jobvl, char *jobvr, long int *n, float *a, 
                         long int *lda, float *wr, float *wi, float *vl, long int *ldvl, float *vr, 
                         long int *ldvr, float *work, long int *lwork, long int *info);

  extern "C"  int dgeev_(char *jobvl, char *jobvr, long int *n, double *a, 
                         long int *lda, double *wr, double *wi, double *vl, long int *ldvl, double *vr, 
                         long int *ldvr, double *work, long int *lwork, long int *info);


  static inline void geev(char jobvl, char jobvr, long int n, float *a, long int lda, float *wr, float *wi, float *vl, long int ldvl, float *vr, long int ldvr, float *work, long int lwork, long int *info) {
    sgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, info);
  }

  static inline void geev(char jobvl, char jobvr, long int n, double *a, long int lda, double *wr, double *wi, double *vl, long int ldvl, double *vr, long int ldvr, double *work, long int lwork, long int *info) {
    dgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, info);
  }
  /// \endcond

  /// Compute the eigenvalues of a Matrix A
  ///
  /// E must be a complex vector and will be resized if necessary to
  /// the appropriate output dimensions based on the dimensions of A.
  template <class MatrixT, class EigenvaluesT> 
  inline void eigen( MatrixT const& A, EigenvaluesT &e ) {
    VW_ASSERT( A.cols()==A.rows(), ArgumentErr() << "Eigendecomposition can only be performed on square matrices." );
    typedef typename MatrixT::value_type real_type;
    const long int lda = A.rows();
    Matrix<real_type> Abuf = transpose(A);
    Vector<real_type> wr_buf( A.cols() );
    Vector<real_type> wi_buf( A.cols() );
    real_type work_size;
    long int n = A.cols(); 
    long int lwork = -1, info;
    geev('N','N',n,&(Abuf(0,0)), lda, &(wr_buf(0)), &(wi_buf(0)), NULL, 1, NULL, 1, &work_size, lwork, &info);
    lwork = (long int)(work_size);
    std::vector<real_type> work( lwork );
    geev('N','N',n,&(Abuf(0,0)), lda, &(wr_buf(0)), &(wi_buf(0)), NULL, 1, NULL, 1, &work[0], lwork, &info);
    if (info < 0) 
      vw_throw( ArgumentErr() << "eigen(): LAPACK driver geev reported an error with argument " << -info << "." );
    if (info > 0) 
      vw_throw( ArgumentErr() << "eigen(): LAPACK driver geev only converged for the first " << info << " eigenvectors." );
    e.set_size( A.cols() );
    for ( unsigned i = 0; i < wr_buf.size(); ++i ) 
      e(i) = std::complex<real_type>(wr_buf(i), wi_buf(i));
  }

  /// Compute the entire eigendecomposition of the matrix A.
  ///
  /// E and V will contain the resulting eigendecomposition.  E will
  /// contain the eigenvalues and must be a complex vector, while V
  /// will contain eigenvectors and must be a complex matrix. E and V
  /// will be resized based on the dimensions of A.
  template <class AMatrixT, class EigenvaluesT, class VMatrixT> 
  inline void eigen( AMatrixT &A, EigenvaluesT &e, VMatrixT &V ) {
    VW_ASSERT( A.cols()==A.rows(), ArgumentErr() << "Eigendecomposition can only be performed on square matrices." );
    typedef typename AMatrixT::value_type real_type;
    const long int lda = A.rows();
    const long int ldvr = A.cols();
    Matrix<real_type> Abuf = transpose(A);
    Matrix<real_type> Vbuf( A.rows(), A.cols() );
    Vector<real_type> wr_buf( A.cols() );
    Vector<real_type> wi_buf( A.cols() );
    real_type work_size;
    long int n = A.cols(), lwork = -1, info;
    geev('N','V',n,&(Abuf(0,0)), lda, &(wr_buf(0)), &(wi_buf(0)), NULL, 1, &(Vbuf(0,0)), ldvr, &work_size, lwork, &info);
    lwork = (long int)(work_size);
    std::vector<real_type> work( lwork );
    geev('N','V',n,&(Abuf(0,0)), lda, &(wr_buf(0)), &(wi_buf(0)), NULL, 1, &(Vbuf(0,0)), ldvr, &work[0], lwork, &info);
    if (info < 0) 
      vw_throw( ArgumentErr() << "eigen(): LAPACK driver geev reported an error with argument " << -info << "." );
    if (info > 0) 
      vw_throw( ArgumentErr() << "eigen(): LAPACK driver geev only converged for the first " << info << " eigenvectors." );
    e.set_size( A.cols() );
    V.set_size( Vbuf.cols(), Vbuf.rows() );
    for ( unsigned i = 0; i < wr_buf.size(); ++i ) {
      e(i) = std::complex<real_type>(wr_buf(i), wi_buf(i));
      // If the eigenvalue is complex, we must tease the real and
      // complex parts out of the Vbuf matrix.
      for ( unsigned r = 0; r < V.rows(); ++r )
        if (wi_buf(i) == 0)
          V(r,i) = Vbuf(i,r);
        else if (wi_buf(i) > 0)
          V(r,i) = std::complex<real_type>(Vbuf(i,r), Vbuf(i+1,r));
        else
          V(r,i) = std::complex<real_type>(Vbuf(i-1,r), -Vbuf(i,r));
    }
  }

  /// \cond INTERNAL
  extern "C"  int sgesdd_(char *jobz, long int *m, long int *n, float *a, 
                          long int *lda, float *s, float *u, long int *ldu, float *vt, long int *ldvt,
                          float *work, long int *lwork, long int *iwork, long int *info);

  extern "C"  int dgesdd_(char *jobz, long int *m, long int *n, double *a,
                          long int *lda, double *s, double *u, long int *ldu, double *vt, long int *ldvt, 
                          double *work, long int *lwork, long int *iwork, long int *info);

  static inline void gesdd(char jobz, long int m, long int n, float *a, long int lda, float *s, float *u, long int ldu, float *vt, long int ldvt, float *work, long int lwork, long int *iwork, long int *info) { 
    sgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, info);
  }

  static inline void gesdd(char jobz, long int m, long int n, double *a, long int lda, double *s, double *u, long int ldu, double *vt, long int ldvt, double *work, long int lwork, long int *iwork, long int *info) { 
    dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, info);
  }
  /// \endcond

  /// Compute the singular values of a Matrix A
  ///
  /// S will be resized if necessary to the appropriate output 
  /// dimensions based on the dimensions of A.
  template <class AMatrixT, class SingularValuesT> 
  inline void svd( AMatrixT const& A, SingularValuesT &s ) {
    typedef typename PromoteType<typename AMatrixT::value_type, typename SingularValuesT::value_type>::type real_type;
    const long int m = A.rows(), n = A.cols();
    const long int minmn = std::min(m,n);
    const long int lda = A.rows();
    Matrix<real_type> Abuf = transpose(A);
    Vector<real_type> sbuf( minmn );
    real_type work_size;
    long int lwork = -1, info;
    std::vector<long int> iwork( 8*minmn );
    gesdd('N', m, n, &(Abuf(0,0)), lda, &(sbuf(0)), NULL, 1, NULL, 1, &work_size, lwork, &iwork[0], &info);
    lwork = (long int)(work_size);
    std::vector<real_type> work( lwork );
    gesdd('N', m, n, &(Abuf(0,0)), lda, &(sbuf(0)), NULL, 1, NULL, 1, &work[0], lwork, &iwork[0], &info);
    if (info < 0) 
      vw_throw( ArgumentErr() << "svd(): LAPACK driver gesdd reported an error with argument " << -info << "." );
    if (info > 0) 
      vw_throw( ArgumentErr() << "svd(): LAPACK driver gesdd did not converge.  Update process failed." );
    s = sbuf;
  }
  
  /// Compute the singular value decomposition of the matrix A.
  ///
  /// U, S, and VT will be resized if necessary to the appropriate output 
  /// dimensions based on the dimensions of A.
  template <class AMatrixT, class UMatrixT, class SingularValuesT, class VTMatrixT> 
  inline void svd( AMatrixT const& A, UMatrixT &U, SingularValuesT &s, VTMatrixT &VT ) {
    typedef typename PromoteType<typename AMatrixT::value_type, typename SingularValuesT::value_type>::type temp_type1;
    typedef typename PromoteType<temp_type1, typename UMatrixT::value_type>::type temp_type2;
    typedef typename PromoteType<temp_type2, typename VTMatrixT::value_type>::type real_type;
    const long int m = A.rows(), n = A.cols();
    const long int minmn = std::min(m,n);
    const long int lda = A.rows();
    Matrix<real_type> Abuf = transpose(A);
    Matrix<real_type> Ubuf( minmn, A.rows() );
    Matrix<real_type> VTbuf( A.cols(), minmn );
    Vector<real_type> sbuf( minmn );
    real_type work_size;
    long int lwork = -1, info;
    std::vector<long int> iwork( 8*minmn );
    long int ldu = m, ldvt = minmn;
    gesdd('S', m, n, &(Abuf(0,0)), lda, &(sbuf(0)), &(Ubuf(0,0)), ldu, &(VTbuf(0,0)), ldvt, &work_size, lwork, &iwork[0], &info);
    lwork = (long int)(work_size);
    std::vector<real_type> work( lwork );
    gesdd('S', m, n, &(Abuf(0,0)), lda, &(sbuf(0)), &(Ubuf(0,0)), ldu, &(VTbuf(0,0)), ldvt, &work[0], lwork, &iwork[0], &info);
    if (info < 0) 
      vw_throw( ArgumentErr() << "svd(): LAPACK driver gesdd reported an error with argument " << -info << "." );
    if (info > 0) 
      vw_throw( ArgumentErr() << "svd(): LAPACK driver gesdd did not converge.  Update process failed." );
    U = transpose(Ubuf);
    VT = transpose(VTbuf);
    s = sbuf;
  }

  /// Compute the singular value decomposition of the matrix A,
  /// including complete orthogonal bases of the domain and range even
  /// when A is rectangular.
  ///
  /// U, S, and VT will be resized if necessary to the appropriate output 
  /// dimensions based on the dimensions of A.
  template <class AMatrixT, class UMatrixT, class SingularValuesT, class VTMatrixT> 
  inline void complete_svd( AMatrixT & A, UMatrixT &U, SingularValuesT &s, VTMatrixT &VT ) {
    typedef typename PromoteType<typename AMatrixT::value_type, typename SingularValuesT::value_type>::type temp_type1;
    typedef typename PromoteType<temp_type1, typename UMatrixT::value_type>::type temp_type2;
    typedef typename PromoteType<temp_type2, typename VTMatrixT::value_type>::type real_type;
    const long int m = A.rows(), n = A.cols();
    const long int minmn = std::min(m,n);
    const long int lda = A.rows();
    Matrix<real_type> Abuf = transpose(A);
    Matrix<real_type> Ubuf( A.rows(), A.rows() );
    Matrix<real_type> VTbuf( A.cols(), A.cols() );
    Vector<real_type> sbuf( minmn );
    real_type work_size;
    long int lwork = -1, info;
    std::vector<long int> iwork( 8*minmn ); 
    long int ldu = m, ldvt = n;
    gesdd('A', m, n, &(Abuf(0,0)), lda, &(sbuf(0)), &(Ubuf(0,0)), ldu, &(VTbuf(0,0)), ldvt, &work_size, lwork, &iwork[0], &info);
    lwork = (long int)(work_size);
    std::vector<real_type> work( lwork );
    gesdd('A', m, n, &(Abuf(0,0)), lda, &(sbuf(0)), &(Ubuf(0,0)), ldu, &(VTbuf(0,0)), ldvt, &work[0], lwork, &iwork[0], &info);
    if (info < 0) 
      vw_throw( ArgumentErr() << "svd(): LAPACK driver gesdd reported an error with argument " << -info << "." );
    if (info > 0) 
      vw_throw( ArgumentErr() << "svd(): LAPACK driver gesdd did not converge.  Update process failed." );
    U = transpose(Ubuf);
    VT = transpose(VTbuf);
    s = sbuf;
  }


  /// Computes the pseudoinverse A* of a real matrix A.
  template <class AMatrixT>
  Matrix<typename AMatrixT::value_type> pseudoinverse( AMatrixT & A, double cond = 0 ) {
    Matrix<typename AMatrixT::value_type> u, vt;
    Vector<typename AMatrixT::value_type> s;
    svd( A, u, s, vt );
    Matrix<typename AMatrixT::value_type> si(s.size(),s.size());
    for( unsigned i=0; i<s.size(); ++i ) {
      for( unsigned j=0; j<s.size(); ++j ) {
        if( i==j ) si(i,j) = ( fabs(s(i)) <= cond*s(0) ) ? 0 : 1/s(i);
        else si(i,j) = 0;
      }
    }
    return transpose(vt)*si*transpose(u);
  }


  /// \cond INTERNAL
  extern "C"  int sgelsd_(const long int* m, const long int* n, const long int* nrhs, float* a, 
                          const long int* lda, float *b, const long int* ldb, float* s, 
                          const float* rcond, long int* rank, float* work, const long int* lwork, 
                          long int* iwork, long int* info);

  extern "C"  int dgelsd_(const long int* m, const long int* n, const long int* nrhs, double* a, 
                          const long int* lda, double *b, const long int* ldb, double* s, 
                          const double* rcond, long int* rank, double* work, const long int* lwork, 
                          long int* iwork, long int* info);

  static inline void gelsd(long int m, long int n, long int nrhs, float *a, long int lda, float *b, long int ldb, float *s, float rcond, long int *rank, float *work, long int lwork, long int* iwork, long int *info) {
    sgelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, iwork, info );
  }

  static inline void gelsd(long int m, long int n, long int nrhs, double *a, long int lda, double *b, long int ldb, double *s, double rcond, long int *rank, double *work, long int lwork, long int* iwork, long int *info) {
    dgelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, iwork, info );
  }
  /// \endcond

  /// Computes the minimum-norm solution to a real linear least squares problem.
  template <class AMatrixT, class BVectorT>
  Vector<typename PromoteType<typename AMatrixT::value_type, typename BVectorT::value_type>::type>
  least_squares( AMatrixT & A, BVectorT & B, double cond = -1 ) {
    typedef typename PromoteType<typename AMatrixT::value_type, typename BVectorT::value_type>::type real_type;
    Matrix<real_type> Abuf = transpose(A);
    const long int m = A.rows(), n = A.cols();
    const int minmn = std::min(m,n), maxmn = std::max(m,n);
    Vector<real_type> Bbuf(maxmn);
    subvector(Bbuf,0,m) = B;
    const long int nrhs = 1, lda = A.rows(), ldb = maxmn;
    std::vector<real_type> s( minmn );
    real_type const rcond = cond;
    long int rank, lwork = -1, info;
    std::vector<long int> iwork( (3*int(log(minmn+1.)/log(2.))+11)*minmn ); // log2(x) = log(x)/log(2)
    real_type work_size;
    gelsd( m, n, nrhs, &(Abuf(0,0)), lda, &(Bbuf(0)), ldb, &s[0], rcond, &rank, &work_size, lwork, &iwork[0], &info );
    lwork = int(work_size);
    std::vector<real_type> work( lwork );
    gelsd( m, n, nrhs, &(Abuf(0,0)), lda, &(Bbuf(0)), ldb, &s[0], rcond, &rank, &work[0], lwork, &iwork[0], &info );
    Bbuf.set_size(n,true);
    return Bbuf;
  }

} // namespace math
} // namespace vw

#endif // __VW_MATH_LINEAR_ALGEBRA_H__
