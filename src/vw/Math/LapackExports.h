// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_MATH_LAPACK_EXPORTS_H__
#define __VW_MATH_LAPACK_EXPORTS_H__

#include <vw/config.h>

namespace vw {
/// Numerical linear algebra and computational geometry.
namespace math {

// There are seemingly two bloodlines of lapack. One decends directly from the
// fortran libs, and one is a c library formed by running f2c on the fortran
// code. The former uses the fortran definition of integer, an int32. The
// latter uses the f2c.h definition, long. The only reasonably safe thing to do
// (lacking standard installation of lapack headers) is to use the type most
// appropriate for the bloodline.

#if (defined(VW_HAVE_PKG_INTEL_LAPACK)       && VW_HAVE_PKG_INTEL_LAPACK==1)
#  include <mkl_lapack.h>
// MKL headers pollute the macro space. #define P4 breaks boost. Clean up after Intel.
#  undef ITP
#  undef NI
#  undef CT
#  undef MNI
#  undef PNR
#  undef DEF
#  undef PIII
#  undef P4
#  undef P4P
#  undef P4M
   typedef MKL_INT f77_int;
#else
#  if (defined(VW_HAVE_PKG_FLAPACK)            && VW_HAVE_PKG_FLAPACK==1) || \
      (defined(VW_HAVE_PKG_SLAPACK)            && VW_HAVE_PKG_SLAPACK==1) || \
      (defined(VW_HAVE_PKG_STANDALONE_FLAPACK) && VW_HAVE_PKG_STANDALONE_FLAPACK==1)

    // fortran-based
    typedef int32  f77_int;

#  elif (defined(VW_HAVE_PKG_CLAPACK)                    && VW_HAVE_PKG_CLAPACK==1) || \
        (defined(VW_HAVE_PKG_STANDALONE_LAPACK_AND_BLAS) && VW_HAVE_PKG_STANDALONE_LAPACK_AND_BLAS==1) || \
        (defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1)

    // f2c-based
    typedef long f77_int;

#   else
// We're only installed if lapack is installed, so something is wrong at this point.
#   error I don't know what lapack you're using
#   endif

  /// \cond INTERNAL
  namespace {
  extern "C" {

    void sgeev_(char    *jobvl, char    *jobvr, f77_int *n,    float *a,
                f77_int *lda,   float   *wr,    float   *wi,   float *vl,
                f77_int *ldvl,  float   *vr,    f77_int *ldvr, float *work,
                f77_int *lwork, f77_int *info);

    void dgeev_(char    *jobvl, char    *jobvr, f77_int *n,    double *a,
                f77_int *lda,   double  *wr,    double  *wi,   double *vl,
                f77_int *ldvl,  double  *vr,    f77_int *ldvr, double *work,
                f77_int *lwork, f77_int *info);


    void sgesdd_(char    *jobz,  f77_int *m,    f77_int *n,    float   *a,
                 f77_int *lda,   float   *s,    float   *u,    f77_int *ldu,
                 float   *vt,    f77_int *ldvt, float   *work, f77_int *lwork,
                 f77_int *iwork, f77_int *info);

    void dgesdd_(char    *jobz,  f77_int *m,    f77_int *n,    double *a,
                 f77_int *lda,   double  *s,    double  *u,    f77_int *ldu,
                 double  *vt,    f77_int *ldvt, double  *work, f77_int *lwork,
                 f77_int *iwork, f77_int *info);


    void sgeqrf_ (f77_int *m,   f77_int *n,    float   *a,     f77_int *lda,
                  float   *tau, float   *work, f77_int *lwork, f77_int *info);

    void dgeqrf_ (f77_int *m,   f77_int *n,    double  *a,     f77_int *lda,
                  double  *tau, double  *work, f77_int *lwork, f77_int *info);

    void sorgqr_(f77_int *m,   f77_int *n,    f77_int *k,     float   *a,   f77_int* lda,
                 float   *tau, float   *work, f77_int *lwork, f77_int *info);

    void dorgqr_(f77_int* m, f77_int* n, f77_int* k, double* a, f77_int* lda,
                 double* tau, double* work, f77_int* lwork, f77_int* info);


    void sgerqf_ (f77_int *m, f77_int *n, float *a, f77_int *lda,
                 float *tau, float *work, f77_int *lwork, f77_int *info);

    void dgerqf_ (f77_int *m, f77_int *n, double *a, f77_int *lda,
                 double *tau, double *work, f77_int *lwork, f77_int *info);

    void sorgrq_(f77_int* m, f77_int* n, f77_int* k, float* a, f77_int* lda,
                float* tau, float* work, f77_int* lwork, f77_int* info);

    void dorgrq_(f77_int* m, f77_int* n, f77_int* k, double* a, f77_int* lda,
                double* tau, double* work, f77_int* lwork, f77_int* info);


    void sgelsd_(f77_int   *m, f77_int   *n, f77_int  *nrhs,
                float      *a, f77_int *lda, float              *b,
                f77_int *ldb, float            *s, float    *rcond,
                f77_int      *rank, float         *work, f77_int  *lwork,
                f77_int      *iwork, f77_int       *info);

    void dgelsd_(f77_int* m, f77_int* n, f77_int* nrhs, double* a,
                f77_int* lda, double *b, f77_int* ldb, double* s,
                double* rcond, f77_int* rank, double* work, f77_int* lwork,
                f77_int* iwork, f77_int* info);


    void sgelss_(f77_int* m, f77_int* n, f77_int* nrhs, float* a,
                f77_int* lda, float *b, f77_int* ldb, float* s,
                float* rcond, f77_int* rank, float* work, f77_int* lwork,
                f77_int* info);

    void dgelss_(f77_int* m, f77_int* n, f77_int* nrhs, double* a,
                f77_int* lda, double *b, f77_int* ldb, double* s,
                double* rcond, f77_int* rank, double* work, f77_int* lwork,
                f77_int* info);
    void sgesv_(f77_int *n, f77_int *nrhs, float *a,
               f77_int *lda, f77_int *ipiv, float *b,
                 f77_int *ldb, f77_int *info);

    void dgesv_(f77_int *n, f77_int *nrhs, double *a,
               f77_int *lda, f77_int *ipiv, double *b,
               f77_int *ldb, f77_int *info);
    void sposv_(char *uplo, f77_int *n, f77_int *nrhs, float *a,
               f77_int *lda, float *b, f77_int *ldb, f77_int *info);

    void dposv_(char *uplo, f77_int *n, f77_int *nrhs, double *a,
                           f77_int *lda, double *b, f77_int *ldb, f77_int *info);
  }
  } // anon namespace
#endif

  static inline void geev(char jobvl, char jobvr, f77_int n, float *a, f77_int lda, float *wr, float *wi, float *vl, f77_int ldvl, float *vr, f77_int ldvr, float *work, f77_int lwork, f77_int *info) {
    sgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, info);
  }

  static inline void geev(char jobvl, char jobvr, f77_int n, double *a, f77_int lda, double *wr, double *wi, double *vl, f77_int ldvl, double *vr, f77_int ldvr, double *work, f77_int lwork, f77_int *info) {
    dgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, info);
  }

  static inline void gesdd(char jobz, f77_int m, f77_int n, float *a, f77_int lda, float *s, float *u, f77_int ldu, float *vt, f77_int ldvt, float *work, f77_int lwork, f77_int *iwork, f77_int *info) {
    sgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, info);
  }

  static inline void gesdd(char jobz, f77_int m, f77_int n, double *a, f77_int lda, double *s, double *u, f77_int ldu, double *vt, f77_int ldvt, double *work, f77_int lwork, f77_int *iwork, f77_int *info) {
    dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, info);
  }

  static inline void sgeqrf(f77_int m, f77_int n, float *a, f77_int lda,
                            float *tau, float *work, f77_int lwork, f77_int *info) {
    sgeqrf_(&m, &n, a, &lda, tau, work, &lwork, info);
  }

  static inline void sgeqrf(f77_int m, f77_int n, double *a, f77_int lda,
                            double *tau, double *work, f77_int lwork, f77_int *info) {
    dgeqrf_(&m, &n, a, &lda, tau, work, &lwork, info);
  }


  static inline void sorgqr(f77_int m, f77_int n, f77_int k, float* a, f77_int lda,
                            float* tau, float* work, f77_int lwork, f77_int* info) {
    sorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, info);
  }

  static inline void sorgqr(f77_int m, f77_int n, f77_int k, double* a, f77_int lda,
                            double* tau, double* work, f77_int lwork, f77_int* info) {
    dorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, info);
  }

  static inline void sgerqf(f77_int m, f77_int n, float *a, f77_int lda,
                            float *tau, float *work, f77_int lwork, f77_int *info) {
    sgerqf_(&m, &n, a, &lda, tau, work, &lwork, info);
  }

  static inline void sgerqf(f77_int m, f77_int n, double *a, f77_int lda,
                            double *tau, double *work, f77_int lwork, f77_int *info) {
    dgerqf_(&m, &n, a, &lda, tau, work, &lwork, info);
  }

  static inline void sorgrq(f77_int m, f77_int n, f77_int k, float* a, f77_int lda,
                            float* tau, float* work, f77_int lwork, f77_int* info) {
    sorgrq_(&m, &n, &k, a, &lda, tau, work, &lwork, info);
  }

  static inline void sorgrq(f77_int m, f77_int n, f77_int k, double* a, f77_int lda,
                            double* tau, double* work, f77_int lwork, f77_int* info) {
    dorgrq_(&m, &n, &k, a, &lda, tau, work, &lwork, info);
  }

  static inline void gelsd(f77_int m, f77_int n, f77_int nrhs, float *a, f77_int lda, float *b, f77_int ldb, float *s, float rcond, f77_int *rank, float *work, f77_int lwork, f77_int* iwork, f77_int *info) {
    sgelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, iwork, info );
  }

  static inline void gelsd(f77_int m, f77_int n, f77_int nrhs, double *a, f77_int lda, double *b, f77_int ldb, double *s, double rcond, f77_int *rank, double *work, f77_int lwork, f77_int* iwork, f77_int *info) {
    dgelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, iwork, info );
  }

  static inline void gelss(f77_int m, f77_int n, f77_int rhs, double *a, f77_int lda, double *b, f77_int ldb, double *s, double rcond, f77_int *rank, double *work, f77_int lwork, f77_int *info) {
    dgelss_( &m, &n, &rhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, info );
  }

  static inline void gelss(f77_int m, f77_int n, f77_int rhs, float *a, f77_int lda, float *b, f77_int ldb, float *s, float rcond, f77_int *rank, float *work, f77_int lwork, f77_int *info) {
    sgelss_( &m, &n, &rhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, info );
  }

  static inline void gesv(f77_int n, f77_int nrhs, float *a, f77_int lda, f77_int *ipiv, float *b, f77_int ldb, f77_int *info) {
    sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, info);
  }

  static inline void gesv(f77_int n, f77_int nrhs, double *a, f77_int lda, f77_int *ipiv, double *b, f77_int ldb, f77_int *info) {
    dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, info);
  }

  static inline void posv(char uplo, f77_int n, f77_int nrhs, float *a, f77_int lda, float *b, f77_int ldb, f77_int *info) {
    sposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, info);
  }

  static inline void posv(char uplo, f77_int n, f77_int nrhs, double *a, f77_int lda, double *b, f77_int ldb, f77_int *info) {
    dposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, info);
  }

  /// \endcond

} // namespace math
} // namespace vw


#endif
