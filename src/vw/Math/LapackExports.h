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

#if (defined(VW_HAVE_PKG_FLAPACK)            && VW_HAVE_PKG_FLAPACK==1) || \
    (defined(VW_HAVE_PKG_SLAPACK)            && VW_HAVE_PKG_SLAPACK==1) || \
    (defined(VW_HAVE_PKG_STANDALONE_FLAPACK) && VW_HAVE_PKG_STANDALONE_FLAPACK==1)

  // fortran-based
  typedef int32  f77_int;

#elif (defined(VW_HAVE_PKG_CLAPACK)                    && VW_HAVE_PKG_CLAPACK==1) || \
      (defined(VW_HAVE_PKG_STANDALONE_LAPACK_AND_BLAS) && VW_HAVE_PKG_STANDALONE_LAPACK_AND_BLAS==1) || \
      (defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1)
  

  // f2c-based
  typedef long f77_int;

#else
// We're only installed if lapack is installed, so something is wrong at this point.
#error I don't know what lapack you're using
#endif

  /// \cond INTERNAL
  namespace {
  extern "C" {

    int sgeev_(char    *jobvl, char *jobvr, f77_int    *n, float    *a,
               f77_int  *lda, float   *wr, float      *wi, float   *vl,
               f77_int *ldvl, float   *vr, f77_int *ldvr, float *work,
               f77_int *lwork, f77_int *info);

    int dgeev_(char    *jobvl, char *jobvr, f77_int    *n, double    *a,
               f77_int  *lda, double  *wr, double     *wi, double   *vl,
               f77_int *ldvl, double  *vr, f77_int *ldvr, double *work,
               f77_int *lwork, f77_int *info);


    int sgesdd_(char     *jobz, f77_int    *m, f77_int *n, float    *a,
                f77_int  *lda, float       *s, float    *u, f77_int *ldu,
                float      *vt, f77_int *ldvt, float *work, f77_int  *lwork,
                f77_int *iwork, f77_int *info);

    int dgesdd_(char     *jobz, f77_int    *m, f77_int  *n, double *a,
                f77_int  *lda, double      *s, double    *u, f77_int *ldu,
                double     *vt, f77_int *ldvt, double *work, f77_int *lwork,
                f77_int *iwork, f77_int *info);


    int sgelsd_(const f77_int   *m, const f77_int   *n, const f77_int  *nrhs,
                float            *a, const f77_int *lda, float              *b,
                const f77_int *ldb, float            *s, const float    *rcond,
                f77_int      *rank, float         *work, const f77_int  *lwork,
                f77_int      *iwork, f77_int       *info);

    int dgelsd_(const f77_int* m, const f77_int* n, const f77_int* nrhs, double* a,
                const f77_int* lda, double *b, const f77_int* ldb, double* s,
                const double* rcond, f77_int* rank, double* work, const f77_int* lwork,
                f77_int* iwork, f77_int* info);


    int sgelss_(const f77_int* m, const f77_int* n, const f77_int* nrhs, float* a,
                const f77_int* lda, float *b, const f77_int* ldb, float* s,
                const float* rcond, f77_int* rank, float* work, const f77_int* lwork,
                f77_int* info);

    int dgelss_(const f77_int* m, const f77_int* n, const f77_int* nrhs, double* a,
                const f77_int* lda, double *b, const f77_int* ldb, double* s,
                const double* rcond, f77_int* rank, double* work, const f77_int* lwork,
                f77_int* info);
    int sgesv_(const f77_int *n, const f77_int *nrhs, float *a,
               const f77_int *lda, f77_int *ipiv, float *b,
                 const f77_int *ldb, f77_int *info);

    int dgesv_(const f77_int *n, const f77_int *nrhs, double *a,
               const f77_int *lda, f77_int *ipiv, double *b,
               const f77_int *ldb, f77_int *info);
    int sposv_(const char *uplo, const f77_int *n, const f77_int *nrhs, float *a,
               const f77_int *lda, float *b, const f77_int *ldb, f77_int *info);

    int dposv_(const char *uplo, const f77_int *n, const f77_int *nrhs, double *a,
                           const f77_int *lda, double *b, const f77_int *ldb, f77_int *info);
  }
  } // anon namespace

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

  static inline void gesv(const f77_int n, const f77_int nrhs, float *a, const f77_int lda, f77_int *ipiv, float *b, const f77_int ldb, f77_int *info) {
    sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, info);
  }

  static inline void gesv(const f77_int n, const f77_int nrhs, double *a, const f77_int lda, f77_int *ipiv, double *b, const f77_int ldb, f77_int *info) {
    dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, info);
  }

  static inline void posv(const char uplo, const f77_int n, const f77_int nrhs, float *a, const f77_int lda, float *b, const f77_int ldb, f77_int *info) {
    sposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, info);
  }

  static inline void posv(const char uplo, const f77_int n, const f77_int nrhs, double *a, const f77_int lda, double *b, const f77_int ldb, f77_int *info) {
    dposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, info);
  }

  /// \endcond

} // namespace math
} // namespace vw


#endif
