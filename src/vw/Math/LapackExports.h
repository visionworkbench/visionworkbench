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


#ifndef __VW_MATH_LAPACK_EXPORTS_H__
#define __VW_MATH_LAPACK_EXPORTS_H__

#include <vw/vw_config.h>
#include <vw/Core/FundamentalTypes.h>

#if (defined(VW_HAVE_PKG_APPLE_LAPACK) && VW_HAVE_PKG_APPLE_LAPACK==1)
  #include <Accelerate/Accelerate.h>
#endif


namespace vw {
/// Numerical linear algebra and computational geometry.
namespace math {

// There are seemingly two bloodlines of lapack. One decends directly from the
// fortran libs, and one is a c library formed by running f2c on the fortran
// code. The former uses the fortran definition of integer, an int32. The
// latter uses the f2c.h definition, long. The only reasonably safe thing to do
// (lacking standard installation of lapack headers) is to use the type most
// appropriate for the bloodline.

#if (defined(VW_HAVE_PKG_INTEL_LAPACK) && VW_HAVE_PKG_INTEL_LAPACK==1)
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
#elif (defined(VW_HAVE_PKG_APPLE_LAPACK) && VW_HAVE_PKG_APPLE_LAPACK==1)
//#  include <Accelerate/Accelerate.h>
   typedef __CLPK_integer f77_int;
#else
#  if (defined(VW_HAVE_PKG_FLAPACK)           ) || \
      (defined(VW_HAVE_PKG_SLAPACK)           ) || \
      (defined(VW_HAVE_PKG_STANDALONE_FLAPACK))

    //#pragma message ( "Using Fortran based LAPACK!" )

    // fortran-based
    typedef int32  f77_int;

#  elif (defined(VW_HAVE_PKG_CLAPACK)                   ) || \
        (defined(VW_HAVE_PKG_STANDALONE_LAPACK_AND_BLAS)) || \
        (defined(VW_HAVE_PKG_LAPACK))

    //#pragma message ( "Using C based LAPACK!" )

    // f2c-based
    typedef long f77_int;

#   else
// We're only installed if lapack is installed, so something is wrong at this point.
#   error I do not know what lapack you are using
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


    void sgesvd_(char    *jobu,  char    *jobvt, f77_int *m,    f77_int *n,
                 float   *a,     f77_int *lda,   float   *s,    float   *u,
                 f77_int *ldu,   float   *vt,    f77_int *ldvt, float   *work,
                 f77_int *lwork, f77_int *info);

    void dgesvd_(char    *jobu,  char    *jobvt, f77_int *m,    f77_int *n,
                 double  *a,     f77_int *lda,   double  *s,    double  *u,
                 f77_int *ldu,   double  *vt,    f77_int *ldvt, double  *work,
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

    void sgetrf_(f77_int *m, f77_int *n, float *a, f77_int *lda,
                 f77_int *ipiv, f77_int *info );
    void dgetrf_(f77_int *m, f77_int *n, double *a, f77_int *lda,
                 f77_int *ipiv, f77_int *info );

    void sgetrs_(char *trans, f77_int* n, f77_int *nrhs, float *a, f77_int *lda,
                 f77_int *ipiv, float *b, f77_int *ldb, f77_int *info );
    void dgetrs_(char *trans, f77_int* n, f77_int *nrhs, double *a, f77_int *lda,
                 f77_int *ipiv, double *b, f77_int *ldb, f77_int *info );
  } // end extern "C"
  } // anon namespace
#endif

  void geev(char jobvl, char jobvr, f77_int n, float *a, f77_int lda, float *wr, float *wi, float *vl, f77_int ldvl, float *vr, f77_int ldvr, float *work, f77_int lwork, f77_int *info);
  void geev(char jobvl, char jobvr, f77_int n, double *a, f77_int lda, double *wr, double *wi, double *vl, f77_int ldvl, double *vr, f77_int ldvr, double *work, f77_int lwork, f77_int *info);

  void gesvd(char    jobu,  char    jobvt, f77_int m,    f77_int n,
             float   *a,     f77_int lda,   float   *s,    float   *u,
             f77_int ldu,   float   *vt,    f77_int ldvt, float   *work,
             f77_int lwork, f77_int *info);
  void gesvd(char    jobu,  char    jobvt, f77_int m,    f77_int n,
             double  *a,     f77_int lda,   double  *s,    double  *u,
             f77_int ldu,   double  *vt,    f77_int ldvt, double  *work,
             f77_int lwork, f77_int *info);

  void gesdd(char jobz, f77_int m, f77_int n, float *a, f77_int lda, float *s, float *u, f77_int ldu, float *vt, f77_int ldvt, float *work, f77_int lwork, f77_int *iwork, f77_int *info);
  void gesdd(char jobz, f77_int m, f77_int n, double *a, f77_int lda, double *s, double *u, f77_int ldu, double *vt, f77_int ldvt, double *work, f77_int lwork, f77_int *iwork, f77_int *info);

  void geqrf(f77_int m, f77_int n, float *a, f77_int lda, float *tau, float *work, f77_int lwork, f77_int *info);
  void geqrf(f77_int m, f77_int n, double *a, f77_int lda, double *tau, double *work, f77_int lwork, f77_int *info);

  void orgqr(f77_int m, f77_int n, f77_int k, float* a, f77_int lda, float* tau, float* work, f77_int lwork, f77_int* info);
  void orgqr(f77_int m, f77_int n, f77_int k, double* a, f77_int lda, double* tau, double* work, f77_int lwork, f77_int* info);

  void gerqf(f77_int m, f77_int n, float *a, f77_int lda, float *tau, float *work, f77_int lwork, f77_int *info);
  void gerqf(f77_int m, f77_int n, double *a, f77_int lda, double *tau, double *work, f77_int lwork, f77_int *info);

  void orgrq(f77_int m, f77_int n, f77_int k, float* a, f77_int lda, float* tau, float* work, f77_int lwork, f77_int* info);
  void orgrq(f77_int m, f77_int n, f77_int k, double* a, f77_int lda, double* tau, double* work, f77_int lwork, f77_int* info);

  void gelsd(f77_int m, f77_int n, f77_int nrhs, float *a, f77_int lda, float *b, f77_int ldb, float *s, float rcond, f77_int *rank, float *work, f77_int lwork, f77_int* iwork, f77_int *info);
  void gelsd(f77_int m, f77_int n, f77_int nrhs, double *a, f77_int lda, double *b, f77_int ldb, double *s, double rcond, f77_int *rank, double *work, f77_int lwork, f77_int* iwork, f77_int *info);

  void gelss(f77_int m, f77_int n, f77_int rhs, float *a, f77_int lda, float *b, f77_int ldb, float *s, float rcond, f77_int *rank, float *work, f77_int lwork, f77_int *info);
  void gelss(f77_int m, f77_int n, f77_int rhs, double *a, f77_int lda, double *b, f77_int ldb, double *s, double rcond, f77_int *rank, double *work, f77_int lwork, f77_int *info);

  void gesv(f77_int n, f77_int nrhs, float *a, f77_int lda, f77_int *ipiv, float *b, f77_int ldb, f77_int *info);
  void gesv(f77_int n, f77_int nrhs, double *a, f77_int lda, f77_int *ipiv, double *b, f77_int ldb, f77_int *info);

  void posv(char uplo, f77_int n, f77_int nrhs, float *a, f77_int lda, float *b, f77_int ldb, f77_int *info);
  void posv(char uplo, f77_int n, f77_int nrhs, double *a, f77_int lda, double *b, f77_int ldb, f77_int *info);

  void getrf(f77_int m, f77_int n, float *a, f77_int lda, f77_int *ipiv, f77_int *info);
  void getrf(f77_int m, f77_int n, double *a, f77_int lda, f77_int *ipiv, f77_int *info);

  void getrs(char trans, f77_int n, f77_int nrhs, float *a, f77_int lda, f77_int *ipiv, float *b, f77_int ldb, f77_int *info);
  void getrs(char trans, f77_int n, f77_int nrhs, double *a, f77_int lda, f77_int *ipiv, double *b, f77_int ldb, f77_int *info);

  /// \endcond

} // namespace math
} // namespace vw


#endif

