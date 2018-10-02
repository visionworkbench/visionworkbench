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


#include <vw/Math/LapackExports.h>
#include <vw/Core/Debugging.h>
#include <vw/Core/Exception.h>

#define CHECK() vw::math::detail::_check_info(info, VW_CURRENT_FUNCTION, __FILE__, __LINE__)

namespace vw {
namespace math {


  namespace detail {
    void _check_info(const f77_int *info, const char* func, const char* file, const int line) {
      if (*info < 0)
        vw::vw_throw( vw::ArgumentErr() << file << ":" << line << "LAPACK reported an error with argument " << -(*info) << " in " << func);
    }
  }

#define D_geev(prefix, type)\
  void geev(char jobvl, char jobvr, f77_int n, type *a, f77_int lda, type *wr, type *wi, type *vl, f77_int ldvl, type *vr, f77_int ldvr, type *work, f77_int lwork, f77_int *info) { \
    prefix ## geev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, info); \
  }

D_geev(s, float);
D_geev(d, double);

#define D_gesvd(prefix, type)\
  void gesvd(char jobu, char jobvt, f77_int m, f77_int n, type *a, f77_int lda, type *s, type *u, f77_int ldu, type *vt, f77_int ldvt, type *work, f77_int lwork, f77_int *info) {\
    prefix ## gesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, info); \
    CHECK(); \
  }

D_gesvd(s, float);
D_gesvd(d, double);

#define D_gesdd(prefix, type)\
  void gesdd(char jobz, f77_int m, f77_int n, type *a, f77_int lda, type *s, type *u, f77_int ldu, type *vt, f77_int ldvt, type *work, f77_int lwork, f77_int *iwork, f77_int *info) {\
    prefix ## gesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, info); \
    CHECK(); \
  }

D_gesdd(s, float);
D_gesdd(d, double);

#define D_geqrf(prefix, type)\
  void geqrf(f77_int m, f77_int n, type *a, f77_int lda, type *tau, type *work, f77_int lwork, f77_int *info) {\
    prefix ## geqrf_(&m, &n, a, &lda, tau, work, &lwork, info); \
    CHECK(); \
  }

D_geqrf(s, float);
D_geqrf(d, double);

#define D_orgqr(prefix, type)\
  void orgqr(f77_int m, f77_int n, f77_int k, type* a, f77_int lda, type* tau, type* work, f77_int lwork, f77_int* info) { \
    prefix ## orgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, info); \
    CHECK(); \
  }

D_orgqr(s, float);
D_orgqr(d, double);

#define D_gerqf(prefix, type)\
  void gerqf(f77_int m, f77_int n, type *a, f77_int lda, type *tau, type *work, f77_int lwork, f77_int *info) { \
    prefix ## gerqf_(&m, &n, a, &lda, tau, work, &lwork, info); \
    CHECK(); \
  }

D_gerqf(s, float);
D_gerqf(d, double);

#define D_orgrq(prefix, type)\
  void orgrq(f77_int m, f77_int n, f77_int k, type* a, f77_int lda, type* tau, type* work, f77_int lwork, f77_int* info) { \
    prefix ## orgrq_(&m, &n, &k, a, &lda, tau, work, &lwork, info); \
    CHECK(); \
  }

D_orgrq(s, float);
D_orgrq(d, double);

#define D_gelsd(prefix, type)\
  void gelsd(f77_int m, f77_int n, f77_int nrhs, type *a, f77_int lda, type *b, f77_int ldb, type *s, type rcond, f77_int *rank, type *work, f77_int lwork, f77_int* iwork, f77_int *info) { \
    prefix ## gelsd_( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, iwork, info ); \
    CHECK(); \
  }

D_gelsd(s, float);
D_gelsd(d, double);

#define D_gelss(prefix, type)\
  void gelss(f77_int m, f77_int n, f77_int rhs, type *a, f77_int lda, type *b, f77_int ldb, type *s, type rcond, f77_int *rank, type *work, f77_int lwork, f77_int *info) { \
    prefix ## gelss_( &m, &n, &rhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, info ); \
    CHECK(); \
  }

D_gelss(s, float);
D_gelss(d, double);

#define D_gesv(prefix, type)\
  void gesv(f77_int n, f77_int nrhs, type *a, f77_int lda, f77_int *ipiv, type *b, f77_int ldb, f77_int *info) { \
    prefix ## gesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, info); \
    CHECK(); \
  }

D_gesv(s, float);
D_gesv(d, double);

#define D_posv(prefix, type)\
  void posv(char uplo, f77_int n, f77_int nrhs, type *a, f77_int lda, type *b, f77_int ldb, f77_int *info) { \
    prefix ## posv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, info);          \
    CHECK(); \
  }

D_posv(s, float);
D_posv(d, double);

#define D_getrf(prefix, type)\
  void getrf(f77_int m, f77_int n, type *a, f77_int lda, f77_int *ipiv, f77_int *info ) { \
    prefix ## getrf_(&m, &n, a, &lda, ipiv, info );                     \
    CHECK();                                                            \
  }

D_getrf(s, float);
D_getrf(d, double);

#define D_getrs(prefix, type)                   \
  void getrs(char trans, f77_int n, f77_int nrhs, type *a, f77_int lda, f77_int *ipiv, type *b, f77_int ldb, f77_int *info ) { \
    prefix ## getrs_(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, info ); \
    CHECK();                                                            \
  }

D_getrs(s, float);
D_getrs(d, double);

}} // vw::math
