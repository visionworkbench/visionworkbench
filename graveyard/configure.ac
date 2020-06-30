dnl __BEGIN_LICENSE__
dnl  Copyright (c) 2006-2013, United States Government as represented by the
dnl  Administrator of the National Aeronautics and Space Administration. All
dnl  rights reserved.
dnl
dnl  The NASA Vision Workbench is licensed under the Apache License,
dnl  Version 2.0 (the "License"); you may not use this file except in
dnl  compliance with the License. You may obtain a copy of the License at
dnl  http://www.apache.org/licenses/LICENSE-2.0
dnl
dnl  Unless required by applicable law or agreed to in writing, software
dnl  distributed under the License is distributed on an "AS IS" BASIS,
dnl  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
dnl  See the License for the specific language governing permissions and
dnl  limitations under the License.
dnl __END_LICENSE__


m4_pattern_forbid([^AX_])
m4_pattern_allow([^AX_CFLAGS$])
m4_pattern_allow([^AX_INCLUDE_DIR$])
m4_pattern_allow([^AX_LIBDIR$])
m4_pattern_allow([^AX_OTHER_LIBDIR$])
m4_pattern_allow([^PKG_[A-Z_]+_CPPFLAGS$])
m4_pattern_allow([^PKG_[A-Z_]+_LDFLAGS$])
m4_pattern_allow([^PKG_[A-Z_]+_LIBS$])
m4_pattern_allow([^PKG_PATHS$])
m4_pattern_allow([^PKG_PATHS_[A-Z_]+$])
m4_pattern_allow([^PKG_BOOST_INCDIR$])
m4_pattern_allow([^PKG_BOOST_LIBDIR$])
m4_pattern_allow([^PKG_VW_LIBS_TEMP$])

##################################################
# A few starting bits
##################################################

AC_PREREQ([2.59])
dnl TODO BEFORE/AFTER RELEASE: update package version
AC_INIT([NASA Vision Workbench], [2.2.0_post], [vision-workbench@lists.nasa.gov], [VisionWorkbench])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE([nostdinc 1.8.5])
AX_CONFIG_HEADER_PREFIX([src/vw/config.h], [VW_])

dnl Fix for early autoconf that doesn't export these
AC_SUBST([abs_top_srcdir])
AC_SUBST([abs_top_builddir])
AC_SUBST([abs_srcdir])
AC_SUBST([abs_builddir])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

m4_divert_text([INIT_PREPARE], [dnl
# load the config.options
if test -f ./config.options; then
  AX_LOG([using config.options << EOF])
  AX_LOG([`cat ./config.options`])
  AX_LOG([EOF])
  . ./config.options
fi
])

##################################################
# Find the necessary programs
##################################################

# Make sure these are set so AC_PROG_CC doesn't initialize them with
# its defaults, since we take care of setting debug and optimization
# levels ourselves.
CFLAGS="$CFLAGS"
CXXFLAGS="$CXXFLAGS"

# Disable static libraries by default
AM_DISABLE_STATIC

# Find the tools
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AX_PROG_AR
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET
AC_LANG(C++)
AX_COMMON_OPTIONS
AC_CHECK_TYPE(ssize_t, int)

# Does the compiler support various fun function attributes?
AX_CHECK_FUNCTION_ATTRIBUTE([noreturn])
AX_CHECK_FUNCTION_ATTRIBUTE([deprecated])
AX_CHECK_FUNCTION_ATTRIBUTE([warn_unused_result])

# Looking for headers
AC_CHECK_HEADERS([unistd.h pwd.h fenv.h ext/stdio_filebuf.h])

# Find some functions
AC_SEARCH_LIBS([mkstemps], [iberty])
AC_CHECK_FUNCS([mkstemps getpwuid getpid feenableexcept])
AC_CHECK_FUNCS([fabsl llabs exp2 log2 tgamma])

# Look for the python stuff
AM_PATH_PYTHON([2.4], [HAVE_PYTHON=yes], [HAVE_PYTHON=no])
AC_PYTHON_DEVEL([HAVE_PYTHON=yes], [HAVE_PYTHON=no])
AC_PROG_SWIG([1.3.29], [HAVE_SWIG=yes], [HAVE_SWIG=no])

AC_SUBST(HAVE_PYTHON)
AC_SUBST(HAVE_SWIG)

if test "$HAVE_SWIG" = "yes"; then
    SWIG_ENABLE_CXX
    if test "$HAVE_PYTHON" = "yes"; then
        SWIG_PYTHON
    fi
fi

# Check for the code coverage programs
AC_CHECK_PROGS(LCOV, lcov, false)
AC_CHECK_PROGS(GENHTML, genhtml, false)

## Check for the Apache Extension Tool
#AC_CHECK_PROGS(APXS, apxs2 apxs, false)
#if test "$APXS" != "false"; then
#  APXS_CFLAGS="`$APXS -q CFLAGS` -I`$APXS -q INCLUDEDIR`"
#  APXS_LDFLAGS="`$APXS -q LDFLAGS`"
#  APR_CONFIG=`$APXS -q APR_CONFIG`
#  APU_CONFIG=`$APXS -q APU_CONFIG`
#  APXS_CFLAGS="$APXS_CFLAGS `$APR_CONFIG --cppflags --includes` `$APU_CONFIG --includes`"
#  APXS_LDFLAGS="$APXS_LDFLAGS `$APR_CONFIG --ldflags --libs` `$APU_CONFIG --ldflags --libs`"
#  APXS_DEFAULT_INSTALLDIR="`$APXS -q LIBEXECDIR`"

#  AX_ARG_WITH(apxs_installdir, [], [mk], [Set apache module install dir])

#  if test -z "$APXS_INSTALLDIR"; then
#    APXS_INSTALLDIR="$APXS_DEFAULT_INSTALLDIR"
#  fi

#  AC_SUBST(APXS_CFLAGS)
#  AC_SUBST(APXS_LDFLAGS)
#fi



##################################################
# Compilation options
##################################################

# For the moment this is just a manual switch, defaulting to off.
# Someday this should be sensibly autodetected.
AX_ARG_ENABLE(sse,           no,  [am-yes cpp-bool], [enable SSE optimizations in some places])
AX_ARG_ENABLE(bounds-check,  no,  [am-yes cpp-bool], [enable image bounds checking (SLOW!)])
AX_ARG_ENABLE(test-coverage, no,  [am-yes],          [turn on code coverage analysis tools])
AX_ARG_ENABLE(config-file,   yes, [cpp-bool],        [~/.vwrc support])
AX_ARG_ENABLE(strict-aliasing, no, [none],           [allow strict aliasing (dangerous)])


##################################################
# Handle options
##################################################

# Check for --enable-test-coverate and add appropriate flags for gcc
if test "$ENABLE_TEST_COVERAGE" = "yes"; then
  if CCACHE_CC=echo $CXX i_am_really_ccache 2>&AS_MESSAGE_LOG_FD | grep -q i_am_really_ccache ; then
    AC_MSG_ERROR([compiler CXX ($CXX) appears to be ccache, which breaks test-coverage])
  fi

  AX_TRY_CPPFLAGS([-fprofile-arcs],  [], [AC_MSG_ERROR([test coverage requires -fprofile-arcs])])
  AX_TRY_CPPFLAGS([-ftest-coverage], [], [AC_MSG_ERROR([test coverage requires -ftest-coverage])])

  echo "*** Enabling test coverage and debugging, and disabling optimization ***"

  dnl These two are necessary and already tested for
  CXXFLAGS="$CXXFLAGS -fprofile-arcs -ftest-coverage"
  LDFLAGS="$LDFLAGS -lgcov"
  ENABLE_OPTIMIZE=no
  ENABLE_DEBUG=yes
fi

# Files test-coverage and profile-arcs generate
AC_SUBST([MOSTLYCLEANFILES], "*.bb *.bbg *.da *.gcov *.gcda *.gcno .*.gcov")

# These need to be here because automake-1.6 is dumb and does not properly
# process AM_CONDITIONALs unless the first argument is a simple string.
AM_CONDITIONAL(ENABLE_SSE, [test x"$ENABLE_SSE" = x"yes"])
AM_CONDITIONAL(ENABLE_BOUNDS_CHECK, [test x"$ENABLE_BOUNDS_CHECK" = x"yes"])
AM_CONDITIONAL(ENABLE_TEST_COVERAGE, [test x"$ENABLE_TEST_COVERAGE" = x"yes"])

if test "$ENABLE_SSE" = "yes"; then
  echo "*** Enabling SSE support ***"
  dnl Add required flag
  CXXFLAGS="$CXXFLAGS -msse4.1"
fi

##################################################
# distribution options
##################################################

AX_ARG_WITH(dist-remove-modules, [Python Geometry], [mk am-set], [remove certain modules from the distribution])

# These need to be here because automake-1.6 is dumb and does not properly
# process AM_CONDITIONALs unless the first argument is a simple string.
AM_CONDITIONAL(DIST_REMOVE_MODULES, [test ! -z "$DIST_REMOVE_MODULES"])



##################################################
# package checks
##################################################

AX_PKG_PTHREADS

AX_PKG_BOOST([BOOST_PROGRAM_OPTIONS BOOST_FILESYSTEM BOOST_THREAD BOOST_GRAPH],
  BOOST_PROGRAM_OPTIONS,
    [AX_PKG_BOOST_LIB(PROGRAM_OPTIONS, [-lboost_program_options], [boost/program_options.hpp], [boost::program_options::options_description("Rawr");])],
  BOOST_SYSTEM,
    [AX_PKG_BOOST_LIB(SYSTEM,          [-lboost_system],          [boost/system/error_code.hpp], [boost::system::error_code();])],
  BOOST_FILESYSTEM,
    [AX_PKG_BOOST_LIB(FILESYSTEM,      [-lboost_filesystem],      [boost/filesystem/path.hpp], [boost::filesystem::path();], [$PKG_BOOST_SYSTEM_LIBS])],
  BOOST_REGEX,
    [AX_PKG_BOOST_LIB(REGEX,           [-lboost_regex],           [boost/regex.hpp], [boost::regex(".*");])],
  BOOST_DATETIME,
    [AX_PKG_BOOST_LIB(DATETIME,        [-lboost_date_time],       [boost/date_time/gregorian/gregorian_types.hpp], [], [$PKG_BOOST_SYSTEM_LIBS])],
  BOOST_THREAD,
    [AX_PKG_BOOST_LIB(THREAD,          [-lboost_thread],          [boost/thread/thread.hpp], [boost::mutex();], [$PKG_BOOST_DATETIME_LIBS], [boost/thread/condition.hpp])],
  BOOST_IOSTREAMS,
    [AX_PKG_BOOST_LIB(IOSTREAMS,       [-lboost_iostreams],       [boost/iostreams/device/mapped_file.hpp], [boost::iostreams::mapped_file_source();])],
  BOOST_GRAPH,
    [AX_PKG_BOOST_LIB(GRAPH,           [],                        [boost/graph/adjacency_list.hpp])]
)

if test x"$HAVE_PKG_BOOST" = "xyes"; then
  AX_PKG_BOOST_CHECK_VERSION([VW_])
  AC_MSG_NOTICE([Using Boost header version $BOOST_VERSION, library suffix $PKG_BOOST_SAME_SUFFIX])
fi

# This provides an easy way for users to override the threads options
# if they are replacing Boost threads with something else for their platform.
AX_PKG(THREADS, [BOOST_THREAD], [], [])

AX_PKG(M, [], [-lm], [math.h])

AX_PKG_ONE_OF(LAPACK,
  APPLE_LAPACK,
    [AX_PKG_APPLE(LAPACK, [], [-framework Accelerate])],
  INTEL_LAPACK,
    [case $host_cpu in
       ia64)   AX_PKG(INTEL_LAPACK, [PTHREADS M], [-lmkl_lapack -lmkl -lmkl_i2p -lguide], [mkl_lapack.h], [], [64],    [dgelss_]);;
       x86_64) AX_PKG(INTEL_LAPACK, [PTHREADS M], [-lmkl_lapack -lmkl -lmkl_def -lguide], [mkl_lapack.h], [], [em64t], [dgelss_]);;
       *)      AX_PKG(INTEL_LAPACK, [PTHREADS M], [-lmkl_lapack -lmkl -lguide],           [mkl_lapack.h], [], [32],    [dgelss_]);;
     esac],
  CLAPACK,
    [AX_PKG(CLAPACK, [], [-lclapack -lblas -lf2c],     [], [], [], [dgelss_])],
  SLAPACK,
    [AX_PKG(SLAPACK, [], [-llapack -lblas],            [], [], [], [dgelss_])],
  FLAPACK,
    [AX_PKG(FLAPACK, [], [-llapack -lblas -lgfortran], [], [], [], [dgelss_])],
  STANDALONE_LAPACK_BLAS,
    [AX_PKG(STANDALONE_BLAS, [], [-lblas], [])
     AX_PKG(STANDALONE_LAPACK, [STANDALONE_BLAS], [-llapack], [], [], [], [dgelss_])
     AX_GROUP_PKG(STANDALONE_LAPACK_AND_BLAS, [STANDALONE_LAPACK STANDALONE_BLAS])],
  STANDALONE_FLAPACK_FBLAS,
    [AX_PKG(STANDALONE_F2C, [], [-lf2c], [])
     AX_PKG(STANDALONE_FBLAS, [STANDALONE_F2C], [-lblas], [])
     AX_PKG(STANDALONE_FLAPACK, [STANDALONE_FBLAS STANDALONE_F2C], [-llapack], [], [], [], [dgelss_])
     AX_GROUP_PKG(STANDALONE_FLAPACK_FBLAS, [STANDALONE_FLAPACK STANDALONE_FBLAS STANDALONE_F2C])])

AX_PKG(Z, [], [-lz], [zlib.h])

AX_PKG(PNG,  [Z], [-lpng], [png.h])
AX_PKG(JPEG, [], [-ljpeg], [stdio.h jpeglib.h]) # jpeglib.h requires stdio.h to have been included

AX_PKG(ILMBASE, [M Z PTHREADS], [-lImath -lHalf -lIex -lIlmThread], [Iex.h], [OpenEXR])
AX_PKG(OPENEXR, [ILMBASE], [-lIlmImf], [ImfOutputFile.h], [OpenEXR])

AX_PKG_ONE_OF(HDF,
    HDF_PLAIN,  [AX_PKG(HDF_PLAIN,  [M JPEG], [-ldf -lmfhdf], [mfhdf.h]       )],
    HDF_SUBDIR, [AX_PKG(HDF_SUBDIR, [M JPEG], [-ldf -lmfhdf], [mfhdf.h], [hdf])]
)

AX_PKG(TIFF, [Z JPEG], [-ltiff], [tiff.h])
if test x"$HAVE_PKG_TIFF" = "xyes"; then
  PKG_CHECK_FUNCTION(TIFF, [TIFFScanlineSize64], [BigTIFF])
  VW_HAS_BIGTIFF="$TIFF_HAS_BIGTIFF"
fi

AC_MSG_CHECKING([if VW has BigTIFF support])
if test x"$VW_HAS_BIGTIFF" = "xyes"; then
  has_bigtiff=1
else
  has_bigtiff=0
fi
AC_DEFINE_UNQUOTED(HAS_BIGTIFF, [$has_bigtiff], [Define to 1 if VW has BigTIFF support])
AC_MSG_RESULT([$VW_HAS_BIGTIFF])

AX_PKG(PROJ4,   [M],          [-lproj],          [proj_api.h])
AX_PKG(GEOTIFF, [PROJ4 JPEG], [-lgeotiff -ltiff] )
AX_PKG(GEOS,    [], [-lgeos -lgeos_c], [geos_c.h])
AX_PKG(GDAL,    [GEOTIFF GEOS],    [-lgdal],          [gdal.h])

# GDAL requires tiff support, but it can be internal or external.
PKG_CHECK_FUNCTION(GDAL, [TIFFScanlineSize], [External_LibTIFF])

if test x"$HAVE_PKG_GDAL" = "xyes"; then
  HAVE_PKG_TIFF="no:gdal has tiff support"
  if test x"$GDAL_HAS_EXTERNAL_LIBTIFF" = "xyes"; then
    PKG_CHECK_FUNCTION(GDAL, [TIFFScanlineSize64], [BigTIFF])
  else
    GDAL_HAS_BIGTIFF=yes
  fi
  VW_HAS_BIGTIFF="$GDAL_HAS_BIGTIFF"
fi

PKG_CHECK_MODULES([CAIROMM], [cairomm-1.0], [
  PKG_CAIROMM_CPPFLAGS="$PKG_CAIROMM_CPPFLAGS $CAIROMM_CFLAGS"
  ], [HAVE_PKG_CAIROMM="no:$CAIROMM_PKG_ERRORS"])
AX_PKG(CAIROMM, [], [$CAIROMM_LIBS], [cairomm/cairomm.h])

#Usage: AX_PKG(<name>, <dependencies>, <libraries>, <headers>[, <relative include path>, <relative lib path>, <required-functions>, <optional-deps>])
AX_PKG(LIBKML,   [BOOST_THREAD BOOST Z], [-lkmlengine -lkmldom -lkmlbase -lminizip], [kml/engine/kmz_file.h])
AX_PKG(OPENCV,     [], [-lopencv_core -lopencv_imgproc -lopencv_features2d -lopencv_xfeatures2d -lopencv_flann],  [ opencv2/core.hpp  opencv2/imgproc/imgproc.hpp opencv2/features2d.hpp opencv2/xfeatures2d.hpp opencv2/objdetect.hpp opencv2/calib3d.hpp])
AX_PKG(FLANN,      [], [-lflann_cpp],[flann/flann.hpp])

# These are here (instead of inside the PKG macro where they belong)
# for backwards compatability with older versions of automake.
AM_CONDITIONAL(HAVE_PKG_LAPACK,     [test "$HAVE_PKG_LAPACK"     = "yes"])
AM_CONDITIONAL(HAVE_PKG_PNG,        [test "$HAVE_PKG_PNG"        = "yes"])
AM_CONDITIONAL(HAVE_PKG_JPEG,       [test "$HAVE_PKG_JPEG"       = "yes"])
AM_CONDITIONAL(HAVE_PKG_TIFF,       [test "$HAVE_PKG_TIFF"       = "yes"])
AM_CONDITIONAL(HAVE_PKG_OPENEXR,    [test "$HAVE_PKG_OPENEXR"    = "yes"])
AM_CONDITIONAL(HAVE_PKG_HDF,        [test "$HAVE_PKG_HDF"        = "yes"])
AM_CONDITIONAL(HAVE_PKG_GDAL,       [test "$HAVE_PKG_GDAL"       = "yes"])
AM_CONDITIONAL(HAVE_PKG_CAIROMM,    [test "$HAVE_PKG_CAIROMM"    = "yes"])
AM_CONDITIONAL(HAVE_PKG_LIBKML,     [test "$HAVE_PKG_LIBKML"     = "yes"])
AM_CONDITIONAL(HAVE_PKG_OPENCV,     [test "$HAVE_PKG_OPENCV"     = "yes"])
AM_CONDITIONAL(HAVE_PKG_FLANN,      [test "$HAVE_PKG_FLANN"      = "yes"])
AM_CONDITIONAL(HAVE_PYTHON,         [test "$HAVE_PYTHON"         = "yes"])


##################################################
# module definitions
##################################################

AX_MODULE(CORE,   [src/vw/Core],   [libvwCore.la],   yes, [],      [BOOST BOOST_PROGRAM_OPTIONS THREADS M], [PTHREADS])
AX_MODULE(MATH,   [src/vw/Math],   [libvwMath.la],   yes, [CORE],  [BOOST_GRAPH],                           [LAPACK FLANN])
AX_MODULE(IMAGE,  [src/vw/Image],  [libvwImage.la],  yes, [MATH],  [OPENCV],                                [])
AX_MODULE(FILEIO, [src/vw/FileIO], [libvwFileIO.la], yes, [IMAGE], [BOOST_FILESYSTEM GDAL],                      [PNG JPEG TIFF Z OPENEXR HDF])
AX_MODULE(VW,     [src/vw],        [libvw.la],       yes, [],      [IMAGE MATH CORE FILEIO],                [])

if test "${MAKE_MODULE_VW}" != "yes"; then
  AC_MSG_ERROR([cannot build the Vision Workbench core modules!])
fi

# remove libvw.la from this. it's a convenience library that just happens to be
# installed, but if we link against it, we create multi-path problems
PKG_VW_LIBS_TEMP=
for i in $PKG_VW_LIBS; do
    case $i in
        *libvw\.*) ;;
        *) PKG_VW_LIBS_TEMP="$PKG_VW_LIBS_TEMP $i" ;;
    esac
done
PKG_VW_LIBS="$PKG_VW_LIBS_TEMP"

AX_MODULE(CAMERA,           [src/vw/Camera],           [libvwCamera.la],           yes, [VW],               [],                      [BOOST_IOSTREAMS])
AX_MODULE(INTERESTPOINT,    [src/vw/InterestPoint],    [libvwInterestPoint.la],    yes, [VW])
AX_MODULE(CARTOGRAPHY,      [src/vw/Cartography],      [libvwCartography.la],      yes, [VW],        [PROJ4 GDAL],            [])
AX_MODULE(MOSAIC,           [src/vw/Mosaic],           [libvwMosaic.la],           yes, [CARTOGRAPHY VW])
AX_MODULE(HDR,              [src/vw/HDR],              [libvwHDR.la],              yes, [CAMERA VW], [LAPACK])
AX_MODULE(STEREO,           [src/vw/Stereo],           [libvwStereo.la],           yes, [CAMERA VW])
AX_MODULE(GEOMETRY,         [src/vw/Geometry],         [libvwGeometry.la],         yes, [VW])
AX_MODULE(BUNDLEADJUSTMENT, [src/vw/BundleAdjustment], [libvwBundleAdjustment.la], yes, [CAMERA CARTOGRAPHY INTERESTPOINT STEREO VW])

AX_MODULE(TOOLS,   [src/vw/tools],  [],     yes, [VW], [BOOST_FILESYSTEM BOOST_PROGRAM_OPTIONS THREADS])
# Would like to delete this module but the googlenasa tools that populate the 
#  Google Earth Mars content actually depend on it!  We don't need to build it 
#  all the time, but don't delete it until the googlenasa stuff dies.
AX_MODULE_PYTHON(  [src/vw/Python], [no])

AX_MODULE(FLOOD_DETECT,   [src/vw/tools],  [],     yes, [VW], [BOOST_FILESYSTEM BOOST_PROGRAM_OPTIONS THREADS])


# These are here (instead of inside the MODULE macro where they belong)
# for backwards compatability with older versions of automake.
AM_CONDITIONAL(MAKE_MODULE_CORE,             [test "$MAKE_MODULE_CORE"   = "yes"])
AM_CONDITIONAL(MAKE_MODULE_MATH,             [test "$MAKE_MODULE_MATH"   = "yes"])
AM_CONDITIONAL(MAKE_MODULE_IMAGE,            [test "$MAKE_MODULE_IMAGE"  = "yes"])
AM_CONDITIONAL(MAKE_MODULE_FILEIO,           [test "$MAKE_MODULE_FILEIO" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_VW,               [test "$MAKE_MODULE_VW"     = "yes"])
AM_CONDITIONAL(MAKE_MODULE_CAMERA,           [test "$MAKE_MODULE_CAMERA" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_MOSAIC,           [test "$MAKE_MODULE_MOSAIC" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_INTERESTPOINT,    [test "$MAKE_MODULE_INTERESTPOINT" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_CARTOGRAPHY,      [test "$MAKE_MODULE_CARTOGRAPHY"   = "yes"])
AM_CONDITIONAL(MAKE_MODULE_HDR,              [test "$MAKE_MODULE_HDR"           = "yes"])
AM_CONDITIONAL(MAKE_MODULE_STEREO,           [test "$MAKE_MODULE_STEREO"        = "yes"])
AM_CONDITIONAL(MAKE_MODULE_GEOMETRY,         [test "$MAKE_MODULE_GEOMETRY"      = "yes"])
AM_CONDITIONAL(MAKE_MODULE_BUNDLEADJUSTMENT, [test "$MAKE_MODULE_BUNDLEADJUSTMENT" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_TOOLS,            [test "$MAKE_MODULE_TOOLS"  = "yes"])
AM_CONDITIONAL(MAKE_MODULE_PYTHON,           [test "$MAKE_MODULE_PYTHON" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_FLOOD_DETECT,     [test "$MAKE_MODULE_FLOOD_DETECT" = "yes"])

# Option to turn off building a bunch of apps we never use anymore
AM_CONDITIONAL(MAKE_APP_LEGACY, [test "$MAKE_APP_LEGACY" = "yes"])

##################################################
# final processing
##################################################

VW_CPPFLAGS="-I\${abs_top_srcdir}/src -I\${abs_top_builddir}/src $OTHER_CPPFLAGS"
VW_LDFLAGS="$OTHER_LDFLAGS"

# Do this last, to make sure it's not overridden
if test x"$ENABLE_STRICT_ALIASING" = x"no"; then
    CXXFLAGS="$CXXFLAGS -fno-strict-aliasing"
fi

AX_LOG([using CPPFLAGS=$CPPFLAGS])
AX_LOG([using VW_CPPFLAGS=$VW_CPPFLAGS])
AX_LOG([using CFLAGS=$CFLAGS])
AX_LOG([using CXXFLAGS=$CXXFLAGS])
AX_LOG([using LDFLAGS=$LDFLAGS])
AX_LOG([using VW_LDFLAGS=$VW_LDFLAGS])
AX_LOG([using VW_ENABLE_EXCEPTIONS=$ENABLE_EXCEPTIONS])

AC_SUBST(VW_CPPFLAGS)
AC_SUBST(VW_LDFLAGS)
AC_SUBST(LCOV_REMOTE_PATH)
AC_SUBST(DEHYDRA_GCC)
AC_SUBST(DEHYDRA_PLUGIN)

##################################################
# makefiles to process
##################################################

AC_CONFIG_FILES([                       \
  Makefile                              \
  src/Makefile                          \
  src/test/Makefile                     \
  src/vw/Makefile                       \
  src/vw/tools/Makefile                 \
  src/vw/Core/Makefile                  \
  src/vw/Core/tests/Makefile            \
  src/vw/Math/Makefile                  \
  src/vw/Math/tests/Makefile            \
  src/vw/Image/Makefile                 \
  src/vw/Image/tests/Makefile           \
  src/vw/FileIO/Makefile                \
  src/vw/FileIO/tests/Makefile          \
  src/vw/Camera/Makefile                \
  src/vw/Camera/tests/Makefile          \
  src/vw/Mosaic/Makefile                \
  src/vw/Mosaic/tests/Makefile          \
  src/vw/Cartography/Makefile           \
  src/vw/Cartography/tests/Makefile     \
  src/vw/HDR/Makefile                   \
  src/vw/HDR/tests/Makefile             \
  src/vw/Stereo/Makefile                \
  src/vw/Stereo/tests/Makefile          \
  src/vw/Geometry/Makefile              \
  src/vw/Geometry/tests/Makefile        \
  src/vw/InterestPoint/Makefile         \
  src/vw/InterestPoint/tests/Makefile   \
  src/vw/Python/vw/Makefile             \
  src/vw/Python/Makefile                \
  src/vw/BundleAdjustment/Makefile      \
  src/vw/BundleAdjustment/tests/Makefile \
])

AC_OUTPUT
