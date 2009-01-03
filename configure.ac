# __BEGIN_LICENSE__
#
# Copyright (C) 2006 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration
# (NASA).  All Rights Reserved.
#
# Copyright 2006 Carnegie Mellon University. All rights reserved.
#
# This software is distributed under the NASA Open Source Agreement
# (NOSA), version 1.3.  The NOSA has been approved by the Open Source
# Initiative.  See the file COPYING at the top of the distribution
# directory tree for the complete NOSA document.
#
# THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
# KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
# LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
# SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
# THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
# DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
#
# __END_LICENSE__

##################################################
# A few starting bits
##################################################

dnl TODO BEFORE RELEASE: update package version
AC_INIT([NASA Vision Workbench], [2.0.alpha5], [vision-workbench@lists.nasa.gov], [VisionWorkbench])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE
AX_CONFIG_HEADER_PREFIX([src/vw/config.h], [VW_])

dnl TODO BEFORE RELEASE: update libtool version
dnl if abi changes in any way:
dnl     current++
dnl     revision=0
dnl     if new abi is backward-compatible:
dnl         age += 1
dnl     else:
dnl         age = 0
dnl else:
dnl     revision++
dnl (see http://sources.redhat.com/autobook/autobook/autobook_91.html)
dnl
dnl LIBTOOL_VERSION=current:revision:age
LIBTOOL_VERSION=1:0:0
AC_SUBST(LIBTOOL_VERSION)

# This should be first, so everything else can use it
AX_ARG_ENABLE(verbose, no, [none], [verbose configuration output])

# load the config.options
if test -f config.options; then
  if test "yes" = "$ENABLE_VERBOSE" ; then
    AC_MSG_NOTICE([sourcing config.options])
    cat config.options
  fi
  source config.options
fi



##################################################
# Find the necessary programs
##################################################

# Make sure these are set so AC_PROG_CC doesn't initialize them with
# its defaults, since we take care of setting debug and optimization
# levels ourselves.
CFLAGS="$CFLAGS"
CXXFLAGS="$CXXFLAGS"

# Find the tools
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AX_PROG_AR
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET
AC_LANG(C++)

# Does the compiler support introspection?
AX_CHECK_INTROSPECTION

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

# We set up the include paths ourselves, so we ask autoconf to get out of the way
DEFAULT_INCLUDES=
AC_SUBST(DEFAULT_INCLUDES)



##################################################
# Compilation options
##################################################

# For the moment this is just a manual switch, defaulting to off.
# Someday this should be sensibly autodetected.
AX_ARG_ENABLE(sse,           no, [am-yes cpp-bool], [enable SSE optimizations in some places (development)])
AX_ARG_ENABLE(exceptions,   yes, [am-yes cpp-bool], [enable the C++ exception mechanism])
AX_ARG_ENABLE(bounds-check,  no, [am-yes cpp-bool], [enable image bounds checking (SLOW!)])
AX_ARG_ENABLE(test-coverage, no, [none],            [turn on code coverage analysis tools])
AX_ARG_ENABLE(debug,         no, [none],            [generate debugging symbols])
AX_ARG_ENABLE(optimize,       3, [none],            [compiler optimization level])
AX_ARG_ENABLE(multi-arch,    [], [none],            [build multi-arch (universal) binaries])
AX_ARG_ENABLE(lib64,       auto, [none],            [force /lib64 instead of /lib])



##################################################
# Handle options
##################################################

# Pass apple gcc options to build a universal binary
for arch in $ENABLE_MULTI_ARCH; do
    AX_CFLAGS="$AX_CFLAGS -arch $arch"
done

# For 64-bit machines, we'll generally want to autodetect if lib64 exists,
# and use it. Sometimes the user will want to force use of lib64 or lib
# instead of the autodetection. Thus we have this command line option.
if test x"$USE_LIB64" = "xauto"; then
  if test -d "/lib64" -o -d "/usr/lib64"; then
    USE_LIB64="yes"
  else
    USE_LIB64="no"
  fi
fi

# Sometimes we have /foo/lib64 and /foo/lib confusion on 64-bit machines,
# so we'll use possibly both if one doesn't appear for a certain
# library path.
if test x"$USE_LIB64" = "xyes"; then
  AX_LIBDIR="lib64"
  AX_OTHER_LIBDIR="lib"
else
  AX_LIBDIR="lib"
  AX_OTHER_LIBDIR="lib64"
fi

# Check for --enable-test-coverate and add appropriate flags for gcc
if test "$ENABLE_TEST_COVERAGE" = "yes"; then
  AX_TRY_CPPFLAGS([-fprofile-arcs],  [ENABLE_TEST_COVERAGE=yes], [ENABLE_TEST_COVERAGE=no])
  AX_TRY_CPPFLAGS([-ftest-coverage], [],                         [ENABLE_TEST_COVERAGE=no])

  if test "$ENABLE_TEST_COVERAGE" = "yes"; then
    echo "*** Enabling test coverage and debugging, and disabling optimization ***"
    CXXFLAGS="$CXXFLAGS -fprofile-arcs -ftest-coverage"
    LDFLAGS="-lgcov"
    ENABLE_OPTIMIZE=no
    ENABLE_DEBUG=yes
  fi
fi

# Files test-coverage and profile-arcs generate
AC_SUBST([MOSTLYCLEANFILES], "*.bb *.bbg *.da *.gcov *.gcda *.gcno")

# These are good if they're supported
AX_TRY_CPPFLAGS([-Wl,-no-undefined], [LDFLAGS="$LDFLAGS -Wl,-no-undefined"])
AX_TRY_CPPFLAGS([-Wl,-as-needed],    [LDFLAGS="$LDFLAGS -Wl,-as-needed"])

# Debugging
case "$ENABLE_DEBUG" in
    yes|1) AX_CFLAGS="$AX_CFLAGS -g -DDEBUG" ;;
    2)     AX_CFLAGS="$AX_CFLAGS -g -DDEBUG -D_GLIBCXX_DEBUG" ;;
    no)    AX_CFLAGS="$AX_CFLAGS -DNDEBUG" ;;
    *)     AC_MSG_ERROR([Unknown debug option: "$ENABLE_DEBUG"]) ;;
esac

# Optimization
case "$ENABLE_OPTIMIZE" in
    yes)     AX_CFLAGS="$AX_CFLAGS -O3" ;;
    3|2|1)   AX_CFLAGS="$AX_CFLAGS -O$ENABLE_OPTIMIZE" ;;
    coreduo) AX_CFLAGS="$AX_CFLAGS -O4 -march=prescott -mtune=prescott -funroll-loops -msse -msse2 -msse3 -mfpmath=sse -ftree-vectorize" ;;
    sse3)    AX_CFLAGS="$AX_CFLAGS -O4 -funroll-loops -msse -msse2 -msse3 -mfpmath=sse -ftree-vectorize" ;;
    no|0)    AC_MSG_WARN([*** The Vision Workbench may not work properly with optimization disabled! ***])
             AX_CFLAGS="$AX_CFLAGS -O0" ;;
    *)       AC_MSG_ERROR([Unknown optimize option: "$ENABLE_OPTIMIZE"]) ;;
esac

CFLAGS="$CFLAGS $AX_CFLAGS"
CXXFLAGS="$CXXFLAGS $AX_CFLAGS"



##################################################
# installation options
##################################################

if test ${prefix} = NONE; then
  if test ! -z ${PREFIX} ; then
    prefix="${PREFIX}"
    AC_MSG_NOTICE([using installation prefix ${prefix}])
  fi
fi



##################################################
# distribution options
##################################################

AX_ARG_WITH(dist-license,                       [], [mk am-set], [special distribution license file to be included as the COPYING file in the distribution])
AX_ARG_WITH(dist-license-summary,               [], [mk am-set], [special distribution license summary file to be included in the headers of source files in the distribution])
AX_ARG_WITH(dist-config-options-default,        [], [mk am-set], [special distribution config.options.default file])
AX_ARG_WITH(dist-remove-modules, [Geometry Python gui], [mk am-set], [remove certain modules from the distribution])



##################################################
# handle distribution options
##################################################

if test ! -z "$DIST_LICENSE" ; then
  AC_MSG_NOTICE([using distribution license file $DIST_LICENSE])
fi

if test ! -z "$DIST_LICENSE_SUMMARY" ; then
  AC_MSG_NOTICE([using distribution license summary file $DIST_LICENSE_SUMMARY])
fi

if test ! -z "$DIST_CONFIG_OPTIONS_DEFAULT" ; then
  if test "yes" = "$DIST_CONFIG_OPTIONS_DEFAULT" ; then
    DIST_CONFIG_OPTIONS_DEFAULT="config/config.options.default"
  fi
  AC_MSG_NOTICE([using distribution config.options.default file: $DIST_CONFIG_OPTIONS_DEFAULT])
fi



##################################################
# package checks
##################################################

PKG_PATHS_FROM_OPT="$PKG_PATHS"
AX_ARG_WITH(pkg-paths, [], [none], [additional search path(s) for packages])
PKG_PATHS="$PKG_PATHS $PKG_PATHS_FROM_OPT ${HOME} ${HOME}/local /sw /opt /opt/local /usr/local /usr/X11R6 /usr"

if test "yes" = "$ENABLE_VERBOSE"; then
  AC_MSG_NOTICE([using PKG_PATHS=$PKG_PATHS])
fi

AX_PKG_PTHREADS

AX_PKG_BOOST
AX_PKG_BOOST_CHECK_VERSION([VW_])

AC_MSG_NOTICE([VW will be built using Boost version $BOOST_VERSION])
AX_PKG_BOOST_LIB(PROGRAM_OPTIONS, [-lboost_program_options], [boost/program_options.hpp])

AX_PKG_ONE_OF(BOOST_FILESYSTEM,
    BOOST_FILESYSTEM_PRE_1_35,  [AX_PKG_BOOST_LIB(FILESYSTEM_PRE_1_35,  [-lboost_filesystem], [boost/filesystem/path.hpp])],
    BOOST_FILESYSTEM_POST_1_35, [AX_PKG_BOOST_LIB(FILESYSTEM_POST_1_35, [-lboost_filesystem -lboost_system], [boost/filesystem/path.hpp])])

AX_PKG_BOOST_LIB(THREAD, [-lboost_thread], [boost/thread/thread.hpp])

# This provides an easy way for users to override the threads options
# if they are replacing Boost threads with something else for their platform.
AX_PKG(THREADS, [BOOST_THREAD], [], [])

AX_PKG_LAPACK

AX_PKG(QHULL, [], [-lqhull], [stdio.h stdlib.h string.h ctype.h math.h qhull/qhull.h])

AX_PKG(GMP, [], [-lgmpxx -lgmp], [gmpxx.h gmp.h])
AX_PKG(MPFR, [], [-lmpfr], [mpfr.h])
AX_PKG(PPL, [GMP], [-lppl], [ppl.hh])
AX_PKG(APRON, [GMP MPFR], [-lpolkag -lapron -litvmpq -litvdbl], [ap_global0.h ap_global1.h pk.h])

AX_PKG(M, [], [-lm], [math.h])
AX_PKG(Z, [], [-lz], [zlib.h])

AX_PKG(PNG, [], [-lpng], [png.h])
AX_PKG(JPEG, [], [-ljpeg], [stdio.h jpeglib.h]) # jpeglib.h requires stdio.h to have been included

AX_PKG(TIFF, [], [-ltiff], [tiff.h])
tiff_has_bigtiff=0
if test x"$HAVE_PKG_TIFF" != x"no"; then
  AX_CHECK_FUNCTIONS(TIFFScanlineSize64, [$PKG_TIFF_LIBS], [tiff_has_bigtiff=1], [], [for BigTIFF in TIFF])
fi
AC_DEFINE_UNQUOTED([TIFF_HAS_BIGTIFF], [$tiff_has_bigtiff], [Define to 1 if libtiff has BigTIFF support])

AX_PKG(ILMBASE, [Z PTHREADS], [-lImath -lHalf -lIex], [Iex.h], [OpenEXR])
AX_PKG(OPENEXR, [ILMBASE], [-lIlmImf], [ImfOutputFile.h], [OpenEXR])

AX_PKG_ONE_OF(HDF,
    HDF_PLAIN,  [AX_PKG(HDF_PLAIN,  [], [-ldf -lmfhdf], [mfhdf.h]       )],
    HDF_SUBDIR, [AX_PKG(HDF_SUBDIR, [], [-ldf -lmfhdf], [mfhdf.h], [hdf])]
)

AX_PKG(X11, [], [-lXext -lX11 -lXmu -lXi -lXt], [])
AX_PKG_GL
#AX_PKG(GL, [X11], [-lGL -lGLU], [GL/gl.h GL/glu.h])
#AX_PKG(GLUT, [X11 GL], [-lglut], [GL/glut.h])
AX_PKG(GLEW, [GL M], [-lGLEW], [GL/glew.h])
AX_PKG(CG, [], [-lCg -lCgGL], [Cg/cg.h])

AX_PKG(GDAL, [], [-lgdal], [gdal.h])
gdal_has_bigtiff=0
if test x"$HAVE_PKG_GDAL" != x"no"; then
  AX_CHECK_FUNCTIONS(TIFFScanlineSize64, [$PKG_GDAL_LIBS], [gdal_has_bigtiff=1], [], [for BigTIFF in GDAL])
fi
AC_DEFINE_UNQUOTED([GDAL_HAS_BIGTIFF], [$gdal_has_bigtiff], [Define to 1 if libgdal has BigTIFF support])

AX_PKG(PROJ4, [], [-lproj], [projects.h])

# Qt with AutoTroll.
AT_WITH_QT([+gui +sql +opengl], [], [], [], [HAVE_PKG_QT_INCLUDE=no])

if test x"$HAVE_PKG_QT_INCLUDE" != x"no"; then
    AT_REQUIRE_QT_VERSION([4.0], [], [HAVE_PKG_QT_INCLUDE=no])
fi

# We seed QT_INCLUDE's cppflags with the ones qmake figured out
ADD_QT_INCLUDE_CPPFLAGS="$ADD_QT_INCLUDE_CPPFLAGS $QT_CPPFLAGS"
AX_PKG(QT_INCLUDE, [], [], [QVector])

if test x"$host_vendor" = "xapple"; then
    AX_PKG(QT,         [QT_INCLUDE], [$QT_LIBS], [])
else
    AX_PKG(QT_LIBS,              [], [$QT_LIBS], [])
    AX_PKG(QT_OPENGL,  [QT_INCLUDE], [],         [QtOpenGL])
    AX_PKG(QT_GUI,     [QT_INCLUDE], [],         [QLabel])
    AX_PKG(QT_SQL,     [QT_INCLUDE], [],         [QtSql])
    AX_GROUP_PKG(QT,   [QT_INCLUDE QT_GUI QT_OPENGL QT_SQL QT_LIBS])
fi

# These are here (instead of inside the PKG macro where they belong)
# for backwards compatability with older versions of automake.
AM_CONDITIONAL(HAVE_PKG_LAPACK, [test "$HAVE_PKG_LAPACK" = "yes"])
AM_CONDITIONAL(HAVE_PKG_QHULL, [test "$HAVE_PKG_QHULL" = "yes"])
AM_CONDITIONAL(HAVE_PKG_GMP, [test "$HAVE_PKG_GMP" = "yes"])
AM_CONDITIONAL(HAVE_PKG_MPFR, [test "$HAVE_PKG_MPFR" = "yes"])
AM_CONDITIONAL(HAVE_PKG_PPL, [test "$HAVE_PKG_PPL" = "yes"])
AM_CONDITIONAL(HAVE_PKG_APRON, [test "$HAVE_PKG_APRON" = "yes"])
AM_CONDITIONAL(HAVE_PKG_PNG, [test "$HAVE_PKG_PNG" = "yes"])
AM_CONDITIONAL(HAVE_PKG_JPEG, [test "$HAVE_PKG_JPEG" = "yes"])
AM_CONDITIONAL(HAVE_PKG_TIFF, [test "$HAVE_PKG_TIFF" = "yes"])
AM_CONDITIONAL(HAVE_PKG_OPENEXR, [test "$HAVE_PKG_OPENEXR" = "yes"])
AM_CONDITIONAL(HAVE_PKG_HDF, [test "$HAVE_PKG_HDF" = "yes"])
AM_CONDITIONAL(HAVE_PKG_GDAL, [test "$HAVE_PKG_GDAL" = "yes"])


##################################################
# module definitions
##################################################

AX_MODULE(CORE, [src/vw/Core], [libvwCore.la], yes, [], [BOOST THREADS M], [PTHREADS])
AX_MODULE(MATH, [src/vw/Math], [libvwMath.la], yes, [CORE], [], [LAPACK])
AX_MODULE(IMAGE, [src/vw/Image], [libvwImage.la], yes, [MATH])
AX_MODULE(FILEIO, [src/vw/FileIO], [libvwFileIO.la], yes, [IMAGE], [], [Z PNG JPEG TIFF OPENEXR HDF GDAL])
AX_MODULE(VW, [src/vw], [libvw.la], yes, [], [CORE MATH IMAGE], [FILEIO])
if test "${MAKE_MODULE_VW}" != "yes"; then
  AC_MSG_ERROR([cannot build the Vision Workbench core modules!])
fi

AX_MODULE(CAMERA, [src/vw/Camera], [libvwCamera.la], yes, [VW])
AX_MODULE(MOSAIC, [src/vw/Mosaic], [libvwMosaic.la], yes, [VW])
AX_MODULE(INTERESTPOINT, [src/vw/InterestPoint], [libvwInterestPoint.la], yes, [VW])
AX_MODULE(CARTOGRAPHY, [src/vw/Cartography], [libvwCartography.la], yes, [VW], [PROJ4], [GDAL])
AX_MODULE(HDR, [src/vw/HDR], [libvwHDR.la], yes, [VW], [CAMERA LAPACK])
AX_MODULE(STEREO, [src/vw/Stereo], [libvwStereo.la], yes, [VW], [])
if test $host_vendor = apple; then
  AX_MODULE(GPU, [src/vw/GPU], [libvwGPU.la], no, [VW], [GL], [CG])
else
  AX_MODULE(GPU, [src/vw/GPU], [libvwGPU.la], no, [VW], [GL GLEW], [CG])
fi
AX_MODULE(GEOMETRY, [src/vw/Geometry], [libvwGeometry.la], yes, [VW], [], [GMP PPL APRON QHULL])

AX_MODULE(TOOLS, [src/vw/tools], [], yes, [VW], [BOOST_FILESYSTEM BOOST_PROGRAM_OPTIONS THREADS])
AX_MODULE(GUI, [src/vw/gui], [], yes, [VW], [QT BOOST_PROGRAM_OPTIONS GL])

AX_MODULE_PYTHON([src/vw/Python], [no])

# These are here (instead of inside the MODULE macro where they belong)
# for backwards compatability with older versions of automake.
AM_CONDITIONAL(MAKE_MODULE_CORE, [test "$MAKE_MODULE_CORE" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_MATH, [test "$MAKE_MODULE_MATH" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_IMAGE, [test "$MAKE_MODULE_IMAGE" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_FILEIO, [test "$MAKE_MODULE_FILEIO" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_VW, [test "$MAKE_MODULE_VW" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_CAMERA, [test "$MAKE_MODULE_CAMERA" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_MOSAIC, [test "$MAKE_MODULE_MOSAIC" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_INTERESTPOINT, [test "$MAKE_MODULE_INTERESTPOINT" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_CARTOGRAPHY, [test "$MAKE_MODULE_CARTOGRAPHY" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_HDR, [test "$MAKE_MODULE_HDR" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_STEREO, [test "$MAKE_MODULE_STEREO" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_GPU, [test "$MAKE_MODULE_GPU" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_GEOMETRY, [test "$MAKE_MODULE_GEOMETRY" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_TOOLS, [test "$MAKE_MODULE_TOOLS" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_GUI, [test "$MAKE_MODULE_GUI" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_PYTHON, [test "$MAKE_MODULE_PYTHON" = "yes"])

##################################################
# library options
##################################################

AX_ARG_WITH(num-threads, [4], [cpp-int], [set the default number of processing threads for multi-threaded operations])
AC_MSG_NOTICE([VW will use $NUM_THREADS THREADS by default])

##################################################
# final processing
##################################################

VW_CPPFLAGS="-I\${top_srcdir}/src $OTHER_CPPFLAGS"
VW_LDFLAGS="$OTHER_LDFLAGS"

if test "yes" = "$ENABLE_VERBOSE" ; then
  AC_MSG_NOTICE([using CPPFLAGS=$CPPFLAGS])
  AC_MSG_NOTICE([using VW_CPPFLAGS=$VW_CPPFLAGS])
  AC_MSG_NOTICE([using CFLAGS=$CFLAGS])
  AC_MSG_NOTICE([using CXXFLAGS=$CXXFLAGS])
  AC_MSG_NOTICE([using LDFLAGS=$LDFLAGS])
  AC_MSG_NOTICE([using VW_LDFLAGS=$VW_LDFLAGS])
  AC_MSG_NOTICE([using VW_ENABLE_EXCEPTIONS=$ENABLE_EXCEPTIONS])
fi

AC_SUBST(VW_CPPFLAGS)
AC_SUBST(VW_LDFLAGS)


##################################################
# makefiles to process
##################################################

AC_CONFIG_FILES([                       \
  Makefile                              \
  src/Makefile                          \
  src/vw/Makefile                       \
  src/vw/tests/Makefile                 \
  src/vw/tools/Makefile                 \
  src/vw/gui/Makefile                   \
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
  src/vw/GPU/Makefile                   \
  src/vw/GPU/tests/Makefile             \
  src/vw/InterestPoint/Makefile         \
  src/vw/InterestPoint/tests/Makefile   \
  src/vw/Geometry/Makefile              \
  src/vw/Geometry/tests/Makefile        \
  src/vw/Python/vw/Makefile             \
  src/vw/Python/Makefile                \
])

AC_OUTPUT
