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

AC_INIT([NASA Vision Workbench], [2.0.alpha5-pre], [vision-workbench@lists.nasa.gov], [VisionWorkbench])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE
AX_CONFIG_HEADER_PREFIX([src/vw/config.h], [VW_])

LIBTOOL_CURRENT=1
LIBTOOL_REVISION=0
LIBTOOL_AGE=0
LIBTOOL_VERSION=$LIBTOOL_CURRENT:$LIBTOOL_REVISION:$LIBTOOL_AGE
AC_SUBST(LIBTOOL_VERSION)


##################################################
#	configuration options
##################################################

# This should be first, so everything else can use it
AC_ARG_ENABLE(verbose,
  AC_HELP_STRING([--enable-verbose], [verbose configuration output @<:@no@:>@]),           
  [ ENABLE_VERBOSE=$enableval ],
  [ if test -z $"ENABLE_VERBOSE"; then ENABLE_VERBOSE=no; fi ]
)

if test -f config.options; then
  if test "yes" = "$ENABLE_VERBOSE" ; then
    AC_MSG_NOTICE([sourcing config.options])
    cat config.options
  fi
  source config.options
fi

AC_ARG_WITH(paths,
  AC_HELP_STRING([--with-paths=PATHS], [additional search path(s) for packages]), 
  [ PKG_PATHS="$withval $PKG_PATHS" ]
)
PKG_PATHS="${PKG_PATHS} ${HOME} ${HOME}/local /sw /opt /opt/local /usr/local /usr/X11R6 /usr"


##################################################
#	compilation options	    
##################################################

# Debugging
AC_ARG_ENABLE(debug,
  AC_HELP_STRING([--enable-debug@<:@=LEVEL@:>@], [generate debugging symbols @<:@no@:>@]),
  [ ENABLE_DEBUG=$enableval ], 
  [ if test -z "$ENABLE_DEBUG"; then ENABLE_DEBUG=no; fi ]
)

if test "yes" = "$ENABLE_DEBUG" || test "1" = "$ENABLE_DEBUG"; then
   AX_CFLAGS="$AX_CFLAGS -g -DDEBUG"
elif test "2" = "$ENABLE_DEBUG"; then
   AX_CFLAGS="$AX_CFLAGS -g -DDEBUG -D_GLIBCXX_DEBUG"
else
   AX_CFLAGS="$AX_CFLAGS -DNDEBUG"
fi

# Target Architecture
AC_ARG_ENABLE(architecture,
  AC_HELP_STRING([--enable-architecture@<:@=ARCH@:>@], [select a target architecture @<:@x86_64@:>@]),
  [ TARGET_ARCHITECTURE=$enableval ],
  [ if test -z "$TARGET_ARCHITECTURE"; then echo; fi ]
)

if test ! -z "$TARGET_ARCHITECTURE"; then
  AC_MSG_NOTICE([Setting target architecture to $TARGET_ARCHITECTURE])
  AX_CFLAGS="$AX_CFLAGS -arch $TARGET_ARCHITECTURE"
else 
  AC_MSG_NOTICE([Setting target architecture to default])
fi 

# For 64-bit machines, we'll generally want to autodetect if lib64 exists,
# and use it. Sometimes the user will want to force use of lib64 or lib
# instead of the autodetection. Thus we have this command line option.
AC_ARG_ENABLE(lib64,
  AC_HELP_STRING([--enable-lib64], [force configure to use the lib64 directory (yes) or the lib directory (no) to find libraries @<:@autodetect@:>@]),
  [ USE_LIB64=$enableval ],
  [ if test -z "$USE_LIB64"; then USE_LIB64="autodetect"; fi ]
)
if test x"$USE_LIB64" = "xautodetect"; then
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

# Optimization
AC_ARG_ENABLE(optimize,
  AC_HELP_STRING([--enable-optimize@<:@=LEVEL@:>@], [compiler optimization level @<:@yes@:>@]), 
  [ ENABLE_OPTIMIZE=$enableval ],
  [ if test -z "$ENABLE_OPTIMIZE"; then ENABLE_OPTIMIZE=yes; fi ]
)

if test "yes" = "$ENABLE_OPTIMIZE" || test "3" = "$ENABLE_OPTIMIZE"; then
   AX_CFLAGS="$AX_CFLAGS -O3"
elif test "2" = "$ENABLE_OPTIMIZE"; then
   AX_CFLAGS="$AX_CFLAGS -O2"
elif test "1" = "$ENABLE_OPTIMIZE"; then
   AX_CFLAGS="$AX_CFLAGS -O1"
elif test "coreduo" = "$ENABLE_OPTIMIZE"; then
   AX_CFLAGS="$AX_CFLAGS -O4 -march=prescott -mtune=prescott -funroll-loops -msse -msse2 -msse3 -mfpmath=sse -ftree-vectorize"
elif test "sse3" = "$ENABLE_OPTIMIZE"; then
   AX_CFLAGS="$AX_CFLAGS -O4 -funroll-loops -msse -msse2 -msse3 -mfpmath=sse -ftree-vectorize"
else
   echo 'WARNING: *** The Vision Workbench may not work properly with optimization disabled! ***'
   AX_CFLAGS="$AX_CFLAGS -O0"
fi

CFLAGS="$CFLAGS -I\${top_srcdir}/src $AX_CFLAGS"
CXXFLAGS="$CXXFLAGS  -I\${top_srcdir}/src $AX_CFLAGS"
if test "x$host_vendor" != "xapple"; then
  LDFLAGS="$LDFLAGS -Wl,-no-undefined -Wl,-as-needed"
fi

AC_ARG_ENABLE(exceptions,
  AC_HELP_STRING([--enable-exceptions@<:@=yes@:>@], [enable the C++ exception mechanism]), 
  [ ENABLE_EXCEPTIONS=$enableval ],
  [ if test -z "$ENABLE_EXCEPTIONS"; then ENABLE_EXCEPTIONS=yes; fi ]
)

if test "$ENABLE_EXCEPTIONS" = "yes"; then
   no_exceptions=0
else
   no_exceptions=1
fi

AM_CONDITIONAL(ENABLE_EXCEPTIONS, [test "$ENABLE_EXCEPTIONS" = "yes"])

AC_DEFINE_UNQUOTED(NO_EXCEPTIONS,
                   [$no_exceptions],
                   [Define to 1 to disable C++ exceptions.])

# For the moment this is just a manual switch, defaulting to off.
# Someday this should be sensibly autodetected.
AC_ARG_ENABLE(sse,
  AC_HELP_STRING([--enable-sse@<:@=no@:>@], [enable SSE optimizations in some places (development)]), 
  [ ENABLE_SSE=$enableval ],
  [ if test -z "$ENABLE_SSE"; then ENABLE_SSE=no; fi ]
)

if test "$ENABLE_SSE" = "yes"; then
   enable_sse=1
else
   enable_sse=0
fi

AC_DEFINE_UNQUOTED(ENABLE_SSE,
                   [$enable_sse],
                   [Define to 1 to enable SSE optimizations.])

AC_ARG_ENABLE(bounds-check,
  AC_HELP_STRING([--enable-bounds_check@<:@=no@:>@], [enable image bounds checking (SLOW!).]),
  [ ENABLE_BOUNDS_CHECK=$enableval ],
  [ if test -z "$ENABLE_BOUNDS_CHECK"; then ENABLE_BOUNDS_CHECK=no; fi ]
)

if test "$ENABLE_BOUNDS_CHECK" = "yes"; then
   enable_bounds_check=1
else
   enable_bounds_check=0
fi

AC_DEFINE_UNQUOTED(IMAGE_BOUNDS_CHECK,
                   [$enable_bounds_check],
                   [Define to 1 to enable image bounds checking (SLOW!).])

# We set up the include paths ourselves, so 
# we ask autoconf to get out of the way.
DEFAULT_INCLUDES=
AC_SUBST(DEFAULT_INCLUDES)


##################################################
#	installation options
##################################################

if test ${prefix} = NONE; then
  if test ! -z ${PREFIX} ; then
    prefix=${PREFIX}
    AC_MSG_NOTICE([using installation prefix ${prefix}])
  fi
fi


##################################################
#	distribution options	    
##################################################

AC_ARG_WITH(dist-license,
  AC_HELP_STRING([--with-dist-license=LICENSE], [special distribution license file to be included as the COPYING file in the distribution]),
  [DIST_LICENSE=$withval]
)
AC_SUBST(DIST_LICENSE)
AM_CONDITIONAL(DIST_LICENSE, [test ! -z "$DIST_LICENSE"])

if test ! -z "$DIST_LICENSE" ; then
  AC_MSG_NOTICE([using distribution license file $DIST_LICENSE])
fi

AC_ARG_WITH(dist-license-summary,
  AC_HELP_STRING([--with-dist-license-summary=SUMMARY], [special distribution license summary file to be included in the headers of source files in the distribution]),
  [DIST_LICENSE_SUMMARY=$withval]
)
AC_SUBST(DIST_LICENSE_SUMMARY)
AM_CONDITIONAL(DIST_LICENSE_SUMMARY, [test ! -z "$DIST_LICENSE_SUMMARY"])

if test ! -z "$DIST_LICENSE_SUMMARY" ; then
  AC_MSG_NOTICE([using distribution license summary file $DIST_LICENSE_SUMMARY])
fi

AC_ARG_WITH(dist-config-options-default,
  AC_HELP_STRING([--with-dist-config-options-default@<:@=FILE@:>@], [special distribution config.options.default file]),
  [DIST_CONFIG_OPTIONS_DEFAULT=$withval]
)
AC_SUBST(DIST_CONFIG_OPTIONS_DEFAULT)
AM_CONDITIONAL(DIST_CONFIG_OPTIONS_DEFAULT, [test ! -z "$DIST_CONFIG_OPTIONS_DEFAULT"])

if test ! -z "$DIST_CONFIG_OPTIONS_DEFAULT" ; then
  if test "yes" = "$DIST_CONFIG_OPTIONS_DEFAULT" ; then
    DIST_CONFIG_OPTIONS_DEFAULT="config/config.options.default"
  fi
  AC_MSG_NOTICE([using distribution config.options.default file: $DIST_CONFIG_OPTIONS_DEFAULT])
fi

AC_ARG_WITH(dist-remove-modules,
  AC_HELP_STRING([--with-dist-remove-modules="MODULE @<:@MODULE...@:>@"], [remove certain modules from the distribution]),
  [DIST_REMOVE_MODULES=$withval], [DIST_REMOVE_MODULES="Geometry Python"]
)
AC_SUBST(DIST_REMOVE_MODULES)
AM_CONDITIONAL(DIST_REMOVE_MODULES, [test ! -z "$DIST_REMOVE_MODULES"])


##################################################
#	package checks	    
##################################################

# Make sure these are set so AC_PROG_CC doesn't initialize them with
# its defaults, since we take care of setting debug and optimization
# levels ourselves.
CFLAGS="$CFLAGS"
CXXFLAGS="$CXXFLAGS"

AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AX_PROG_AR
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET

AM_PATH_PYTHON([2.4], [HAVE_PYTHON=yes], [HAVE_PYTHON=no])
AC_PYTHON_DEVEL([HAVE_PYTHON=yes], [HAVE_PYTHON=no])
AC_PROG_SWIG([1.3.29], [HAVE_SWIG=yes], [HAVE_SWIG=no])

if test "$HAVE_SWIG" = "yes"; then
    SWIG_ENABLE_CXX
    if test "$HAVE_PYTHON" = "yes"; then
        SWIG_PYTHON
    fi
fi

AC_SUBST(HAVE_PYTHON)
AC_SUBST(HAVE_SWIG)

AC_LANG(C++)

AX_CHECK_INTROSPECTION

if test "yes" = "$ENABLE_VERBOSE"; then
  AC_MSG_NOTICE([using PKG_PATHS=$PKG_PATHS])
fi

AX_PKG_PTHREADS

AX_PKG_BOOST
AX_PKG_BOOST_CHECK_VERSION([VW_])

AC_MSG_NOTICE([VW will be built using Boost version $BOOST_VERSION])
AX_PKG_BOOST_LIB(PROGRAM_OPTIONS, [-lboost_program_options], [boost/program_options.hpp])
AX_PKG_BOOST_LIB(FILESYSTEM_PRE_1_35, [-lboost_filesystem], [boost/filesystem/path.hpp])
if test "$HAVE_PKG_BOOST_FILESYSTEM_PRE_1_35" = "no"; then
  echo "Checking for BOOST_FILESYSTEM 1.35 or later..."
  AX_PKG_BOOST_LIB(FILESYSTEM, [-lboost_filesystem -lboost_system], [boost/filesystem/path.hpp])
else 
  AX_PKG(BOOST_FILESYSTEM, [BOOST_FILESYSTEM_PRE_1_35], [], [])
fi
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
AX_PKG(ILMBASE, [Z PTHREADS], [-lImath -lHalf -lIex], [Iex.h], [OpenEXR])
AX_PKG(OPENEXR, [ILMBASE], [-lIlmImf], [ImfOutputFile.h], [OpenEXR])

AX_PKG(HDF_PLAIN, [], [-ldf -lmfhdf], [mfhdf.h])
if test "$HAVE_PKG_HDF_PLAIN" = "no"; then
    AX_PKG(HDF, [], [-ldf -lmfhdf], [mfhdf.h], [hdf])
else
    AX_PKG(HDF, [HDF_PLAIN], [], [])
fi

AX_PKG(X11, [], [-lXext -lX11 -lXmu -lXi -lXt], [])
AX_PKG_GL
#AX_PKG(GL, [X11], [-lGL -lGLU], [GL/gl.h GL/glu.h])
#AX_PKG(GLUT, [X11 GL], [-lglut], [GL/glut.h])
AX_PKG(GLEW, [GL M], [-lGLEW], [GL/glew.h])
AX_PKG(CG, [], [-lCg -lCgGL], [Cg/cg.h])

AX_PKG(GDAL, [], [-lgdal], [gdal.h])
AX_PKG(PROJ4, [], [-lproj], [projects.h])

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
#	module definitions
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

AX_MODULE(TOOLS, [src/vw/tools], [], yes, [], [])

AX_MODULE_PYTHON([no])

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
AM_CONDITIONAL(MAKE_MODULE_PYTHON, [test "$MAKE_MODULE_PYTHON" = "yes"])

##################################################
#	library options	    
##################################################

if test -z $NUM_THREADS ; then NUM_THREADS=4 ; fi
AH_TEMPLATE([NUM_THREADS],
            [Define to the default number of threads that the Vision 
             Workbench should spawn for multi-threaded operations.])
AC_ARG_WITH(num-threads,
  AC_HELP_STRING([--with-num-threads[=4]], [set the default number of processing threads for multi-threaded operations]),
  [NUM_THREADS=$withval])
AC_DEFINE_UNQUOTED([NUM_THREADS],[$NUM_THREADS])
AC_MSG_NOTICE([VW will use $NUM_THREADS THREADS by default])

##################################################
#	final processing	    
##################################################

VW_CPPFLAGS="$OTHER_CPPFLAGS"
VW_LDFLAGS="$OTHER_LDFLAGS"

if test "yes" = "$ENABLE_VERBOSE" ; then
  AC_MSG_NOTICE([using CPPFLAGS=$CPPFLAGS])
  AC_MSG_NOTICE([using VW_CPPFLAGS=$VW_CPPFLAGS])
  AC_MSG_NOTICE([using CFLAGS=$CFLAGS])
  AC_MSG_NOTICE([using CXXFLAGS=$CXXFLAGS])
  AC_MSG_NOTICE([using LDFLAGS=$LDFLAGS])
  AC_MSG_NOTICE([using VW_LDFLAGS=$VW_LDFLAGS])
  AC_MSG_NOTICE([using VW_NO_EXCEPTIONS=$no_exceptions])
fi

AC_SUBST(VW_CPPFLAGS)
AC_SUBST(VW_LDFLAGS)


##################################################
#	makefiles to process	    
##################################################

AC_CONFIG_FILES([                       \
  Makefile                              \
  src/Makefile                          \
  src/vw/Makefile                       \
  src/vw/tests/Makefile                 \
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
  src/vw/GPU/Makefile                   \
  src/vw/GPU/tests/Makefile             \
  src/vw/InterestPoint/Makefile         \
  src/vw/InterestPoint/tests/Makefile   \
  src/vw/Geometry/Makefile              \
  src/vw/Geometry/tests/Makefile        \
  Python/vw/Makefile                    \
  Python/Makefile                       \
])

AC_OUTPUT
