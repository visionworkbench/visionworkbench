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


AC_INIT([NASA Vision Workbench], [0.2], [vwb@nx.arc.nasa.gov], [VisionWorkbench])
AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE
AX_CONFIG_HEADER_PREFIX([src/vw/config.h], [VW_])


##################################################
#	configuration options
##################################################

# This should be first, so everything else can use it
AC_ARG_ENABLE(verbose,
  AC_HELP_STRING([--enable-verbose], [verbose configuration output @<:@no@:>@]),           
  [ ENABLE_VERBOSE=$enableval ],
  [ if test -z $"ENABLE_VERBOSE"; then ENABLE_VERBOSE=no; fi ]
)

if ! test -f config.options && test -f config/config.options.default; then
  if test "yes" = "$ENABLE_VERBOSE" ; then
    AC_MSG_NOTICE([creating default config.options file])
  fi
  cp config/config.options.default config.options
fi

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

CPPFLAGS="$CPPFLAGS -I\${top_srcdir}/src"

AC_ARG_ENABLE(debug,
  AC_HELP_STRING([--enable-debug], [generate debugging symbols @<:@no@:>@]),
  [ ENABLE_DEBUG=$enableval ], 
  [ if test -z "$ENABLE_DEBUG"; then ENABLE_DEBUG=no; fi ]
)

if test "yes" = "$ENABLE_DEBUG"; then
   AX_CFLAGS="$AX_CFLAGS -g -DDEBUG"
else
   AX_CFLAGS="$AX_CFLAGS -DNDEBUG"
fi

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
   AX_CFLAGS="$AX_CFLAGS -O4 -march=prescott -mtune=prescott -funroll-loops -fexpensive-optimizations -msse -msse2 -msse3 -mfpmath=sse -malign-double"
elif test "sse3" = "$ENABLE_OPTIMIZE"; then
   AX_CFLAGS="$AX_CFLAGS -O4 -funroll-loops -fexpensive-optimizations -msse -msse2 -msse3 -mfpmath=sse -malign-double"
else
   AX_CFLAGS="$AX_CFLAGS -O0"
fi

CFLAGS="$CFLAGS $AX_CFLAGS"
CXXFLAGS="$CXXFLAGS $AX_CFLAGS"

if test "yes" = "$ENABLE_VERBOSE" ; then
  AC_MSG_NOTICE([using CPPFLAGS=$CPPFLAGS])
  AC_MSG_NOTICE([using CFLAGS=$CFLAGS])
  AC_MSG_NOTICE([using CXXFLAGS=$CXXFLAGS])
fi

# We set up the include paths ourselves, so 
# we ask autoconf to get out of the way.
DEFAULT_INCLUDES=
AC_SUBST(DEFAULT_INCLUDES)

##################################################
#	installation options
##################################################

AC_ARG_ENABLE(custom-arch-name,
  AC_HELP_STRING([--enable-custom-arch-name=NAME], [custom architecture name, used for generating the default intallation prefix]),
  [ENABLE_CUSTOM_ARCH_NAME=$enableval]
)

if test ${prefix} = NONE; then
  if test ! -z ${ENABLE_CUSTOM_ARCH_NAME} ; then
    prefix=`pwd`/build/${ENABLE_CUSTOM_ARCH_NAME}
  else 
    prefix=`pwd`/build/${host}
  fi
fi

AC_MSG_NOTICE([using installation prefix ${prefix}])


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
  AC_HELP_STRING([--with-dist-remove-modules="RELPATH @<:@RELPATH...@:>@"], [remove certain modules from the distribution]),
  [DIST_REMOVE_MODULES=$withval]
)
AC_SUBST(DIST_REMOVE_MODULES)
AM_CONDITIONAL(DIST_REMOVE_MODULES, [test ! -z "$DIST_REMOVE_MODULES"])


##################################################
#	package checks	    
##################################################

AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AX_PROG_AR
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET

if test "yes" = "$ENABLE_VERBOSE"; then
  AC_MSG_NOTICE([using PKG_PATHS=$PKG_PATHS])
fi

AX_PKG_PTHREADS

AX_PKG_BOOST
AX_PKG_BOOST_LIB(PROGRAM_OPTIONS, [-lboost_program_options], [boost/program_options.hpp])
AX_PKG_BOOST_LIB(FILESYSTEM, [-lboost_filesystem], [boost/filesystem/path.hpp])
AX_PKG_BOOST_LIB(THREAD, [-lboost_thread], [boost/thread/thread.hpp])

AX_PKG_LAPACK
AM_CONDITIONAL(HAVE_PKG_LAPACK, [test "$HAVE_PKG_LAPACK" = "yes"])

AX_PKG(M, [], [-lm], [math.h])
AX_PKG(Z, [], [-lz], [zlib.h])
AX_PKG(PNG, [], [-lpng], [png.h])
AX_PKG(TIFF, [], [-ltiff], [tiff.h])
AX_PKG(JPEG, [], [-ljpeg], [jpeglib.h])

AX_PKG(X11, [], [-lXext -lX11 -lXmu -lXi -lXt], [])
AX_PKG(GL, [X11], [-lGL], [GL/gl.h])
AX_PKG(GLU, [GL], [-lGLU], [GL/glu.h])
AX_PKG(GLUT, [GL GLU], [-lglut], [GL/glut.h])
AX_PKG(FREEGLUT, [GL GLU], [-lglut], [GL/glut.h GL/freeglut_ext.h])
AX_PKG(GLEW, [X11 GL GLU M], [-lGLEW -lSM -lICE], [GL/glew.h])
AX_PKG(CG, [], [-lCg -lCgGL], [Cg/cg.h])

AX_PKG(OPENEXR, [Z PTHREADS], [-lImath -lHalf -lIex -lIlmImf], [ImfOutputFile.h], [OpenEXR])
AX_PKG(GDAL, [], [-lgdal], [gdal.h])
AX_PKG(PROJ4, [], [-lproj], [projects.h])

AX_PKG(VCL, [], [-lvcl], [vcl_algorithm.h], [vcl])
AX_PKG(VXL, [VCL], [-lvgl_algo -lvil_algo -lvnl_algo -lvul -lvgl -lvnl -lvil -lvcl -lnetlib], 
            [vil/vil_image_view.h], [core])

# These are here (instead of inside the PKG macro where they belong)
# for backwards compatability with older versions of automake.
AM_CONDITIONAL(HAVE_PKG_PNG, [test "$HAVE_PKG_PNG" = "yes"])
AM_CONDITIONAL(HAVE_PKG_TIFF, [test "$HAVE_PKG_TIFF" = "yes"])
AM_CONDITIONAL(HAVE_PKG_JPEG, [test "$HAVE_PKG_JPEG" = "yes"])
AM_CONDITIONAL(HAVE_PKG_OPENEXR, [test "$HAVE_PKG_OPENEXR" = "yes"])

##################################################
#	module definitions
##################################################

# These module source directories may or may not be provided.
AX_MODULE_SRC(STEREO, [src/vw/Stereo])

AX_MODULE(CORE, [src/vw/Core/libvwCore.la], yes, [M BOOST])
AX_MODULE(MATH, [src/vw/Math/libvwMath.la], yes, [CORE], [LAPACK])
AX_MODULE(IMAGE, [src/vw/Image/libvwImage.la], yes, [MATH], [VXL])
AX_MODULE(FILEIO, [src/vw/FileIO/libvwFileIO.la], yes, [IMAGE], [OPENEXR Z PNG TIFF JPEG])
AX_MODULE(VW, [src/vw/libvw.la], yes, [CORE MATH IMAGE], [FILEIO])
if test ${MAKE_MODULE_VW} != "yes"; then
  AC_MSG_ERROR([cannot build the Vision Workbench core!])
fi

AX_MODULE(CAMERA, [src/vw/Camera/libvwCamera.la], yes, [IMAGE MATH])
AX_MODULE(MOSAIC, [src/vw/Mosaic/libvwMosaic.la], yes, [VW PTHREADS BOOST_PROGRAM_OPTIONS BOOST_FILESYSTEM])
AX_MODULE(STEREO, [src/vw/Stereo/libvwStereo.la], $HAVE_PKG_STEREO_SRC, [STEREO_SRC VW PTHREADS BOOST_THREAD])
AX_MODULE(CARTOGRAPHY, [src/vw/Cartography/libvwCartography.la], yes, [VW GDAL PROJ4])

#AX_MODULE(INTERESTPOINT, InterestPoint, no, [VW])
#AX_MODULE(HDR, HDR, no, [VW])
#AX_MODULE(GPU, GPU, no, [GL GLU GLUT GLEW VW], [CG PTHREADS])
#AX_MODULE(DISPLAY, Display, no, [BOOST_THREAD GL FREEGLUT VW])
#AX_MODULE(TEXTURE, Texture, no, [VW])

# These are here (instead of inside the MODULE macro where they belong)
# for backwards compatability with older versions of automake.
AM_CONDITIONAL(MAKE_MODULE_CORE, [test "$MAKE_MODULE_CORE" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_MATH, [test "$MAKE_MODULE_MATH" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_IMAGE, [test "$MAKE_MODULE_IMAGE" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_FILEIO, [test "$MAKE_MODULE_FILEIO" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_CAMERA, [test "$MAKE_MODULE_CAMERA" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_VW, [test "$MAKE_MODULE_VW" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_MOSAIC, [test "$MAKE_MODULE_MOSAIC" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_STEREO, [test "$MAKE_MODULE_STEREO" = "yes"])
AM_CONDITIONAL(MAKE_MODULE_CARTOGRAPHY, [test "$MAKE_MODULE_CARTOGRAPHY" = "yes"])

#AM_CONDITIONAL(MAKE_MODULE_INTERESTPOINT, [test "$MAKE_MODULE_INTERESTPOINT" = "yes"])
#AM_CONDITIONAL(MAKE_MODULE_HDR, [test "$MAKE_MODULE_HDR" = "yes"])
#AM_CONDITIONAL(MAKE_MODULE_GPU, [test "$MAKE_MODULE_GPU" = "yes"])
#AM_CONDITIONAL(MAKE_MODULE_DISPLAY, [test "$MAKE_MODULE_DISPLAY" = "yes"])
#AM_CONDITIONAL(MAKE_MODULE_CAMERA, [test "$MAKE_MODULE_CAMERA" = "yes"])
#AM_CONDITIONAL(MAKE_MODULE_TEXTURE, [test "$MAKE_MODULE_TEXTURE" = "yes"])


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
])

if test $HAVE_PKG_STEREO_SRC = "yes"; then
  AC_CONFIG_FILES([                       \
    src/vw/Stereo/Makefile                \
    src/vw/Stereo/tests/Makefile          \
  ])
fi

AC_OUTPUT
