dnl __BEGIN_LICENSE__
dnl 
dnl Copyright (C) 2006 United States Government as represented by the
dnl Administrator of the National Aeronautics and Space Administration
dnl (NASA).  All Rights Reserved.
dnl 
dnl Copyright 2006 Carnegie Mellon University. All rights reserved.
dnl 
dnl This software is distributed under the NASA Open Source Agreement
dnl (NOSA), version 1.3.  The NOSA has been approved by the Open Source
dnl Initiative.  See the file COPYING at the top of the distribution
dnl directory tree for the complete NOSA document.
dnl 
dnl THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
dnl KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
dnl LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
dnl SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
dnl A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
dnl THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
dnl DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
dnl 
dnl __END_LICENSE__

dnl Usage: AX_CONFIG_HEADER_PREFIX(<filename>, <prefix>)
dnl Generates a configuration header file, adding a prefix to all symbols.
dnl This is a two-step process.  First we generate the usual header file 
dnl with a filename ending in ".pre".  Then we process that file, adding 
dnl the prefix to all symbolx, and copy it into the final location if it 
dnl has changed.
AC_DEFUN([AX_CONFIG_HEADER_PREFIX],
[
  AM_CONFIG_HEADER([$1.pre],
  [
    echo "/* $1.  Generated from $1.pre by config.status.  */" > "$1.new"
    echo "#ifndef __$2_CONFIG_H__" >> "$1.new"
    echo "#define __$2_CONFIG_H__" >> "$1.new"
    sed -e 's/#define /#define $2/' -e 's/#undef /#undef $2/' < "$1.pre" >> "$1.new"
    echo "#endif // __$2_CONFIG_H__" >> "$1.new"
    if test -f "$1" && diff "$1" "$1.new" > /dev/null ; then
      echo "config.status: $1 is unchanged"
    else
      echo "config.status: creating $1"
      cp "$1.new" "$1"
    fi
    rm -f "$1.new"
  ])
])

dnl Usage: AX_PROG_AR
dnl Provides a workaround for Mac OS X's unusual ar behavior, so that 
dnl it's possible to build static convenience libraries using libtool 
dnl on that platform.  Basically, if we're on a Mac we introduce a 
dnl wrapper shell script for ar that detects when we're creating an 
dnl empty library and creates it by hand.  In all other cases it just 
dnl relays the arguments to the user's AR.
AC_DEFUN([AX_PROG_AR],
[
  AC_REQUIRE([AM_AUX_DIR_EXPAND])
  AC_MSG_CHECKING([whether to use ar.sh wrapper script])
  if test $host_vendor = "apple"; then
    AC_MSG_RESULT([yes])
    if test -z "$AR" ; then
      ax_ar_prog=ar;
    else
      ax_ar_prog=$AR;
    fi
    AR=$am_aux_dir/ar.sh
    cat > $AR <<_AX_EOF
#!/bin/sh
if test -z "[\${3}]" ; then
  echo '!<arch>' > [\${2}]
else
  $ax_ar_prog [\${*}]
fi
_AX_EOF
    chmod +x $am_aux_dir/ar.sh
  else
    AC_MSG_RESULT([no])
  fi
])

dnl Usage: AX_FIND_FILES(<filenames>, <search paths>)
dnl Looks to see if all the given filenames (relative paths) are accessible 
dnl from one of the given base paths.  Returns the path or the empty string 
dnl in ${ax_find_files_path}.
AC_DEFUN([AX_FIND_FILES],
[
  ax_find_files_path="";
  for path in $2; do
    ax_find_files_passed=yes
    for filename in $1; do
      pathname="$path/$filename"
      if test "$ENABLE_VERBOSE" = "yes"; then
        AC_MSG_CHECKING([for ${pathname}])
      fi
      ax_find_files_paths=`ls $pathname 2>/dev/null`
      if test ! -z "$ax_find_files_paths" ; then
        if test "$ENABLE_VERBOSE" = "yes"; then
          AC_MSG_RESULT([found])
        fi
      else
        if test "$ENABLE_VERBOSE" = "yes"; then
          AC_MSG_RESULT([not found])
        fi
        ax_find_files_passed=no
        break
      fi
    done
    if test "$ax_find_files_passed" = "yes"; then
      ax_find_files_path="$path"
      break
    fi
  done
])


dnl Usage: AX_PKG(<name>, <dependencies>, <libraries>, <headers>[, <relative include path>])
AC_DEFUN([AX_PKG],
[
  AC_ARG_WITH(translit($1,`A-Z',`a-z'),
    AC_HELP_STRING([--with-]translit($1,`A-Z',`a-z'), [enable searching for the $1 package @<:@auto@:>@]),
    [ HAVE_PKG_$1=$withval ]
  )

  if test x"$ENABLE_VERBOSE" = "xyes"; then
    AC_MSG_CHECKING([for package $1 in current paths])
  else
    AC_MSG_CHECKING([for package $1])
  fi

  AC_LANG_ASSERT(C++)

  # We can skip searching if we're already at "no"
  if test "no" = "$HAVE_PKG_$1"; then
    AC_MSG_RESULT([no (disabled by user)])

  else
    # Test for and inherit libraries from dependencies
    PKG_$1_LIBS="$3"
    for x in $2; do
      ax_pkg_have_dep=HAVE_PKG_${x}
      if test "${!ax_pkg_have_dep}" = "yes"; then
        ax_pkg_dep_libs="PKG_${x}_LIBS"
        PKG_$1_LIBS="$PKG_$1_LIBS ${!ax_pkg_dep_libs}"
        unset ax_pkg_dep_libs
      else
        unset PKG_$1_LIBS
        HAVE_PKG_$1="no"
        break
      fi
    done

    if test "x$HAVE_PKG_$1" = "xno" ; then
      AC_MSG_RESULT([no (needs $x)])

    # We skip the search if the user has been explicit about "yes"
    elif test "x$HAVE_PKG_$1" = "xyes" ; then
      AC_MSG_RESULT([yes (using user-supplied flags)])

      # Add package-specific cppflags/ldflags to what we have so far.
      # This is necessary for systems that have one version of a library 
      # installed in a default location, and another that they want to use 
      # instead in a separate location.
      if test -n ${PKG_$1_CPPFLAGS} ; then
        VW_CPPFLAGS="${PKG_$1_CPPFLAGS} ${VW_CPPFLAGS}"
      fi
      if test -n ${PKG_$1_LDFLAGS} ; then
        VW_LDFLAGS="${PKG_$1_LDFLAGS} ${VW_LDFLAGS}"
      fi

    # Otherwise we look for a path that contains the needed headers and libraries
    else

      if test "x$ENABLE_VERBOSE" = "yes"; then
	    AC_MSG_RESULT([searching...])
      fi

      if test -n "${HAVE_PKG_$1}" && test "${HAVE_PKG_$1}" != "yes" && test "${HAVE_PKG_$1}" != "no"; then
        PKG_PATHS_$1=${HAVE_PKG_$1}
      else
        PKG_PATHS_$1=${PKG_PATHS}
      fi

      HAVE_PKG_$1=no

      ax_pkg_old_libs="$LIBS"
      LIBS="$PKG_$1_LIBS $LIBS"
      for path in none $PKG_PATHS_$1; do
	    ax_pkg_old_cppflags=$CPPFLAGS
	    ax_pkg_old_ldflags=$LDFLAGS
	    ax_pkg_old_vw_cppflags=$VW_CPPFLAGS
	    ax_pkg_old_vw_ldflags=$VW_LDFLAGS
	    echo > conftest.h
	    for header in $4 ; do
	      echo "#include <$header>" >> conftest.h
	    done
	    CPPFLAGS="$ax_pkg_old_cppflags $VW_CPPFLAGS"
	    LDFLAGS="$ax_pkg_old_ldflags $VW_LDFLAGS"
	    if test "$path" != "none"; then
	      if test x"$ENABLE_VERBOSE" = "xyes"; then
	        AC_MSG_CHECKING([for package $1 in $path])
	      fi
          if test -z "$5"; then
            VW_CPPFLAGS="-I$path/include $VW_CPPFLAGS"
	      else
	        VW_CPPFLAGS="-I$path/include/$5 $VW_CPPFLAGS"
          fi
	      CPPFLAGS="$ax_pkg_old_cppflags $VW_CPPFLAGS"
	      AC_LINK_IFELSE(
	        AC_LANG_PROGRAM([#include "conftest.h"],[]),
	        [ HAVE_PKG_$1=yes ; AC_MSG_RESULT([yes]) ; break ] )
          # Sometimes we'll have /foo/lib64 and /foo/lib confusion on
          # 64-bit machines, so accept both if one doesn't appear.
          if test -d $path/${AX_LIBDIR}; then
	        VW_LDFLAGS="-L$path/${AX_LIBDIR} $VW_LDFLAGS"
          elif test x"${AX_LIBDIR}" = "xlib64"; then
            VW_LDFLAGS="-L$path/${AX_OTHER_LIBDIR} $VW_LDFLAGS"
          fi
          LDFLAGS="$ax_pkg_old_ldflags $VW_LDFLAGS"
        fi
        AC_LINK_IFELSE(
          AC_LANG_PROGRAM([#include "conftest.h"],[]),
          [ HAVE_PKG_$1=yes ; AC_MSG_RESULT([yes]) ; break ] )
        if test x"$ENABLE_VERBOSE" = "xyes"; then
          AC_MSG_RESULT([no])
        fi
        CPPFLAGS=$ax_pkg_old_cppflags
        LDFLAGS=$ax_pkg_old_ldflags
        VW_CPPFLAGS=$ax_pkg_old_vw_cppflags
        VW_LDFLAGS=$ax_pkg_old_vw_ldflags
      done
      CPPFLAGS=$ax_pkg_old_cppflags
      LDFLAGS=$ax_pkg_old_ldflags
      LIBS=$ax_pkg_old_libs

      if test "x$HAVE_PKG_$1" = "xno" -a "x$ENABLE_VERBOSE" != "xyes"; then
	    AC_MSG_RESULT([no (not found)])
      fi

    fi

  fi
  if test "${HAVE_PKG_$1}" = "yes" ; then
    ax_have_pkg_bool=1
  else
    ax_have_pkg_bool=0
    PKG_$1_LIBS=
  fi
  AC_DEFINE_UNQUOTED([HAVE_PKG_$1],
                     [$ax_have_pkg_bool],
                     [Define to 1 if the $1 package is available.])

  AC_SUBST(PKG_$1_LIBS)
  AC_SUBST(HAVE_PKG_$1)

  if test x"$ENABLE_VERBOSE" == "xyes"; then
    AC_MSG_NOTICE([HAVE_PKG_$1 = ${HAVE_PKG_$1}])
    AC_MSG_NOTICE([PKG_$1_LIBS= $PKG_$1_LIBS])
    AC_MSG_NOTICE([VW_CPPFLAGS= $VW_CPPFLAGS])
    AC_MSG_NOTICE([VW_LDFLAGS= $VW_LDFLAGS])
  fi
])


# Usage: AX_PKG_BOOST
AC_DEFUN([AX_PKG_BOOST],
[
  AC_MSG_CHECKING(for package BOOST)
  if test "$ENABLE_VERBOSE" = "yes"; then
    AC_MSG_RESULT([])
  fi

  AC_LANG_ASSERT(C++)

  if test -n "${HAVE_PKG_BOOST}" && test "${HAVE_PKG_BOOST}" != "yes" && test "${HAVE_PKG_BOOST}" != "no"; then
    PKG_PATHS_BOOST=${HAVE_PKG_BOOST}
    unset HAVE_PKG_BOOST
  else
    PKG_PATHS_BOOST=${PKG_PATHS}
  fi

  # Skip testing if the user has overridden
  if test -z ${HAVE_PKG_BOOST}; then

    PKG_BOOST_LIBS=
    HAVE_PKG_BOOST=no

    for ax_boost_base_path in $PKG_PATHS_BOOST; do
      # First look for a system-style installation
      if test "$ENABLE_VERBOSE" = "yes"; then
        AC_MSG_CHECKING([for system-style boost in ${ax_boost_base_path}])
      fi
      if test -d "${ax_boost_base_path}/include/boost" ; then
        PKG_BOOST_INCDIR="${ax_boost_base_path}/include"
        PKG_BOOST_LIBDIR="${ax_boost_base_path}/${AX_LIBDIR}"
        # In case it's not in lib64 despite specifying lib64...
        if test ! -d $PKG_BOOST_LIBDIR -a x"${AX_LIBDIR}" = "xlib64"; then
          PKG_BOOST_LIBDIR="${ax_boost_base_path}/${AX_OTHER_LIBDIR}"
        fi
        HAVE_PKG_BOOST="yes"
        if test "$ENABLE_VERBOSE" = "yes"; then
          AC_MSG_RESULT([found])
        fi
        break
      else
        if test "$ENABLE_VERBOSE" = "yes"; then
          AC_MSG_RESULT([not found])
        fi
      fi
      # Next look for a default-style installation
      if test "$ENABLE_VERBOSE" = "yes"; then
        AC_MSG_CHECKING([for default-style boost in ${ax_boost_base_path}])
      fi
      for ax_boost_inc_path in `ls -d ${ax_boost_base_path}/include/boost-* 2> /dev/null` ; do
        # At the moment we greedily accept the first one we find, regardless of version
        PKG_BOOST_INCDIR="${ax_boost_inc_path}"
        PKG_BOOST_LIBDIR="${ax_boost_base_path}/lib"
        HAVE_PKG_BOOST="yes"
        if test "$ENABLE_VERBOSE" = "yes"; then
          AC_MSG_RESULT([found])
        fi
        break 2
      done
      if test "$ENABLE_VERBOSE" = "yes"; then
        AC_MSG_RESULT([not found])
      fi
    done
  fi

  if test "${HAVE_PKG_BOOST}" = "yes" ; then
    ax_pkg_old_vw_cppflags=$VW_CPPFLAGS
    ax_pkg_old_vw_ldflags=$VW_LDFLAGS
    ax_pkg_old_cppflags=$CPPFLAGS
    ax_pkg_old_ldflags=$LDFLAGS
    ax_pkg_old_libs=$LIBS
    while true ; do
      # First see if the current paths are sufficient
      if test "x${ENABLE_VERBOSE}" = "xyes" ; then
	AC_MSG_CHECKING([whether current paths are sufficient...])
      fi
      AC_LINK_IFELSE( AC_LANG_PROGRAM([#include <boost/version.hpp>],[]), [ax_result=yes], [ax_result=no] )
      if test "x${ENABLE_VERBOSE}" = "xyes" ; then
	AC_MSG_RESULT([$ax_result])
      fi
      if test "$ax_result" = "yes" ; then break ; fi
      # Try it with just the include path
      VW_CPPFLAGS="-I${PKG_BOOST_INCDIR} $VW_CPPFLAGS"
      CPPFLAGS="$ax_pkg_old_cppflags $VW_CPPFLAGS"
      if test "x${ENABLE_VERBOSE}" = "xyes" ; then
	AC_MSG_CHECKING([whether adding the include path is sufficient...])
      fi
      AC_LINK_IFELSE( AC_LANG_PROGRAM([#include <boost/version.hpp>],[]), [ax_result=yes], [ax_result=no] )
      if test "x${ENABLE_VERBOSE}" = "xyes" ; then
	AC_MSG_RESULT([$ax_result])
      fi
      if test "$ax_result" = "yes" ; then break ; fi
      # Finally, try it with the linker path
      VW_LDFLAGS="-L${PKG_BOOST_LIBDIR} $VW_LDFLAGS"
      LDFLAGS="$ax_pkg_old_ldflags $VW_LDFLAGS"
      if test "x${ENABLE_VERBOSE}" = "xyes" ; then
	AC_MSG_CHECKING([whether adding the include and linker paths works...])
      fi
      AC_LINK_IFELSE( AC_LANG_PROGRAM([#include <boost/version.hpp>],[]), [ax_result=yes], [ax_result=no] )
      if test "x${ENABLE_VERBOSE}" = "xyes" ; then
	AC_MSG_RESULT([$ax_result])
      fi
      if test "$ax_result" = "yes" ; then break ; fi
      # The detected version of boost seems to be invalid!
      HAVE_PKG_BOOST="no"
      VW_CPPFLAGS="$ax_pkg_old_vw_cppflags"
      VW_LDFLAGS="$ax_pkg_old_vw_ldflags"
      unset PKG_BOOST_INCDIR
      unset PKG_BOOST_LIBDIR
      break
    done
  fi
  CPPFLAGS="$ax_pkg_old_cppflags"
  LDFLAGS="$ax_pkg_old_ldflags"

  if test "${HAVE_PKG_BOOST}" = "yes" ; then
    ax_have_pkg_bool=1
  else
    ax_have_pkg_bool=0
  fi
  AC_DEFINE_UNQUOTED([HAVE_PKG_BOOST],
                     [$ax_have_pkg_bool],
                     [Define to 1 if the BOOST package is available.])

  AC_SUBST(HAVE_PKG_BOOST)

  AC_LANG_CONFTEST(
  [AC_LANG_PROGRAM([[
#include <iostream>
#include <boost/version.hpp>
#define STR2(s) #s
#define STR(s) STR2(s)
]],[[
std::cout << STR(BOOST_VERSION);
]])])
  $CXX $VW_CPPFLAGS -I${PKG_BOOST_INCDIR} -o conftest conftest.$ac_ext
  BOOST_VERSION=`./conftest`
  AC_DEFINE_UNQUOTED([BOOST_VERSION],
                     [$BOOST_VERSION],
                     [The version of Boost with which the Vision Workbench was built.])

  AH_VERBATIM([_VW_CHECK_BOOST_VERSION],
[// Check to make sure the user is using the same version of Boost
// headers that the Vision Workbench was built with.
#include <boost/version.hpp>
#if BOOST_VERSION != VW_BOOST_VERSION
#error You are using a different version of Boost than you used to build the Vision Workbench!
#endif
])

  if test "$ENABLE_VERBOSE" = "yes"; then
    AC_MSG_NOTICE([HAVE_PKG_BOOST= $HAVE_PKG_BOOST])
    AC_MSG_NOTICE([VW_CPPFLAGS= $VW_CPPFLAGS])
    AC_MSG_NOTICE([VW_LDFLAGS= $VW_LDFLAGS])
  else
    AC_MSG_RESULT([$HAVE_PKG_BOOST])
  fi

])


# Usage: AX_PKG_LAPACK
#
# TODO: Add support for other sources of LAPACK and BLAS, such as
# ATLAS and the Intel math libraries.  For people who don't have any
# of these third party libraries installed, we need to fall back to
# compiling BLAS and LAPACK ourselves.
AC_DEFUN([AX_PKG_LAPACK],
[

  # If we are running MacOS X, we can use Apple's vecLib framework to
  # provide us with LAPACK and BLAS routines.
  if test $host_vendor = apple; then
    AC_MSG_CHECKING(for package LAPACK)
    if test "$ENABLE_VERBOSE" = "yes"; then
      AC_MSG_RESULT([])
    fi

    HAVE_PKG_LAPACK="yes"
    PKG_LAPACK_LIBS="$VW_LDFLAGS -framework vecLib"

    if test "$ENABLE_VERBOSE" = "yes"; then
      AC_MSG_RESULT([found])
    fi
    if test "${HAVE_PKG_LAPACK}" = "yes" ; then
      ax_have_pkg_bool=1
    else
      ax_have_pkg_bool=0
    fi
    AC_DEFINE_UNQUOTED([HAVE_PKG_LAPACK],
                       [$ax_have_pkg_bool],
                       [Define to 1 if the LAPACK package is available.])

    AC_SUBST(HAVE_PKG_LAPACK)

    if test "$ENABLE_VERBOSE" = "yes"; then
      AC_MSG_NOTICE([HAVE_PKG_LAPACK = ${HAVE_PKG_LAPACK}])
      AC_MSG_NOTICE([PKG_LAPACK_LIBS = ${PKG_LAPACK_LIBS}])
      AC_MSG_NOTICE([VW_CPPFLAGS = ${VW_CPPFLAGS}])
      AC_MSG_NOTICE([VW_LDFLAGS = ${VW_LDFLAGS}])
    else
      AC_MSG_RESULT([${HAVE_PKG_LAPACK}])
    fi  

  # For all other platforms, we search for static LAPACK libraries
  # in the conventional manner.
  else
    # First check for CLAPACK
    AX_PKG(CLAPACK, [], [-lclapack -lblas -lf2c], [])
    if test "$HAVE_PKG_CLAPACK" = "no"; then
      # Otherwise check for standard LAPACK
      AC_MSG_NOTICE(["CLAPACK not found, trying standard LAPACK."])
      AX_PKG(SLAPACK, [], [-llapack -lblas], [])

      if test "$HAVE_PKG_SLAPACK" = "no"; then
        # Some newer boxes require -lgfortran, so try that too
        AC_MSG_NOTICE(["trying standard LAPACK with -lgfortran."])
        AX_PKG(FLAPACK, [], [-llapack -lblas -lgfortran], [])

        if test "$HAVE_PKG_FLAPACK" = "no"; then
          # On some systems, BLAS and LAPACK are installed in different places
          AC_MSG_NOTICE(["trying to find BLAS and LAPACK seperately."])
          AX_PKG(STANDALONE_BLAS, [], [-lblas], [])
          AX_PKG(STANDALONE_LAPACK, [], [-llapack], [])
          AX_PKG(STANDALONE_LAPACK_AND_BLAS, [STANDALONE_LAPACK STANDALONE_BLAS], [], [])

	  if test "$HAVE_PKG_STANDALONE_LAPACK_AND_BLAS" = "no"; then
            # On some systems, F2C, FBLAS and FLAPACK are installed in different places
            AC_MSG_NOTICE(["trying to find F2C, FBLAS, and FLAPACK seperately."])
            AX_PKG(STANDALONE_F2C, [], [-lf2c], [])
            AX_PKG(STANDALONE_FBLAS, [STANDALONE_F2C], [-lblas], [])
            AX_PKG(STANDALONE_FLAPACK, [STANDALONE_F2C], [-llapack], [])
            AX_PKG(LAPACK, [STANDALONE_FLAPACK STANDALONE_FBLAS STANDALONE_F2C], [], [])
	  else 
            AX_PKG(LAPACK, [STANDALONE_LAPACK_AND_BLAS], [], [])
	  fi # FBLAS and FLAPACK
        else
          AX_PKG(LAPACK, [FLAPACK], [], [])
        fi # BLAS and LAPACK
      else
        AX_PKG(LAPACK, [SLAPACK], [], [])
      fi # fortan LAPACK
    else
      AX_PKG(LAPACK, [CLAPACK], [], [])
    fi
  fi
])

# Usage: AX_PKG_GL
#
# TODO: Add support for other sources of GL and BLAS, such as
# ATLAS and the Intel math libraries.  For people who don't have any
# of these third party libraries installed, we need to fall back to
# compiling BLAS and GL ourselves.
AC_DEFUN([AX_PKG_GL],
[

  # If we are running MacOS X, we can use Apple's vecLib framework to
  # provide us with GL and BLAS routines.
  if test $host_vendor = apple; then
    AC_MSG_CHECKING(for package GL)
    if test "$ENABLE_VERBOSE" = "yes"; then
      AC_MSG_RESULT([])
    fi

    HAVE_PKG_GL="yes"
    PKG_GL_LIBS="$VW_LDFLAGS -framework OpenGL -framework GLUT"

    if test "$ENABLE_VERBOSE" = "yes"; then
      AC_MSG_RESULT([found])
    fi
    if test "${HAVE_PKG_GL}" = "yes" ; then
      ax_have_pkg_bool=1
    else
      ax_have_pkg_bool=0
    fi
    AC_DEFINE_UNQUOTED([HAVE_PKG_GL],
                       [$ax_have_pkg_bool],
                       [Define to 1 if the GL package is available.])

    AC_SUBST(HAVE_PKG_GL)

    if test "$ENABLE_VERBOSE" = "yes"; then
      AC_MSG_NOTICE([HAVE_PKG_GL = ${HAVE_PKG_GL}])
      AC_MSG_NOTICE([PKG_GL_LIBS = ${PKG_GL_LIBS}])
      AC_MSG_NOTICE([VW_CPPFLAGS = ${VW_CPPFLAGS}])
      AC_MSG_NOTICE([VW_LDFLAGS = ${VW_LDFLAGS}])
    else
      AC_MSG_RESULT([${HAVE_PKG_GL}])
    fi  

  # For all other platforms, we search for static GL libraries
  # in the conventional manner.
  else
    AX_PKG(GL, [X11], [-lGL -lGLU -lglut], [GL/gl.h GL/glu.h GL/glut.h])
  fi

])

dnl Here's a new version of AX_PKG_BOOST_LIB designed to find the 
dnl multithreaded boost libraries and boost libraries that are just weirdly
dnl named in general. Boost libraries follow a weird naming convention 
dnl that makes our old logic not work. You can't just add -mt to the old 
dnl library you're looking for, because the -compiler part comes first. 
dnl IE, the non-multithreaded library would be named libboost_X-gcc41.so, 
dnl and the multithreaded library would be named libboost_X-gcc41-mt.so.
dnl
dnl For that reason, we've added an environment variable:
dnl BOOST_LIBRARIES_SUFFIX. The function here tries to find a version of 
dnl Boost with the string in said variable somewhere inside the Boost 
dnl library names, but after the initial name of the library (specified 
dnl as the second parameter to this function). A blank value will give 
dnl normal behavior.
# Usage: AX_PKG_BOOST_LIB(<name>, <libraries>, <header>)
AC_DEFUN([AX_PKG_BOOST_LIB],
[
  AC_MSG_CHECKING(for package BOOST_$1)
  if test "$ENABLE_VERBOSE" = "yes"; then
    AC_MSG_RESULT([])
  fi

  AC_LANG_ASSERT(C++)

  # Skip testing if the user has overridden
  if test -z ${HAVE_PKG_BOOST_$1}; then

    HAVE_PKG_BOOST_$1=no

    # Check for general Boost presence
    if test "x${HAVE_PKG_BOOST}" = "xyes" ; then
      # Check for required headers
      AX_FIND_FILES([$3],[${PKG_BOOST_INCDIR}])
      if test ! -z "$ax_find_files_path" ; then
        # Check for required libraries with no suffix, aside from the one 
        # given by environment variable.
        AX_FIND_FILES([`echo $2 | sed "s/-l\([[^[:space:]]]*\)/lib\1${BOOST_LIBRARIES_SUFFIX}.*/g"`],[$PKG_BOOST_LIBDIR])
        if test ! -z "$ax_find_files_path" ; then
          HAVE_PKG_BOOST_$1="yes"
          PKG_BOOST_$1_LIBS="$2${BOOST_LIBRARIES_SUFFIX}"
        else
          # Check for required libraries with some suffix. We have to check
          # for both a suffix before ${BOOST_LIBRARIES_SUFFIX} (pre-suffix) 
          # and a suffix after (for example) the -mt (post-suffix), because 
          # boost likes to stick the name of the compiler before the -mt.
          # Extremely annoying.

          ax_pkg_boost_lib=`echo $2 | awk '{print [$]1}' | sed 's/-l\([[^[:space:]-]]*\).*/lib\1/g'`
          ax_pkg_boost_file=`ls ${PKG_BOOST_LIBDIR}/${ax_pkg_boost_lib}-*${BOOST_LIBRARIES_SUFFIX}* | head -n 1 | sed "s,^${PKG_BOOST_LIBDIR}/\(.*\),\1,"`

          # The pre-suffix.
          ax_pkg_boost_presuffix=`echo ${ax_pkg_boost_file} | sed "s/${ax_pkg_boost_lib}\([[^.]]*\)${BOOST_LIBRARIES_SUFFIX}.*/\1/"`

          # The post-suffix.
          ax_pkg_boost_postsuffix=`echo ${ax_pkg_boost_file} | sed "s/${ax_pkg_boost_lib}${ax_pkg_boost_presuffix}${BOOST_LIBRARIES_SUFFIX}\([[^.]]*\).*/\1/"`

          AX_FIND_FILES([`echo $2 | sed "s/-l\([[^[:space:]]]*\)/lib\1${ax_pkg_boost_presuffix}${BOOST_LIBRARIES_SUFFIX}${ax_pkg_boost_postsuffix}.*/g"`],[$PKG_BOOST_LIBDIR])
          if test ! -z $ax_find_files_path ; then
            HAVE_PKG_BOOST_$1="yes"
            PKG_BOOST_$1_LIBS=`echo $2 | sed "s/[[^ ]]*/&${ax_pkg_boost_presuffix}${BOOST_LIBRARIES_SUFFIX}${ax_pkg_boost_postsuffx}/g"`
          fi
        fi
      fi
    fi
  fi

  ax_pkg_old_vw_cppflags=$VW_CPPFLAGS
  ax_pkg_old_vw_ldflags=$VW_LDFLAGS
  ax_pkg_old_cppflags=$CPPFLAGS
  ax_pkg_old_ldflags=$LDFLAGS
  ax_pkg_old_libs=$LIBS
  while true ; do
    echo > conftest.h
    for header in $3 ; do
      echo "#include <$header>" >> conftest.h
    done
    # First see if the current paths are sufficient
    if test "x${ENABLE_VERBOSE}" = "xyes" ; then
      AC_MSG_CHECKING([whether current paths are sufficient...])
    fi
    CPPFLAGS="$ax_pkg_old_cppflags $VW_CPPFLAGS"
    LDFLAGS="$ax_pkg_old_ldflags $VW_LDFLAGS"
    LIBS="$PKG_BOOST_$1_LIBS $ax_pkg_old_libs"
    AC_LINK_IFELSE( AC_LANG_PROGRAM([#include "conftest.h"],[]), [ax_result=yes], [ax_result=no] )
    if test "x${ENABLE_VERBOSE}" = "xyes" ; then
      AC_MSG_RESULT([$ax_result])
    fi
    if test "$ax_result" = "yes" ; then break ; fi
    # Try it with just the include path
    if test "x${ENABLE_VERBOSE}" = "xyes" ; then
      AC_MSG_CHECKING([whether adding the include path is sufficient...])
    fi
    VW_CPPFLAGS="-I${PKG_BOOST_INCDIR} $VW_CPPFLAGS"
    CPPFLAGS="$ax_pkg_old_cppflags $VW_CPPFLAGS"
    AC_LINK_IFELSE( AC_LANG_PROGRAM([#include "conftest.h"],[]), [ax_result=yes], [ax_result=no] )
    if test "x${ENABLE_VERBOSE}" = "xyes" ; then
      AC_MSG_RESULT([$ax_result])
    fi
    if test "$ax_result" = "yes" ; then break ; fi
    # Finally, try it with the linker path
    if test "x${ENABLE_VERBOSE}" = "xyes" ; then
      AC_MSG_CHECKING([whether adding the include and linker paths works...])
    fi
    VW_LDFLAGS="-L${PKG_BOOST_LIBDIR} $VW_LDFLAGS"
    LDFLAGS="$ax_pkg_old_ldflags $VW_LDFLAGS"
    AC_LINK_IFELSE( AC_LANG_PROGRAM([#include "conftest.h"],[]), [ax_result=yes], [ax_result=no] )
    if test "x${ENABLE_VERBOSE}" = "xyes" ; then
      AC_MSG_RESULT([$ax_result])
    fi
    if test "$ax_result" = "yes" ; then break ; fi
    # The detected version of boost seems to be invalid!
    HAVE_PKG_BOOST_$1="no"
    VW_CPPFLAGS="$ax_pkg_old_vw_cppflags"
    VW_LDFLAGS="$ax_pkg_old_vw_ldflags"
    break
  done

  CPPFLAGS="$ax_pkg_old_cppflags"
  LDFLAGS="$ax_pkg_old_ldflags"
  LIBS="$ax_pkg_old_libs"

  if test "${HAVE_PKG_BOOST_$1}" = "yes" ; then
    ax_have_pkg_bool=1
  else
    ax_have_pkg_bool=0
  fi
  AC_DEFINE_UNQUOTED([HAVE_PKG_BOOST_$1],
                     [$ax_have_pkg_bool],
                     [Define to 1 if the BOOST_$1 package is available.])

  AC_SUBST(HAVE_PKG_BOOST_$1)
  AC_SUBST(PKG_BOOST_$1_LIBS)

  if test "$ENABLE_VERBOSE" = "yes"; then
    AC_MSG_NOTICE([HAVE_PKG_BOOST_$1 = ${HAVE_PKG_BOOST_$1}])
    AC_MSG_NOTICE([PKG_BOOST_$1_LIBS= $PKG_BOOST_$1_LIBS])
  else
    AC_MSG_RESULT([${HAVE_PKG_BOOST_$1}])
  fi
])

dnl Usage: AX_PKG_PTHREADS
AC_DEFUN([AX_PKG_PTHREADS],
[
  AC_MSG_CHECKING(for package PTHREADS)
  if test "$ENABLE_VERBOSE" = "yes"; then
    AC_MSG_RESULT([])
  fi

  AC_LANG_PUSH(C)
  HAVE_PKG_PTHREADS=no

  ax_pkg_pthreads_cppflags_options="none -pthread"
  ax_pkg_pthreads_ldflags_options="-lpthread none"

  for ax_pkg_pthreads_ldflags in $ax_pkg_pthreads_ldflags_options; do
    if test "$ax_pkg_pthreads_ldflags" = "none" ; then
      PKG_PTHREADS_LDFLAGS=""
    else
      PKG_PTHREADS_LDFLAGS=$ax_pkg_pthreads_ldflags
    fi
    for ax_pkg_pthreads_cppflags in $ax_pkg_pthreads_cppflags_options; do
      if test "$ax_pkg_pthreads_cppflags" = "none" ; then
        PKG_PTHREADS_CPPFLAGS=""
      else
        PKG_PTHREADS_CPPFLAGS=$ax_pkg_pthreads_cppflags
      fi

      ax_pkg_pthreads_save_CFLAGS="$CFLAGS"
      ax_pkg_pthreads_save_LDFLAGS="$LDFLAGS"
      CFLAGS="$CFLAGS $PKG_PTHREADS_CPPFLAGS"
      LDFLAGS="$PKG_PTHREADS_LDFLAGS $LDFLAGS"

      if test "$ENABLE_VERBOSE" = "yes" ; then
        AC_MSG_CHECKING([whether pthreads work with flags: \"$CFLAGS\" : \"$LDFLAGS\"])
      fi

      AC_TRY_LINK([#include <pthread.h>],
                  [pthread_t th; pthread_create(0,0,0,0);],
                  [HAVE_PKG_PTHREADS=yes])

      CFLAGS="$ax_pkg_pthreads_save_CFLAGS"
      LDFLAGS="$ax_pkg_pthreads_save_LDFLAGS"

      if test "$ENABLE_VERBOSE" = "yes" ; then
        AC_MSG_RESULT($HAVE_PKG_PTHREADS)
      fi

      if test "$HAVE_PKG_PTHREADS" = "yes"; then
        break 2;
      fi
    done
  done

  AC_LANG_POP(C)

  if test "$HAVE_PKG_PTHREADS" = "yes" ; then
    CFLAGS="$CFLAGS $PKG_PTHREADS_CPPFLAGS"
    CXXFLAGS="$CXXFLAGS $PKG_PTHREADS_CPPFLAGS"
    PKG_PTHREADS_LIBS="$PKG_PTHREADS_LDFLAGS"
  fi

  if test "${HAVE_PKG_PTHREADS}" = "yes" ; then
    ax_have_pkg_bool=1
  else
    ax_have_pkg_bool=0
  fi
  AC_DEFINE_UNQUOTED([HAVE_PKG_PTHREADS],
                     [$ax_have_pkg_bool],
                     [Define to 1 if the PTHREADS package is available.])

  AC_SUBST(HAVE_PKG_PTHREADS)
  AC_SUBST(PKG_PTHREADS_LIBS)

  if test "$ENABLE_VERBOSE" = "yes"; then
    AC_MSG_NOTICE([HAVE_PKG_PTHREADS = ${HAVE_PKG_PTHREADS}])
    AC_MSG_NOTICE([PKG_PTHREADS_LIBS = ${PKG_PTHREADS_LIBS}])
    AC_MSG_NOTICE([CFLAGS= $CFLAGS])
    AC_MSG_NOTICE([CXXFLAGS= $CXXFLAGS])
  else
    AC_MSG_RESULT([${HAVE_PKG_PTHREADS}])
  fi
])


# Usage: AX_MODULE(<name>, <directory>, <library>, <default>, <prerequisites>, <required dependencies>[, <optional dependencies>])
AC_DEFUN([AX_MODULE],
[
  # Silently ignore modules that don't exist in this distribution
  if test -d $2 ; then

    HAVE_PKG_$1_SRC=yes

    AC_ARG_ENABLE([module-]translit($1,`A-Z',`a-z'),
      AC_HELP_STRING([--enable-module-]translit($1,`A-Z',`a-z'), [enable the $1 module @<:@$4@:>@]), 
      [ ENABLE_MODULE_$1=$enableval ],
      [ if test x"$ENABLE_MODULE_$1" = x; then ENABLE_MODULE_$1=`/bin/echo -n $4 | tr [A-Z] [a-z]` ; fi ]
    )

    AC_MSG_CHECKING([whether to build module $1])
    ax_module_enable=$ENABLE_MODULE_$1

    if test "$ax_module_enable" != "yes" ; then
      AC_MSG_RESULT([no (disabled)])
    fi

    ax_libs=""

    # Check for prerequisites
    if test "$ax_module_enable" = "yes" ; then
      for ax_dependency in $5 ; do
        ax_dependency_have="HAVE_PKG_${ax_dependency}"
        if test x"${!ax_dependency_have}" = "xyes"; then
          ax_dep_libs="PKG_${ax_dependency}_LIBS"
          ax_libs="${ax_libs} ${!ax_dep_libs}"
	else
          AC_MSG_RESULT([no])
          AC_MSG_NOTICE([warning: unable to build requested module $1 (needs ${ax_dependency})!])
          ax_module_enable=no;
          break;
        fi
      done
    fi

    # Check for required dependencies
    if test "$ax_module_enable" = "yes" ; then
      for ax_dependency in $6 ; do
        ax_dependency_have="HAVE_PKG_${ax_dependency}"
        if test x"${!ax_dependency_have}" = "xyes"; then
          ax_dep_libs="PKG_${ax_dependency}_LIBS"
          ax_libs="${ax_libs} ${!ax_dep_libs}"
        else
          AC_MSG_RESULT([no])
          AC_MSG_NOTICE([warning: unable to build requested module $1 (needs ${ax_dependency})!])
          ax_module_enable=no;
          break;
        fi
      done
    fi

    if test "$ax_module_enable" = "yes" ; then
      # Check for optional dependencies
      for ax_dependency in $7 ; do
        ax_dependency_have="HAVE_PKG_${ax_dependency}"
        if test x"${!ax_dependency_have}" = "xyes"; then
          ax_dep_libs="PKG_${ax_dependency}_LIBS"
          ax_libs="${ax_libs} ${!ax_dep_libs}"
        fi
      done

      # Set up the variables
      MODULE_$1_LIBS=$ax_libs
      PKG_$1_LIBS="$ax_libs \$(top_srcdir)/$2/$3"
      AC_MSG_RESULT([yes])
    fi
  
  else
    HAVE_PKG_$1_SRC=no
    ax_module_enable=no
    MODULE_$1_LIBS=
    PKG_$1_LIBS=
  fi

  AC_SUBST(MODULE_$1_LIBS)
  AC_SUBST(PKG_$1_LIBS)

  HAVE_PKG_$1=${ax_module_enable}
  MAKE_MODULE_$1=${ax_module_enable}
  AC_SUBST(MAKE_MODULE_$1)

  if test "${HAVE_PKG_$1}" = "yes" ; then
    ax_have_pkg_bool=1
  else
    ax_have_pkg_bool=0
  fi
  AC_DEFINE_UNQUOTED(HAVE_PKG_$1,
                     [$ax_have_pkg_bool],
                     [Define to 1 if the $1 module is available.])

  if test "$ENABLE_VERBOSE" = "yes" && test "$HAVE_PKG_$1_SRC" = "yes" ; then
    AC_MSG_NOTICE(MAKE_MODULE_$1 = ${MAKE_MODULE_$1})
    AC_MSG_NOTICE(HAVE_PKG_$1 = ${HAVE_PKG_$1})
    AC_MSG_NOTICE(MODULE_$1_LIBS = ${MODULE_$1_LIBS})
    AC_MSG_NOTICE(PKG_$1_LIBS = ${PKG_$1_LIBS})
  fi

  #  We're putting these in configure.ac manually by now, for 
  #  backwards compatability with older versions of automake.
  #  AM_CONDITIONAL([MAKE_MODULE_$1], [test "$MAKE_MODULE_$1" = "yes"])
])

# This checks for various introspection functions
AC_DEFUN([AX_CHECK_INTROSPECTION],
[
    AC_CHECK_HEADERS(execinfo.h cxxabi.h typeinfo dlfcn.h)

    AC_CHECK_FUNCS([__cxa_current_exception_type __cxa_demangle backtrace])

    AC_SEARCH_LIBS(dladdr, dl, [AC_DEFINE(HAVE_DLADDR,1,Define if you have dladdr())])
])

AC_DEFUN([AC_PROG_SWIG],[
        AC_PATH_PROG([SWIG],[swig])
        if test -z "$SWIG" ; then
                AC_MSG_WARN([cannot find 'swig' program. You should look at http://www.swig.org])
                SWIG='echo "Error: SWIG is not installed. You should look at http://www.swig.org" ; false'
        elif test -n "$1" ; then
                AC_MSG_CHECKING([for SWIG version])
                [swig_version=`$SWIG -version 2>&1 | grep 'SWIG Version' | sed 's/.*\([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/g'`]
                AC_MSG_RESULT([$swig_version])
                if test -n "$swig_version" ; then
                        # Calculate the required version number components
                        [required=$1]
                        [required_major=`echo $required | sed 's/[^0-9].*//'`]
                        if test -z "$required_major" ; then
                                [required_major=0]
                        fi
                        [required=`echo $required | sed 's/[0-9]*[^0-9]//'`]
                        [required_minor=`echo $required | sed 's/[^0-9].*//'`]
                        if test -z "$required_minor" ; then
                                [required_minor=0]
                        fi
                        [required=`echo $required | sed 's/[0-9]*[^0-9]//'`]
                        [required_patch=`echo $required | sed 's/[^0-9].*//'`]
                        if test -z "$required_patch" ; then
                                [required_patch=0]
                        fi
                        # Calculate the available version number components
                        [available=$swig_version]
                        [available_major=`echo $available | sed 's/[^0-9].*//'`]
                        if test -z "$available_major" ; then
                                [available_major=0]
                        fi
                        [available=`echo $available | sed 's/[0-9]*[^0-9]//'`]
                        [available_minor=`echo $available | sed 's/[^0-9].*//'`]
                        if test -z "$available_minor" ; then
                                [available_minor=0]
                        fi
                        [available=`echo $available | sed 's/[0-9]*[^0-9]//'`]
                        [available_patch=`echo $available | sed 's/[^0-9].*//'`]
                        if test -z "$available_patch" ; then
                                [available_patch=0]
                        fi
                        if test $available_major -ne $required_major \
                                -o $available_minor -ne $required_minor \
                                -o $available_patch -lt $required_patch ; then
                                AC_MSG_WARN([SWIG version >= $1 is required.  You have $swig_version.  You should look at http://www.swig.org])
                                SWIG='echo "Error: SWIG version >= $1 is required.  You have '"$swig_version"'.  You should look at http://www.swig.org" ; false'
                        else
                                AC_MSG_NOTICE([SWIG executable is '$SWIG'])
                                SWIG_LIB=`$SWIG -swiglib`
                                AC_MSG_NOTICE([SWIG library directory is '$SWIG_LIB'])
                        fi
                else
                        AC_MSG_WARN([cannot determine SWIG version])
                        SWIG='echo "Error: Cannot determine SWIG version.  You should look at http://www.swig.org" ; false'
                fi
        fi
        AC_SUBST([SWIG_LIB])
])

AC_DEFUN([SWIG_ENABLE_CXX],[
        AC_REQUIRE([AC_PROG_SWIG])
        AC_REQUIRE([AC_PROG_CXX])
        SWIG="$SWIG -c++"
])

AC_DEFUN([SWIG_PYTHON],[
        AC_REQUIRE([AC_PROG_SWIG])
        AC_REQUIRE([AC_PYTHON_DEVEL])
        test "x$1" != "xno" || swig_shadow=" -noproxy"
        AC_SUBST([SWIG_PYTHON_OPT],[-python$swig_shadow])
        AC_SUBST([SWIG_PYTHON_CPPFLAGS],[$PYTHON_CPPFLAGS])
])

AC_DEFUN([AC_PYTHON_DEVEL],[
	#
	# Allow the use of a (user set) custom python version
	#
	AC_ARG_VAR([PYTHON_VERSION],[The installed Python
		version to use, for example '2.3'. This string
		will be appended to the Python interpreter
		canonical name.])

	AC_PATH_PROG([PYTHON],[python[$PYTHON_VERSION]])
	if test -z "$PYTHON"; then
	   AC_MSG_ERROR([Cannot find python$PYTHON_VERSION in your system path])
	   PYTHON_VERSION=""
	fi

	#
	# Check for a version of Python >= 2.1.0
	#
	AC_MSG_CHECKING([for a version of Python >= '2.1.0'])
	ac_supports_python_ver=`$PYTHON -c "import sys, string; \
		ver = string.split(sys.version)[[0]]; \
		print ver >= '2.1.0'"`
	if test "$ac_supports_python_ver" != "True"; then
		if test -z "$PYTHON_NOVERSIONCHECK"; then
			AC_MSG_RESULT([no])
			AC_MSG_FAILURE([
This version of the AC@&t@_PYTHON_DEVEL macro
doesn't work properly with versions of Python before
2.1.0. You may need to re-run configure, setting the
variables PYTHON_CPPFLAGS, PYTHON_LDFLAGS, PYTHON_SITE_PKG,
PYTHON_EXTRA_LIBS and PYTHON_EXTRA_LDFLAGS by hand.
Moreover, to disable this check, set PYTHON_NOVERSIONCHECK
to something else than an empty string.
])
		else
			AC_MSG_RESULT([skip at user request])
		fi
	else
		AC_MSG_RESULT([yes])
	fi

	#
	# if the macro parameter ``version'' is set, honour it
	#
	if test -n "$1"; then
		AC_MSG_CHECKING([for a version of Python $1])
		ac_supports_python_ver=`$PYTHON -c "import sys, string; \
			ver = string.split(sys.version)[[0]]; \
			print ver $1"`
		if test "$ac_supports_python_ver" = "True"; then
	   	   AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
			AC_MSG_ERROR([this package requires Python $1.
If you have it installed, but it isn't the default Python
interpreter in your system path, please pass the PYTHON_VERSION
variable to configure. See ``configure --help'' for reference.
])
			PYTHON_VERSION=""
		fi
	fi

	#
	# Check if you have distutils, else fail
	#
	AC_MSG_CHECKING([for the distutils Python package])
	ac_distutils_result=`$PYTHON -c "import distutils" 2>&1`
	if test -z "$ac_distutils_result"; then
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
		AC_MSG_ERROR([cannot import Python module "distutils".
Please check your Python installation. The error was:
$ac_distutils_result])
		PYTHON_VERSION=""
	fi

	#
	# Check for Python include path
	#
	AC_MSG_CHECKING([for Python include path])
	if test -z "$PYTHON_CPPFLAGS"; then
		python_path=`$PYTHON -c "import distutils.sysconfig; \
           		print distutils.sysconfig.get_python_inc();"`
		if test -n "${python_path}"; then
		   	python_path="-I$python_path"
		fi
		PYTHON_CPPFLAGS=$python_path
	fi
	AC_MSG_RESULT([$PYTHON_CPPFLAGS])
	AC_SUBST([PYTHON_CPPFLAGS])

	#
	# Check for Python library path
	#
	AC_MSG_CHECKING([for Python library path])
	if test -z "$PYTHON_LDFLAGS"; then
		# (makes two attempts to ensure we've got a version number
		# from the interpreter)
		py_version=`$PYTHON -c "from distutils.sysconfig import *; \
			from string import join; \
			print join(get_config_vars('VERSION'))"`
		if test "$py_version" == "[None]"; then
			if test -n "$PYTHON_VERSION"; then
				py_version=$PYTHON_VERSION
			else
				py_version=`$PYTHON -c "import sys; \
					print sys.version[[:3]]"`
			fi
		fi

		PYTHON_LDFLAGS=`$PYTHON -c "from distutils.sysconfig import *; \
			from string import join; \
			print '-L' + get_python_lib(0,1), \
		      	'-lpython';"`$py_version
	fi
	AC_MSG_RESULT([$PYTHON_LDFLAGS])
	AC_SUBST([PYTHON_LDFLAGS])

	#
	# Check for site packages
	#
	AC_MSG_CHECKING([for Python site-packages path])
	if test -z "$PYTHON_SITE_PKG"; then
		PYTHON_SITE_PKG=`$PYTHON -c "import distutils.sysconfig; \
		        print distutils.sysconfig.get_python_lib(0,0);"`
	fi
	AC_MSG_RESULT([$PYTHON_SITE_PKG])
	AC_SUBST([PYTHON_SITE_PKG])

	#
	# libraries which must be linked in when embedding
	#
	AC_MSG_CHECKING(python extra libraries)
	if test -z "$PYTHON_EXTRA_LIBS"; then
	   PYTHON_EXTRA_LIBS=`$PYTHON -c "import distutils.sysconfig; \
                conf = distutils.sysconfig.get_config_var; \
                print conf('LOCALMODLIBS'), conf('LIBS')"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LIBS])
	AC_SUBST(PYTHON_EXTRA_LIBS)

	#
	# linking flags needed when embedding
	#
	AC_MSG_CHECKING(python extra linking flags)
	if test -z "$PYTHON_EXTRA_LDFLAGS"; then
		PYTHON_EXTRA_LDFLAGS=`$PYTHON -c "import distutils.sysconfig; \
			conf = distutils.sysconfig.get_config_var; \
			print conf('LINKFORSHARED')"`
	fi
	AC_MSG_RESULT([$PYTHON_EXTRA_LDFLAGS])
	AC_SUBST(PYTHON_EXTRA_LDFLAGS)

	#
	# final check to see if everything compiles alright
	#
	AC_MSG_CHECKING([consistency of all components of python development environment])
	AC_LANG_PUSH([C])
	# save current global flags
	LIBS="$ac_save_LIBS $PYTHON_LDFLAGS"
	CPPFLAGS="$ac_save_CPPFLAGS $PYTHON_CPPFLAGS"
	AC_TRY_LINK([
		#include <Python.h>
	],[
		Py_Initialize();
	],[pythonexists=yes],[pythonexists=no])

	AC_MSG_RESULT([$pythonexists])

        if test ! "$pythonexists" = "yes"; then
	   AC_MSG_ERROR([
  Could not link test program to Python. Maybe the main Python library has been
  installed in some non-standard library path. If so, pass it to configure,
  via the LDFLAGS environment variable.
  Example: ./configure LDFLAGS="-L/usr/non-standard-path/python/lib"
  ============================================================================
   ERROR!
   You probably have to install the development version of the Python package
   for your distribution.  The exact name of this package varies among them.
  ============================================================================
	   ])
	  PYTHON_VERSION=""
	fi
	AC_LANG_POP
	# turn back to default flags
	CPPFLAGS="$ac_save_CPPFLAGS"
	LIBS="$ac_save_LIBS"

	#
	# all done!
	#
])
