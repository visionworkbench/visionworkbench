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

  AC_LANG_SAVE
  AC_LANG(C++)

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

    # Otherwise we look for a path that contains the needed headers and libraries
    else

      if test "x$ENABLE_VERBOSE" = "yes"; then
	AC_MSG_RESULT([searching...])
      fi

      HAVE_PKG_$1=no

      ax_pkg_old_libs=$LIBS
      LIBS=$PKG_$1_LIBS $LIBS
      for path in none $PKG_PATHS; do
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
	  VW_LDFLAGS="-L$path/lib $VW_LDFLAGS"
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

  AC_LANG_RESTORE
])


# Usage: AX_PKG_BOOST
AC_DEFUN([AX_PKG_BOOST],
[
  AC_MSG_CHECKING(for package BOOST)
  if test "$ENABLE_VERBOSE" = "yes"; then
    AC_MSG_RESULT([])
  fi

  AC_LANG_SAVE
  AC_LANG(C++)

  # Skip testing if the user has overridden
  if test -z ${HAVE_PKG_BOOST}; then

    PKG_BOOST_LIBS=
    HAVE_PKG_BOOST=no

    for ax_boost_base_path in $PKG_PATHS; do
      # First look for a system-style installation
      if test "$ENABLE_VERBOSE" = "yes"; then
        AC_MSG_CHECKING([for system-style boost in ${ax_boost_base_path}])
      fi
      if test -d "${ax_boost_base_path}/include/boost" ; then
        PKG_BOOST_INCDIR="${ax_boost_base_path}/include"
        PKG_BOOST_LIBDIR="${ax_boost_base_path}/lib"
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

  AC_LANG_RESTORE
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
dnl Thus we have this function, which is more flexible. It takes an 
dnl additional parameter at the end, which is a "semi-suffix" -- text that 
dnl appears at either the end of the library filename (but before the .so),
dnl or somewhere in the middle (in case there's another suffix).
# Usage: AX_PKG_BOOST_LIB(<name>, <libraries>, <header>, <semi-suffix>)
AC_DEFUN([AX_PKG_BOOST_LIB],
[
  AC_MSG_CHECKING(for package BOOST_$1)
  if test "$ENABLE_VERBOSE" = "yes"; then
    AC_MSG_RESULT([])
  fi

  AC_LANG_SAVE
  AC_LANG(C++)

  # Skip testing if the user has overridden
  if test -z ${HAVE_PKG_BOOST_$1}; then

    HAVE_PKG_BOOST_$1=no

    # Check for general Boost presence
    if test "x${HAVE_PKG_BOOST}" = "xyes" ; then
      # Check for required headers
      AX_FIND_FILES([$3],[${PKG_BOOST_INCDIR}])
      if test ! -z "$ax_find_files_path" ; then
        # Check for required libraries with no suffix
        AX_FIND_FILES([`echo $2 | sed 's/-l\([[^[:space:]]]*\)/lib\1$4.*/g'`],[$PKG_BOOST_LIBDIR])
        if test ! -z "$ax_find_files_path" ; then
          HAVE_PKG_BOOST_$1="yes"
          PKG_BOOST_$1_LIBS="$2"
        else
          # Check for required libraries with some suffix. We have to check
          # for both a suffix before $4 (pre-suffix) and a suffix 
          # after (for example) the -mt (post-suffix), because boost likes 
          # to stick the name of the compiler before the -mt. Extremely
          # annoying.

          ax_pkg_boost_lib=`echo $2 | awk '{print [$]1}' | sed 's/-l\([[^[:space:]-]]*\).*/lib\1/g'`
          ax_pkg_boost_file=`ls ${PKG_BOOST_LIBDIR}/${ax_pkg_boost_lib}-*$4* | head -n 1 | sed "s,^${PKG_BOOST_LIBDIR}/\(.*\),\1,"`

          # The pre-suffix.
          ax_pkg_boost_presuffix=`echo ${ax_pkg_boost_file} | sed "s/${ax_pkg_boost_lib}\(.*\)-$4.*/\1/"`
          if ! test -z $ax_pkg_boost_presuffix; then
            # Add - to the end of the presuffix. We might want to make 
            # this a little more extensible (allow other characters 
            # besides -).
            ax_pkg_boost_presuffix=${ax_pkg_boost_presuffix}-
          fi

          # The post-suffix.
          ax_pkg_boost_postsuffix=`echo ${ax_pkg_boost_file} | sed "s/${ax_pkg_boost_lib}.*-$4-\([[^.]]*\).*/\1/"`
          if ! test -z $ax_pkg_boost_postsuffix; then
            # Add - to the start of the postsuffix. We might want to make 
            # this a little more extensible (allow other characters 
            # besides -).
            ax_pkg_boost_postsuffix=-${ax_pkg_boost_postsuffix}
          fi

          AX_FIND_FILES([`echo $2 | sed "s/-l\([[^[:space:]]]*\)/lib\1${ax_pkg_boost_presuffix}$4${ax_pkg_boost_postsuffix}.*/g"`],[$PKG_BOOST_LIBDIR])
          if test ! -z $ax_find_files_path ; then
            HAVE_PKG_BOOST_$1="yes"
            PKG_BOOST_$1_LIBS=`echo $2 | sed "s/[[^ ]]*/&${ax_pkg_boost_presuffix}$4${ax_pkg_boost_postsuffx}/g"`
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

  AC_LANG_RESTORE
])

dnl Usage: AX_PKG_PTHREADS
AC_DEFUN([AX_PKG_PTHREADS],
[
  AC_MSG_CHECKING(for package PTHREADS)
  if test "$ENABLE_VERBOSE" = "yes"; then
    AC_MSG_RESULT([])
  fi

  AC_LANG_SAVE
  AC_LANG_C
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

  AC_LANG_RESTORE

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
