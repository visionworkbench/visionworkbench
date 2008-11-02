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
