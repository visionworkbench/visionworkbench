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
