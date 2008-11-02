dnl usage: AX_MODULE_PYTHON(<default>)
AC_DEFUN([AX_MODULE_PYTHON],
[
  AC_ARG_ENABLE([module-python],
    AC_HELP_STRING([--enable-module-python], [enable the python bindings @<:@$1@:>@]),
    [ ENABLE_MODULE_PYTHON=$enableval ],
    [ if test x"$ENABLE_MODULE_PYTHON" = x""; then ENABLE_MODULE_PYTHON=`/bin/echo -n $1 | tr [A-Z] [a-z]` ; fi ]
  )

  AC_MSG_CHECKING([whether you want python bindings])
  ax_module_enable=$ENABLE_MODULE_PYTHON

  if test "$ax_module_enable" != "yes" ; then
    AC_MSG_RESULT([no (disabled)])
  else
    AC_MSG_RESULT([yes])

    AM_PATH_PYTHON([2.4])
    AC_PROG_SWIG([1.3.35])

    ax_module_enable=no

    if test -z "$PYTHON" || test -z "$SWIG"; then
      AC_MSG_RESULT([no])
      AC_MSG_NOTICE([warning: python bindings need swig and python])
    else
      SWIG_ENABLE_CXX
      SWIG_PYTHON
      AC_PYTHON_MODULE([numpy])

      if test "$HAVE_PYMOD_NUMPY" != "yes"; then
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([warning: python bindings need numpy])
      else
        AC_MSG_CHECKING([for numpy include path])
        numpy_include=`$PYTHON -c 'import numpy; print numpy.get_include();'`
        NUMPY_CPPFLAGS="-I$numpy_include"
        AC_MSG_RESULT([$NUMPY_CPPFLAGS])
        ax_module_enable=yes
      fi
    fi
  fi

  AC_MSG_CHECKING([whether to build python bindings])
  if test "$ax_module_enable" != "yes" ; then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
  fi

  AC_SUBST([NUMPY_CPPFLAGS])

  HAVE_PKG_PYTHON=${ax_module_enable}
  MAKE_MODULE_PYTHON=${ax_module_enable}
  AC_SUBST(MAKE_MODULE_PYTHON)

  if test "${HAVE_PKG_PYTHON}" = "yes" ; then
    ax_have_pkg_bool=1
  else
    ax_have_pkg_bool=0
  fi

  AC_DEFINE_UNQUOTED(HAVE_PKG_PYTHON,
                     [$ax_have_pkg_bool],
                     [Define to 1 if the PYTHON module is available.])

  if test "$ENABLE_VERBOSE" = "yes"; then
    AC_MSG_NOTICE(MAKE_MODULE_PYTHON = ${MAKE_MODULE_PYTHON})
    AC_MSG_NOTICE(HAVE_PKG_PYTHON = ${HAVE_PKG_PYTHON})
  fi
])


