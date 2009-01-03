dnl AX_LOAD_DEPS(<pkg>, <deps>[, <missing>])
dnl copies deps into the pkg CPPFLAGS/LIBS
dnl missing could be a shell var name in which to store
dnl the names of the missing deps
AC_DEFUN([AX_LOAD_DEPS],
 [AS_VAR_PUSHDEF([pkg_c], [PKG_]m4_toupper([[$1]])[_CPPFLAGS])
  AS_VAR_PUSHDEF([pkg_l], [PKG_]m4_toupper([[$1]])[_LIBS])
  AS_VAR_PUSHDEF([missing], m4_default([$3], [ax_load_deps_missing_deps]))

  m4_foreach_w(dep, m4_toupper([[$2]]),
   [AS_IF([test x"$HAVE_PKG_]dep[" != "xyes"], [missing="$missing dep"],
      [pkg_c="$pkg_c $PKG_]dep[_CPPFLAGS"
       pkg_l="$pkg_l $PKG_]dep[_LIBS"])])

  AS_VAR_POPDEF([missing])
  AS_VAR_POPDEF([pkg_l])
  AS_VAR_POPDEF([pkg_c])
])
