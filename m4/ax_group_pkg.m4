dnl AX_GROUP_PKG(<pkg> [, <dep1> .. <depn>])
dnl A stripped-down AX_PKG that copies the vars from one pkg to another without
dnl re-checking. As long as all deps were found, pkg will be marked found,
dnl too. note: no commas after the first one separating pkg from dep1!
dnl
dnl Also, to help other macros: if there are no deps, pkg is marked NOT FOUND.

AC_DEFUN([AX_GROUP_PKG],
[AC_MSG_CHECKING([for package $1])
  m4_define([pkg], m4_toupper([[$1]]))
  AS_VAR_PUSHDEF([have_pkg], [HAVE_PKG_]pkg)
  AS_VAR_PUSHDEF([pkg_c], [PKG_]pkg[_CPPFLAGS])
  AS_VAR_PUSHDEF([pkg_l], [PKG_]pkg[_LIBS])
  AS_VAR_PUSHDEF([missing], [ax_group_pkg_missing_deps])
  AS_VAR_PUSHDEF([bool], [ax_group_pkg_have_pkg_bool])

  AC_ARG_WITH(m4_tolower([[$1]]),
    AC_HELP_STRING([--with-]m4_tolower([[$1]]), [enable searching for the $1 package @<:@auto@:>@]),
    [ have_pkg=$withval ]
  )

  m4_if([$2], [], [have_pkg=[no_deps]],
  [m4_foreach_w(dep, m4_toupper([[$2]]),
    [AS_IF([test x"$HAVE_PKG_]dep[" != "xyes"], [missing="$missing dep"],
      [pkg_c="$pkg_c $PKG_]dep[_CPPFLAGS"
       pkg_l="$pkg_l $PKG_]dep[_LIBS"])])])

  AS_IF(
    [test x"$have_pkg" = "xno"],      [AS_VAR_SET([bool], 0); AC_MSG_RESULT([no (disabled by user)])],
    [test x"$have_pkg" = "xno_deps"], [AS_VAR_SET([bool], 0); AC_MSG_RESULT([no])],
    [test -z "$missing"],             [AS_VAR_SET([bool], 1); AC_MSG_RESULT([yes])],
    [bool=0; AC_MSG_RESULT([no ([missing] $missing)]) ])

  AS_IF( [test x"$bool" = "x1"], [have_pkg=yes], [pkg_c=""; pkg_l=""; have_pkg=no])

  AC_DEFINE_UNQUOTED(have_pkg, $bool, [Define to 1 if the pkg package is available])
  AC_SUBST(pkg_c)
  AC_SUBST(pkg_l)
  AC_SUBST(have_pkg)

  AS_VAR_POPDEF([bool])
  AS_VAR_POPDEF([missing])
  AS_VAR_POPDEF([pkg_l])
  AS_VAR_POPDEF([pkg_c])
  AS_VAR_POPDEF([have_pkg])
  m4_undefine([pkg])
])
