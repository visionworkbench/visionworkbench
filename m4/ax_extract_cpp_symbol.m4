dnl __BEGIN_LICENSE__
dnl __END_LICENSE__

dnl Usage: AX_EXTRACT_CPP_SYMBOL([symbol], [headers], [ifyes], [ifno])
dnl         It the variable $output will contain the extracted value
dnl example: AX_EXTRACT_CPP_SYMBOL([BOOST_VERSION], [BOOST_VERSION="$output"], [BOOST_VERSION=nope], [#include <boost/version.hpp>])
AC_DEFUN([AX_EXTRACT_CPP_SYMBOL],
[
    AC_REQUIRE_CPP()dnl
    AC_CHECK_PROGS([SED], [sed gsed])
AC_LANG_CONFTEST(
    [AC_LANG_SOURCE([$2
#define __ac_extract_cpp_symbol_delimiter "__ac_extract_cpp_symbol_delimiter"
__ac_extract_cpp_symbol_delimiter $1 __ac_extract_cpp_symbol_delimiter])])

AS_VAR_PUSHDEF([output], [ax_extract_cpp_symbol_$1])dnl
if (eval "$ac_cpp conftest.$ac_ext") >conftest.out 2>&AC_FD_CC; then
    output="`${SED} -n -e 's/^.*"__ac_extract_cpp_symbol_delimiter" \(.*\) "__ac_extract_cpp_symbol_delimiter".*$/\1/p' conftest.out 2>/dev/null`"
    if test x"${output}" != x"$1"; then
        ifelse([$3], , :, [$3])
    ifelse([$4], , , [else
    $4
])
    fi
    ifelse([$4], , , [else
    $4
])dnl
fi
AS_VAR_POPDEF([output])
rm -f conftest*
])
