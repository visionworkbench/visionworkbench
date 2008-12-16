# These are backports from autoconf 2.60.
# They can be removed when 2.60 is everywhere.

# _AS_VERSION_COMPARE_PREPARE
# ---------------------------
# Output variables for comparing version numbers.
m4_ifdef([_AS_VERSION_COMPARE_PREPARE], [], [
  AC_DEFUN([_AS_VERSION_COMPARE_PREPARE],
  [[as_awk_strverscmp='
    # Use only awk features that work with 7th edition Unix awk (1978).
    # My, what an old awk you have, Mr. Solaris!
    END {
      while (length(v1) || length(v2)) {
        # Set d1 to be the next thing to compare from v1, and likewise for d2.
        # Normally this is a single character, but if v1 and v2 contain digits,
        # compare them as integers and fractions as strverscmp does.
        if (v1 ~ /^[0-9]/ && v2 ~ /^[0-9]/) {
     # Split v1 and v2 into their leading digit string components d1 and d2,
     # and advance v1 and v2 past the leading digit strings.
     for (len1 = 1; substr(v1, len1 + 1) ~ /^[0-9]/; len1++) continue
     for (len2 = 1; substr(v2, len2 + 1) ~ /^[0-9]/; len2++) continue
     d1 = substr(v1, 1, len1); v1 = substr(v1, len1 + 1)
     d2 = substr(v2, 1, len2); v2 = substr(v2, len2 + 1)
     if (d1 ~ /^0/) {
       if (d2 ~ /^0/) {
         # Compare two fractions.
         while (d1 ~ /^0/ && d2 ~ /^0/) {
           d1 = substr(d1, 2); len1--
           d2 = substr(d2, 2); len2--
         }
         if (len1 != len2 && ! (len1 && len2 && substr(d1, 1, 1) == substr(d2, 1, 1))) {
           # The two components differ in length, and the common prefix
           # contains only leading zeros.  Consider the longer to be less.
           d1 = -len1
           d2 = -len2
         } else {
           # Otherwise, compare as strings.
           d1 = "x" d1
           d2 = "x" d2
         }
       } else {
         # A fraction is less than an integer.
         exit 1
       }
     } else {
       if (d2 ~ /^0/) {
         # An integer is greater than a fraction.
         exit 2
       } else {
         # Compare two integers.
         d1 += 0
         d2 += 0
       }
     }
        } else {
     # The normal case, without worrying about digits.
     if (v1 == "") d1 = v1; else { d1 = substr(v1, 1, 1); v1 = substr(v1,2) }
     if (v2 == "") d2 = v2; else { d2 = substr(v2, 1, 1); v2 = substr(v2,2) }
        }
        if (d1 < d2) exit 1
        if (d1 > d2) exit 2
      }
    }
']])])

# AS_VERSION_COMPARE(VERSION-1, VERSION-2,
#                    [ACTION-IF-LESS], [ACTION-IF-EQUAL], [ACTION-IF-GREATER])
# -----------------------------------------------------------------------------
# Compare two strings possibly containing shell variables as version strings.
m4_ifdef([AS_VERSION_COMPARE], [], [
  AC_DEFUN([AS_VERSION_COMPARE],
  [AS_REQUIRE([_$0_PREPARE])dnl
  as_arg_v1=$1
  as_arg_v2=$2
  dnl This usage is portable even to ancient awk,
  dnl so don't worry about finding a "nice" awk version.
  awk "$as_awk_strverscmp" v1="$as_arg_v1" v2="$as_arg_v2" /dev/null
  case $? in
  1) $3;;
  0) $4;;
  2) $5;;
  esac[]dnl
])])
