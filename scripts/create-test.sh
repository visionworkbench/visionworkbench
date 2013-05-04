#!/bin/sh
# __BEGIN_LICENSE__
#  Copyright (c) 2006-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NASA Vision Workbench is licensed under the Apache License,
#  Version 2.0 (the "License"); you may not use this file except in
#  compliance with the License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__


#set -x

begins_with() {
    if test x"${1#$2}" = "x$1"; then
        return 1 #false!
    else
        return 0 #true!
    fi
}

ends_with() {
    if test x"${1%$2}" = "x$1"; then
        return 1 #false!
    else
        return 0 #true!
    fi
}

if test "$#" -lt "1"; then
    echo "Usage: $0 <test-suite-name> [<module-name>]"
    echo "module-name is inferred from pwd if possible"
    exit 1
fi

test_name="$1"

wd="`pwd`"
srcdir="$(dirname $(cd ${0%/*}/../ && echo $PWD/${0##*/}))"
rel="${wd##$srcdir}"

if test -z "$2"; then
    if ends_with $rel /tests; then
        module_name="${rel%/tests}"
        module_name="${module_name##*/}"
        test_path="$wd"
    elif begins_with $rel /src/vw/; then
        module_name="${rel#/src/vw/}"
        module_name="${module_name%%/*}"
        test_path="$srcdir/src/vw/$module_name/tests"
    else
        echo "Couldn't infer module name from current directory"
        echo "$wd"
        exit 1
    fi
else
    module_name="$2"
    test_path="$srcdir/src/vw/$module_name/tests"
fi

if test -d "$test_path"; then :; else
    echo "Directory $test_path does not exist. Please create it first."
    exit 1
fi

full_test="$test_path/Test${test_name}.h"

if test -e "$full_test"; then
    echo "Refusing to overwrite test"
    echo "$full_test"
    exit 1
fi

echo "creating ${full_test}"
cat <<EOF > ${full_test}
#include <cxxtest/TestSuite.h>
#include <vw/${module_name}.h>

class ${test_name} : public CxxTest::TestSuite
{public:

   void testOne()
   {
      TS_ASSERT_EQUALS( 1 + 1, 2 );
      TS_ASSERT_EQUALS( 1 + 2, 2 );
   }

};
EOF


# rel = */x/tests -> module[x], create test in pwd (since src/vw/tests matches)
# rel = src/vw/x  -> module[x], create test in $srcdir/src/vw/x/tests
# else module-name is required
