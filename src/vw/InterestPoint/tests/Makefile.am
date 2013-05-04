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


########################################################################
# sources
########################################################################

if MAKE_MODULE_INTERESTPOINT

TestMatcher_SOURCES   = TestMatcher.cxx
TestIntegral_SOURCES  = TestIntegral.cxx
TestBoxFilter_SOURCES = TestBoxFilter.cxx
TestInterestData_SOURCES = TestInterestData.cxx

TESTS = TestMatcher TestIntegral TestBoxFilter TestInterestData

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_INTERESTPOINT_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = sub.png noisy_gradient_60.png
#CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
