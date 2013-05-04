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

if MAKE_MODULE_CORE

TestCache_SOURCES            = TestCache.cxx
TestCompoundTypes_SOURCES    = TestCompoundTypes.cxx
TestExceptions_SOURCES       = TestExceptions.cxx
TestFunctors_SOURCES         = TestFunctors.cxx
TestFundamentalTypes_SOURCES = TestFundamentalTypes.cxx
TestLog_SOURCES              = TestLog.cxx
TestSettings_SOURCES         = TestSettings.cxx
TestThreadPool_SOURCES       = TestThreadPool.cxx
TestThreadQueue_SOURCES      = TestThreadQueue.cxx
TestThread_SOURCES           = TestThread.cxx
TestTypeDeduction_SOURCES    = TestTypeDeduction.cxx

TESTS = \
  TestCache \
  TestCompoundTypes \
  TestExceptions \
  TestFunctors \
  TestFundamentalTypes \
  TestLog \
  TestSettings \
  TestThread \
  TestThreadPool \
  TestThreadQueue \
  TestTypeDeduction

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_CORE_LIBS@

check_PROGRAMS = $(TESTS)
CLEANFILES = log.txt

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
