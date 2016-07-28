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

if MAKE_MODULE_STEREO

TestAlgorithms_SOURCES    = TestAlgorithms.cxx
TestCorrelate_SOURCES     = TestCorrelate.cxx
TestCorrelationView_SOURCES = TestCorrelationView.cxx
TestCorrelation_SOURCES   = TestCorrelation.cxx
TestCostFunctions_SOURCES = TestCostFunctions.cxx
TestDisparity_SOURCES     = TestDisparity.cxx
TestPreFilter_SOURCES     = TestPreFilter.cxx
TestPyramidCorrelationView_SOURCES = TestPyramidCorrelationView.cxx
TestStereoModel_SOURCES   = TestStereoModel.cxx
TestSubPixel_SOURCES      = TestSubPixel.cxx
TestSGM_SOURCES = TestSGM.cxx

TESTS = \
	TestAlgorithms \
	TestCorrelate \
	TestCorrelationView \
	TestCorrelation \
	TestCostFunctions \
	TestDisparity \
	TestPreFilter \
	TestPyramidCorrelationView \
	TestStereoModel \
	TestSubPixel \
	TestSGM

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_STEREO_LIBS@

check_PROGRAMS = $(TESTS)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
