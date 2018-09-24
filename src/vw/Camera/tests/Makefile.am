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

if MAKE_MODULE_CAMERA

TestCAHVModel_SOURCES             = TestCAHVModel.cxx
TestCAHVORModel_SOURCES           = TestCAHVORModel.cxx
TestCAHVOREModel_SOURCES          = TestCAHVOREModel.cxx
TestCameraGeometry_SOURCES        = TestCameraGeometry.cxx
TestExifData_SOURCES              = TestExifData.cxx
TestExtrinsics_SOURCES            = TestExtrinsics.cxx
TestPinholeModel_SOURCES          = TestPinholeModel.cxx
TestAdjustedCamera_SOURCES        = TestAdjustedCamera.cxx
TestLinescanModel_SOURCES         = TestLinescanModel.cxx

#TestLensDistortion_SOURCES       = TestLensDistortion.oldtest
#TestCameraTransform_SOURCES      = TestCameraTransform.oldtest

TESTS = TestCAHVModel TestCAHVORModel TestCAHVOREModel  \
        TestCameraGeometry TestExifData TestExtrinsics  \
        TestLinescanModel TestPinholeModel               \
        TestAdjustedCamera

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_CAMERA_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = test.jpg cahvor.txt

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
