# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
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
TestLinearPushbroomModel_SOURCES  = TestLinearPushbroomModel.cxx
TestPinholeModel_SOURCES          = TestPinholeModel.cxx
TestPinholeModelCalibrate_SOURCES = TestPinholeModelCalibrate.cxx
TestAdjustedCamera_SOURCES        = TestAdjustedCamera.cxx

#TestLensDistortion_SOURCES       = TestLensDistortion.oldtest
#TestCameraTransform_SOURCES      = TestCameraTransform.oldtest

TESTS = TestCAHVModel TestCAHVORModel TestCAHVOREModel        \
        TestCameraGeometry TestExifData                       \
        TestLinearPushbroomModel TestPinholeModel             \
        TestPinholeModelCalibrate TestAdjustedCamera

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
