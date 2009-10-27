# __BEGIN_LICENSE__
# Copyright (C) 2006-2009 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CAMERA

TestCAHVORModel_SOURCES           = TestCAHVORModel.oldtest
TestCAHVModel_SOURCES             = TestCAHVModel.oldtest
TestLinearPushbroomModel_SOURCES  = TestLinearPushbroomModel.oldtest
TestPinholeModel_SOURCES          = TestPinholeModel.oldtest
TestPinholeModelCalibrate_SOURCES = TestPinholeModelCalibrate.oldtest
TestBundleAdjust_SOURCES          = TestBundleAdjust.oldtest
TestExifData_SOURCES              = TestExifData.oldtest
#TestLensDistortion_SOURCES       = TestLensDistortion.oldtest
#TestCameraTransform_SOURCES      = TestCameraTransform.oldtest

TESTS = TestCAHVORModel TestLinearPushbroomModel TestPinholeModel TestPinholeModelCalibrate \
	TestBundleAdjust TestExifData TestCAHVModel # TestLensDistortion TestCameraTransform

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_CAMERA_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = test.jpg cahvor.txt
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests_old.am
