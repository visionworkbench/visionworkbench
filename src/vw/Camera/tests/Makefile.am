# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CAMERA

nodist_TestCAHVORModel_SOURCES = TestCAHVORModel.cxx
nodist_TestCAHVModel_SOURCES = TestCAHVModel.cxx
nodist_TestLinearPushbroomModel_SOURCES = TestLinearPushbroomModel.cxx
nodist_TestPinholeModel_SOURCES = TestPinholeModel.cxx
nodist_TestPinholeModelCalibrate_SOURCES = TestPinholeModelCalibrate.cxx
nodist_TestBundleAdjust_SOURCES = TestBundleAdjust.cxx
nodist_TestExifData_SOURCES = TestExifData.cxx
#nodist_TestLensDistortion_SOURCES = TestLensDistortion.cxx
#nodist_TestCameraTransform_SOURCES = TestCameraTransform.cxx

TESTS = TestCAHVORModel TestLinearPushbroomModel TestPinholeModel TestPinholeModelCalibrate \
	TestBundleAdjust TestExifData TestCAHVModel # TestLensDistortion TestCameraTransform

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_CAMERA_LIBS@ @PKG_VW_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(TESTS:%=%.h) test.jpg
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
