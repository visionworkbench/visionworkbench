# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CAMERA

TestBundleAdjustment_SOURCES      = TestBundleAdjustment.cxx
TestCAHVModel_SOURCES             = TestCAHVModel.cxx
TestCAHVORModel_SOURCES           = TestCAHVORModel.cxx
TestCameraGeometry_SOURCES        = TestCameraGeometry.cxx
TestExifData_SOURCES              = TestExifData.cxx
TestLinearPushbroomModel_SOURCES  = TestLinearPushbroomModel.cxx
TestPinholeModel_SOURCES          = TestPinholeModel.cxx
TestPinholeModelCalibrate_SOURCES = TestPinholeModelCalibrate.cxx

#TestLensDistortion_SOURCES       = TestLensDistortion.oldtest
#TestCameraTransform_SOURCES      = TestCameraTransform.oldtest

TESTS = TestBundleAdjustment TestCAHVModel TestCAHVORModel TestCameraGeometry \
        TestExifData TestLinearPushbroomModel TestPinholeModel                \
        TestPinholeModelCalibrate

#TESTS = TestCAHVORModel TestLinearPushbroomModel TestPinholeModel TestPinholeModelCalibrate \
#        TestExifData TestCAHVModel TestCameraGeometry TestBundleAdjustment TestLensDistortion TestCameraTransform

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_CAMERA_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = test.jpg cahvor.txt
#CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
