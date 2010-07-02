# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CARTOGRAPHY

TestGeoReference_SOURCES           = TestGeoReference.cxx
TestGeoTransform_SOURCES           = TestGeoTransform.cxx
TestPointImageManipulation_SOURCES = TestPointImageManipulation.cxx
TestToastTransform_SOURCES         = TestToastTransform.cxx
TestCameraBBox_SOURCES             = TestCameraBBox.cxx

TESTS = TestGeoReference TestGeoTransform TestPointImageManipulation   \
        TestToastTransform TestCameraBBox

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_CARTOGRAPHY_LIBS@ @PKG_CAMERA_LIBS@

check_PROGRAMS = $(TESTS)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
