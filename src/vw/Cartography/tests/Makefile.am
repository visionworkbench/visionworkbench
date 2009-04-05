# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CARTOGRAPHY

nodist_TestGeoReference_SOURCES = TestGeoReference.cxx
nodist_TestGeoTransform_SOURCES = TestGeoTransform.cxx
nodist_TestPointImageManipulation_SOURCES = TestPointImageManipulation.cxx
nodist_TestToastTransform_SOURCES = TestToastTransform.cxx

TESTS = TestGeoReference TestGeoTransform TestPointImageManipulation	\
	TestToastTransform

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_CARTOGRAPHY_LIBS@ @PKG_VW_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(TESTS:%=%.h)
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
