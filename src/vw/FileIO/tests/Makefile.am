# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_FILEIO

nodist_TestDiskImageResource_SOURCES = TestDiskImageResource.cxx
nodist_TestDiskImageView_SOURCES = TestDiskImageView.cxx
nodist_TestEndianness_SOURCES = TestEndianness.cxx

TESTS = TestDiskImageResource TestDiskImageView TestEndianness

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_FILEIO_LIBS@ @PKG_IMAGE_LIBS@ @PKG_MATH_LIBS@ @PKG_CORE_LIBS@

check_PROGRAMS = $(TESTS)

EXTRA_DIST = $(TESTS:%=%.h) rgb2x2.png rgb2x2.tif rgb2x2.jpg png16.png

CLEANFILES = $(TESTS:%=%.cxx) tmp.png tmp.tif rwtest.* test-png16.png

include $(top_srcdir)/config/rules.mak
