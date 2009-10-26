# __BEGIN_LICENSE__
# Copyright (C) 2006-2009 United States Government as represented by
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
nodist_TestBlockFileIO_SOURCES = TestBlockFileIO.cxx
nodist_TestExtTester_SOURCES = TestExtTester.cxx

TESTS = TestDiskImageResource TestDiskImageView TestEndianness TestBlockFileIO TestExtTester

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_FILEIO_LIBS@

check_PROGRAMS = $(TESTS)

EXTRA_DIST = $(TESTS:%=%.h) rgb2x2.png rgb2x2.tif rgb2x2.jpg png16.png mural.png mural.jpg

CLEANFILES = $(TESTS:%=%.cxx) tmp.png tmp.tif rwtest.* test-png16.png cropped.mural.* mural.tif

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests_old.am
