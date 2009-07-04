# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_FILEIO

TESTS = TestBlockFileIO TestExtTester

nodist_TestBlockFileIO_SOURCES = TestBlockFileIO.cxx
nodist_TestExtTester_SOURCES = TestExtTester.cxx

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS) @PKG_BOOST_FILESYSTEM_CPPFLAGS@
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_VW_LIBS@ @PKG_BOOST_FILESYSTEM_LIBS@

check_PROGRAMS = $(TESTS)

EXTRA_DIST = $(TESTS:%=%.h) mural.png mural.jpg

CLEANFILES = $(TESTS:%=%.cxx) cropped.mural.* mural.tif

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests_old.am
