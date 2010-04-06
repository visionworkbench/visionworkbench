# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_PLATE

TestBlobManager_SOURCES       = TestBlobManager.oldtest
TestBlobIO_SOURCES            = TestBlobIO.oldtest
TestBlobIterator_SOURCES      = TestBlobIterator.oldtest
TestIndexPage_SOURCES         = TestIndexPage.oldtest
TestLocalIndex_SOURCES        = TestLocalIndex.oldtest
TestAmqp_SOURCES              = TestAmqp.oldtest
TestTileManipulation_SOURCES  = TestTileManipulation.oldtest
#TestPlateIndexServer_SOURCES = TestPlateIndexServer.oldtest

TESTS = TestBlobManager TestBlobIO TestBlobIterator TestLocalIndex \
		TestIndexPage TestAmqp TestTileManipulation

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_IMAGE_LIBS@ @PKG_MATH_LIBS@ @PKG_CORE_LIBS@ @PKG_PLATE_LIBS@

check_PROGRAMS = $(TESTS)
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests_old.am
