# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_PLATE

TestBlobManager_SOURCES       = TestBlobManager.cxx
TestBlobIO_SOURCES            = TestBlobIO.cxx
TestIndexPage_SOURCES         = TestIndexPage.cxx
TestLocalIndex_SOURCES        = TestLocalIndex.cxx
TestAmqp_SOURCES              = TestAmqp.cxx
TestTileManipulation_SOURCES  = TestTileManipulation.cxx
#TestPlateIndexServer_SOURCES = TestPlateIndexServer.cxx

TESTS = TestBlobManager TestBlobIO TestLocalIndex \
        TestIndexPage TestAmqp TestTileManipulation

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_IMAGE_LIBS@ @PKG_MATH_LIBS@ @PKG_CORE_LIBS@ @PKG_PLATE_LIBS@

check_PROGRAMS = $(TESTS)
CLEANFILES =

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
