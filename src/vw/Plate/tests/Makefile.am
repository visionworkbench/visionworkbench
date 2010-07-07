# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

check_SCRIPTS  =

if MAKE_MODULE_PLATE

TestBlobManager_SOURCES       = TestBlobManager.cxx
TestBlobIO_SOURCES            = TestBlobIO.cxx
TestIndexPage_SOURCES         = TestIndexPage.cxx
TestLocalIndex_SOURCES        = TestLocalIndex.cxx
TestAmqp_SOURCES              = TestAmqp.cxx
TestTileManipulation_SOURCES  = TestTileManipulation.cxx
TestModPlate_SOURCES          = TestModPlate.cxx
TestHTTPUtils_SOURCES         = TestHTTPUtils.cxx

if HAVE_PYTHON
check_SCRIPTS  += TestModPlateRegression.py
endif

check_PROGRAMS = TestBlobManager TestBlobIO TestLocalIndex TestIndexPage TestAmqp \
				 TestTileManipulation TestModPlate TestHTTPUtils


if MAKE_MODPLATE
MORE_FLAGS = -DVW_HAVE_APACHE=1
endif

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ @APXS_CFLAGS@ $(MORE_FLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @APXS_LDFLAGS@ @PKG_IMAGE_LIBS@ @PKG_MATH_LIBS@ @PKG_CORE_LIBS@ @PKG_PLATE_LIBS@

TESTS = $(check_PROGRAMS) $(check_SCRIPTS)
CLEANFILES =

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
