# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_FILEIO

TestDiskImageResource_SOURCES = TestDiskImageResource.cxx
TestDiskImageView_SOURCES     = TestDiskImageView.cxx
TestEndianness_SOURCES        = TestEndianness.cxx
TestBlockFileIO_SOURCES       = TestBlockFileIO.cxx
TestGDALFeatures_SOURCES      = TestGDALFeatures.cxx
TestMemoryImageResource_SOURCES = TestMemoryImageResource.cxx

TESTS = \
  TestBlockFileIO \
  TestDiskImageResource \
  TestDiskImageView \
  TestMemoryImageResource \
  TestEndianness \
  TestGDALFeatures

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_FILEIO_LIBS@

check_PROGRAMS = $(TESTS)

EXTRA_DIST = rgb2x2.png rgb2x2.tif rgb2x2.jpg png16.png mural.png mural.jpg

CLEANFILES = tmp.png tmp.tif rwtest.* test-png16.png cropped.mural.* mural.tif nodata.tif

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
