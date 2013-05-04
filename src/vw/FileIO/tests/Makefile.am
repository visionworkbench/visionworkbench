# __BEGIN_LICENSE__
#  Copyright (c) 2006-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NASA Vision Workbench is licensed under the Apache License,
#  Version 2.0 (the "License"); you may not use this file except in
#  compliance with the License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
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
TestTemporaryFile_SOURCES    = TestTemporaryFile.cxx

TESTS = \
  TestBlockFileIO \
  TestDiskImageResource \
  TestDiskImageView \
  TestMemoryImageResource \
  TestTemporaryFile \
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

EXTRA_DIST = mural.jpg mural.png png16.png rgb2x2.jpg rgb2x2.png rgb2x2.tif rgb4x4_alpha.png rgb4x4_alpha.tif rgb4x4f_alpha.tif rgb4x4f_band.tif rgb4x4_halfalpha.png rgb4x4_halfalpha.tif rgb2x2.exr

CLEANFILES = tmp.png tmp.tif rwtest.* test-png16.png cropped.mural.* mural.tif nodata.tif

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
