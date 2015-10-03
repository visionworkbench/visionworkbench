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

check_SCRIPTS  =

if MAKE_MODULE_PLATE

protocol_headers = TestRpcService.pb.h
protocol_sources = TestRpcService.pb.cc

BUILT_SOURCES = $(protocol_sources)

TestBlobIO_SOURCES            = TestBlobIO.cxx
TestBlobManager_SOURCES       = TestBlobManager.cxx
TestDatastore_SOURCES         = TestDatastore.cxx
TestHTTPUtils_SOURCES         = TestHTTPUtils.cxx
TestIndexPage_SOURCES         = TestIndexPage.cxx
TestLocalIndex_SOURCES        = TestLocalIndex.cxx
TestModPlate_SOURCES          = TestModPlate.cxx
TestPlateManager_SOURCES      = TestPlateManager.cxx
TestRpc_SOURCES               = TestRpc.cxx $(protocol_sources)
TestRpcChannel_SOURCES        = TestRpcChannel.cxx
TestSnapshotManager_SOURCES   = TestSnapshotManager.cxx
TestTileManipulation_SOURCES  = TestTileManipulation.cxx
TestTransactions_SOURCES      = TestTransactions.cxx

if HAVE_PYTHON
check_SCRIPTS  += TestModPlateRegression.py
endif

check_PROGRAMS = \
  TestBlobIO \
  TestBlobManager \
  TestDatastore \
  TestHTTPUtils \
  TestIndexPage \
  TestLocalIndex \
  TestModPlate \
  TestPlateManager \
  TestRpc \
  TestRpcChannel \
  TestSnapshotManager \
  TestTileManipulation \
  TestTransactions

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

TESTS = $(check_PROGRAMS)
CLEANFILES =

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
include $(top_srcdir)/thirdparty/protobuf.mak
