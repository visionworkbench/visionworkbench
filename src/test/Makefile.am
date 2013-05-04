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

if MAKE_MODULE_CARTOGRAPHY
if HAVE_PKG_GDAL
geotif_generate_SOURCES = geotif-generate.cc
geotif_generate_CPPFLAGS = @VW_CPPFLAGS@ @PKG_CARTOGRAPHY_CPPFLAGS@
geotif_generate_LDADD = @PKG_CARTOGRAPHY_LIBS@
noinst_PROGRAMS = geotif-generate
endif
endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@

EXTRA_DIST = Helpers.h test_main.cc image2qtree_test.py
clean-local:
	-rm -rf image2qtree.out

includedir = $(prefix)/include/test
include $(top_srcdir)/config/rules.mak
