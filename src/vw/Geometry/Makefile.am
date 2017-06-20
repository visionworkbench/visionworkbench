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

if MAKE_MODULE_GEOMETRY


include_HEADERS = Shape.h SpatialTree.h PointListIO.h Sphere.h Box.h ATrans.h Frame.h TreeNode.h FrameTreeNode.h FrameStore.h FrameHandle.h geomUtils.h edgeUtils.h dPoly.h cutPoly.h baseUtils.h

libvwGeometry_la_SOURCES = SpatialTree.cc FrameTreeNode.cc FrameStore.cc geomUtils.cc edgeUtils.cc cutPoly.cc dPoly.cc
libvwGeometry_la_LIBADD = @MODULE_GEOMETRY_LIBS@

lib_LTLIBRARIES = libvwGeometry.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Geometry

include $(top_srcdir)/config/rules.mak
