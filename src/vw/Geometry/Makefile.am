# __BEGIN_LICENSE__
# Copyright (C) 2006-2009 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_GEOMETRY


include_HEADERS = Shape.h SpatialTree.h PointListIO.h Sphere.h Box.h ATrans.h Frame.h TreeNode.h FrameTreeNode.h

libvwGeometry_la_SOURCES = SpatialTree.cc FrameTreeNode.cpp
libvwGeometry_la_LIBADD = @MODULE_GEOMETRY_LIBS@

lib_LTLIBRARIES = libvwGeometry.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Geometry

include $(top_srcdir)/config/rules.mak
