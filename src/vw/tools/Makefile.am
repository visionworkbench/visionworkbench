# __BEGIN_LICENSE__
# 
# Copyright (C) 2006 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration
# (NASA).  All Rights Reserved.
# 
# Copyright 2006 Carnegie Mellon University. All rights reserved.
# 
# This software is distributed under the NASA Open Source Agreement
# (NOSA), version 1.3.  The NOSA has been approved by the Open Source
# Initiative.  See the file COPYING at the top of the distribution
# directory tree for the complete NOSA document.
# 
# THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
# KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
# LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
# SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
# THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
# DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
# 
# __END_LICENSE__

########################################################################
# sources
########################################################################

if ENABLE_EXCEPTIONS

if MAKE_MODULE_CARTOGRAPHY
if HAVE_PKG_GDAL
if MAKE_MODULE_MOSAIC
image2kml_prog = image2kml
image2kml_SOURCES = image2kml.cc
image2kml_LDFLAGS = @PKG_VW_LDFLAGS@ @PKG_CARTOGRAPHY_LDFLAGS@ @PKG_MOSAIC_LDFLAGS@ @PKG_BOOST_THREAD_LDFLAGS@

image2tms_prog = image2tms
image2tms_SOURCES = image2tms.cc
image2tms_LDFLAGS = @PKG_VW_LDFLAGS@ @PKG_CARTOGRAPHY_LDFLAGS@ @PKG_MOSAIC_LDFLAGS@ @PKG_BOOST_THREAD_LDFLAGS@

geoblend_prog = geoblend
geoblend_SOURCES = geoblend.cc
geoblend_LDFLAGS = @PKG_VW_LDFLAGS@ @PKG_CARTOGRAPHY_LDFLAGS@ @PKG_MOSAIC_LDFLAGS@ @PKG_BOOST_THREAD_LDFLAGS@
endif

georef_prog = georef
georef_SOURCES = georef.cc
georef_LDFLAGS = @PKG_VW_LDFLAGS@ @PKG_CARTOGRAPHY_LDFLAGS@ @PKG_BOOST_THREAD_LDFLAGS@
endif
endif

if HAVE_PKG_JPEG2K
jp2dump_prog = jp2dump
jp2dump_SOURCES = jp2dump.cc
jp2dump_LDFLAGS = @PKG_VW_LDFLAGS@
endif

bin_PROGRAMS = $(image2kml_prog) $(image2tms_prog) $(geoblend_prog) $(georef_prog) $(jp2dump_prog)

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @PKG_VW_CPPFLAGS@
AM_LDFLAGS = @PKG_BOOST_PROGRAM_OPTIONS_LDFLAGS@

includedir = $(prefix)/include/vw

include $(top_srcdir)/config/rules.mak
