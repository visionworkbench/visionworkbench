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

# 'make clean' should get rid of any QT meta-object compiled files.
CLEANFILES = *.moc.cc

# For simplicity, these example applications depend on exceptions
if ENABLE_EXCEPTIONS

# Graphical tools based on the core modules
vwv_SOURCES = vwv.cc vwv_MainWindow.cc vwv_GlPreviewWidget.cc	\
	      vwv_MainWindow.h vwv_GlPreviewWidget.h		\
	      vwv_TextureCache.cc

# These are automatically-generated sources?
nodist_vwv_SOURCES = vwv_MainWindow.moc.cc  vwv_GlPreviewWidget.moc.cc

vwv_LDADD    = @PKG_GUI_LIBS@

bin_PROGRAMS = vwv

noinst_PROGRAMS = 

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@

includedir = $(prefix)/include/vw

include $(top_srcdir)/config/rules.mak
