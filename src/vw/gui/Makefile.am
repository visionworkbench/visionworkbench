# __BEGIN_LICENSE__
# Copyright (C) 2006-2009 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
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
	      vwv_TextureCache.cc TileGenerator.cc

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
