# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
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
vwv_SOURCES = vwv.cc MainWindow.cc GlPreviewWidget.cc MainWindow.h	\
	      GlPreviewWidget.h TileGenerator.h TextureCache.cc	\
	      TileGenerator.cc

# These are automatically-generated sources?
nodist_vwv_SOURCES = MainWindow.moc.cc GlPreviewWidget.moc.cc	\
		TileGenerator.moc.cc

vwv_LDADD    = @PKG_GUI_LIBS@ @PKG_PLATE_LIBS@

bin_PROGRAMS = vwv

noinst_PROGRAMS =

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ @PKG_GUI_CPPFLAGS@ @PKG_PLATE_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@

includedir = $(prefix)/include/vw

include $(top_srcdir)/config/rules.mak
