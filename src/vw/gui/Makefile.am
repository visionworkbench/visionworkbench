# __BEGIN_LICENSE__
# Copyright (C) 2006-2011 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

# 'make clean' should get rid of any QT meta-object compiled files.
CLEANFILES = *.moc.cc

if ENABLE_EXCEPTIONS

noinst_HEADERS =                \
  GlPreviewWidget.h             \
  ImageTileGenerator.h          \
  MainWindow.h                  \
  PlatefileTileGenerator.h      \
  TestPatternTileGenerator.h    \
  TextureCache.h                \
  TileGenerator.h               \
  Tree.h                        \
  WebTileGenerator.h

vwv_SOURCES =                   \
  $(BUILT_SOURCES)              \
  GlPreviewWidget.cc            \
  ImageTileGenerator.cc         \
  MainWindow.cc                 \
  TestPatternTileGenerator.cc   \
  TextureCache.cc               \
  TileGenerator.cc              \
  vwv.cc

BUILT_SOURCES =                 \
  GlPreviewWidget.moc.cc        \
  MainWindow.moc.cc

vwv_LDADD    = @PKG_GUI_LIBS@

if MAKE_MODULE_PLATE
  vwv_LDADD += @PKG_PLATE_LIBS@
  vwv_SOURCES += PlatefileTileGenerator.cc \
                 WebTileGenerator.cc
  BUILT_SOURCES += WebTileGenerator.moc.cc
endif

bin_PROGRAMS = vwv

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ @PKG_GUI_CPPFLAGS@
if MAKE_MODULE_PLATE
  AM_CPPFLAGS += @PKG_PLATE_CPPFLAGS@
endif
AM_LDFLAGS = @VW_LDFLAGS@

includedir = $(prefix)/include/vw

include $(top_srcdir)/config/rules.mak
