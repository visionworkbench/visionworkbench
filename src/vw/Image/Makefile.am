# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_IMAGE

include_HEADERS = \
  Algorithms.h \
  BlockProcessor.h \
  BlockRasterize.h \
  Convolution.h \
  EdgeExtend.h \
  EdgeExtension.h \
  Filter.h \
  Filter.tcc \
  ImageIO.h \
  ImageMath.h \
  ImageResource.h \
  ImageResourceImpl.h \
  ImageResourceStream.h \
  ImageResourceView.h \
  ImageViewBase.h \
  ImageView.h \
  ImageViewRef.h \
  Interpolation.h \
  Manipulation.h \
  MaskViews.h \
  Palette.h \
  PerPixelAccessorViews.h \
  PerPixelViews.h \
  PixelAccessors.h \
  PixelIterator.h \
  PixelMask.h \
  PixelMath.h \
  PixelTypeInfo.h \
  PixelTypes.h \
  SparseImageCheck.h \
  Statistics.h \
  Transform.h \
  UtilityViews.h \
  ViewImageResource.h

libvwImage_la_SOURCES = \
  Filter.cc \
  ImageResource.cc \
  ImageResourceStream.cc \
  Interpolation.cc \
  PixelTypeInfo.cc

libvwImage_la_LIBADD = @MODULE_IMAGE_LIBS@

lib_LTLIBRARIES = libvwImage.la

if HAVE_PKG_OPENCV
include_HEADERS += ImageResourceOpenCV.h
libvwImage_la_SOURCES += ImageResourceOpenCV.cc
endif

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Image

include $(top_srcdir)/config/rules.mak
