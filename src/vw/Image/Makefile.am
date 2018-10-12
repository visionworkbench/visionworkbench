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

if MAKE_MODULE_IMAGE

include_HEADERS = \
  AlgorithmFunctions.h \
  AlgorithmFunctions.tcc \
  Algorithms.h \
  Algorithms.tcc \
  Algorithms2.h \
  AntiAliasing.h \
  BlobIndex.h \
  BlockProcessor.h \
  BlockRasterize.h \
  CensusTransform.h \
  Convolution.h \
  EdgeExtension.h \
  EdgeExtension.tcc \
  ErodeView.h \
  Filter.h \
  Filter.tcc \
  Fourier.h \
  ImageIO.h \
  ImageMath.h \
  ImageResource.h \
  ImageResourceImpl.h \
  ImageResourceStream.h \
  ImageResourceView.h \
  ImageViewBase.h \
  ImageView.h \
  ImageViewRef.h \
  InpaintView.h \
  Interpolation.h \
  Manipulation.h \
  Manipulation.tcc \
  MaskViews.h \
  MaskViews.tcc \
  PerPixelAccessorViews.h \
  PerPixelViews.h \
  PerPixelViews.tcc \
  PixelAccessors.h \
  PixelIterator.h \
  PixelMask.h \
  PixelMath.h \
  PixelTypeInfo.h \
  PixelTypes.h \
  SparseImageCheck.h \
  SparseView.h \
  Statistics.h \
  Statistics.tcc \
  Transform.h \
  UtilityViews.h \
  ViewImageResource.h \
  WindowAlgorithms.h

libvwImage_la_SOURCES = \
  BlobIndex.cc \
  Filter.cc \
  ImageResource.cc \
  ImageResourceStream.cc \
  Interpolation.cc \
  Transform.cc \
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
