# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

StandardShaders.cc:
	perl StandardShaders.pl

CLEANFILES = StandardShaders.cc

if MAKE_MODULE_GPU

include_HEADERS = Setup.h Utilities.h StandardShaders.h GPUProgram.h	\
  TexAlloc.h TexObj.h GPUImage.h Filter.h Manipulation.h Transform.h	\
  GenericShaders.h ImageMath.h Statistics.h Algorithms.h		\
  Interpolation.h EdgeExtension.h Expressions.h

libvwGPU_la_SOURCES = Setup.cc Utilities.cc StandardShaders.cc	\
  GPUProgram.cc TexAlloc.cc TexObj.cc GPUImage.cc Filter.cc	\
  Manipulation.cc Transform.cc GenericShaders.cc ImageMath.cc	\
  Statistics.cc Algorithms.cc Expressions.cc

libvwGPU_la_LIBADD = @MODULE_GPU_LIBS@

lib_LTLIBRARIES = libvwGPU.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@
EXTRA_DIST= \
	StandardShaders.pl \
	StandardShaders/VWGPU_Algorithms/clamp-1i2f.cg \
	StandardShaders/VWGPU_Algorithms/clamp-1i2f.glsl \
	StandardShaders/VWGPU_Algorithms/copy-1i0f.cg \
	StandardShaders/VWGPU_Algorithms/copy-1i0f.glsl \
	StandardShaders/VWGPU_Algorithms/fill-1i4f.cg \
	StandardShaders/VWGPU_Algorithms/fill-1i4f.glsl \
	StandardShaders/VWGPU_Algorithms/normalize-1i3f.cg \
	StandardShaders/VWGPU_Algorithms/normalize-1i3f.glsl \
	StandardShaders/VWGPU_Algorithms/threshold-1i3f-r.cg \
	StandardShaders/VWGPU_Algorithms/threshold-1i3f-rgba.cg \
	StandardShaders/VWGPU_Algorithms/threshold-1i3f-rgba.glsl \
	StandardShaders/VWGPU_Algorithms/threshold-1i3f-rgb.cg \
	StandardShaders/VWGPU_Algorithms/threshold-1i3f-rgb.glsl \
	StandardShaders/VWGPU_Algorithms/threshold-1i3f-r.glsl \
	StandardShaders/VWGPU_Filter/convolution.cg \
	StandardShaders/VWGPU_Filter/convolution-columns.cg \
	StandardShaders/VWGPU_Filter/convolution-columns.glsl \
	StandardShaders/VWGPU_Filter/convolution.glsl \
	StandardShaders/VWGPU_Filter/convolution-rows.cg \
	StandardShaders/VWGPU_Filter/convolution-rows.glsl \
	StandardShaders/VWGPU_ImageMath/abs-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/abs-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/acos-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/acos-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/acosh-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/acosh-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/add-IF-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/add-IF-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/add-II-2i0f.cg \
	StandardShaders/VWGPU_ImageMath/add-II-2i0f.glsl \
	StandardShaders/VWGPU_ImageMath/asin-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/asin-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/asinh-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/asinh-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/atan-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/atan-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/atan2-FI-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/atan2-FI-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/atan2-IF-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/atan2-IF-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/atan2-II-2i0f.cg \
	StandardShaders/VWGPU_ImageMath/atan2-II-2i0f.glsl \
	StandardShaders/VWGPU_ImageMath/atanh-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/atanh-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/cbrt-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/cbrt-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/ceil-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/ceil-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/copysign-FI-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/copysign-FI-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/copysign-IF-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/copysign-IF-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/copysign-II-2i0f.cg \
	StandardShaders/VWGPU_ImageMath/copysign-II-2i0f.glsl \
	StandardShaders/VWGPU_ImageMath/cos-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/cos-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/cosh-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/cosh-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/divide-FI-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/divide-FI-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/divide-IF-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/divide-IF-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/divide-II-2i0f.cg \
	StandardShaders/VWGPU_ImageMath/divide-II-2i0f.glsl \
	StandardShaders/VWGPU_ImageMath/exp-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/exp-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/exp2-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/exp2-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/expm1-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/expm1-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/fdim-2i0f.glsl \
	StandardShaders/VWGPU_ImageMath/fdim-FI-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/fdim-FI-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/fdim-IF-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/fdim-IF-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/fdim-II-2i0f.cg \
	StandardShaders/VWGPU_ImageMath/fdim-II-2i0f.glsl \
	StandardShaders/VWGPU_ImageMath/floor-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/floor-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/hypot-FI-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/hypot-FI-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/hypot-IF-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/hypot-IF-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/hypot-II-2i0f.cg \
	StandardShaders/VWGPU_ImageMath/hypot-II-2i0f.glsl \
	StandardShaders/VWGPU_ImageMath/log10-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/log10-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/log-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/log-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/log1p-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/log1p-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/log2-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/log2-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/multiply-IF-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/multiply-IF-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/multiply-II-2i0f.cg \
	StandardShaders/VWGPU_ImageMath/multiply-II-2i0f.glsl \
	StandardShaders/VWGPU_ImageMath/pow-FI-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/pow-FI-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/pow-IF-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/pow-IF-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/pow-II-2i0f.cg \
	StandardShaders/VWGPU_ImageMath/pow-II-2i0f.glsl \
	StandardShaders/VWGPU_ImageMath/round-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/round-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/sin-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/sin-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/sinh-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/sinh-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/sqrt-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/sqrt-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/subtract-FI-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/subtract-FI-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/subtract-IF-1i1f.cg \
	StandardShaders/VWGPU_ImageMath/subtract-IF-1i1f.glsl \
	StandardShaders/VWGPU_ImageMath/subtract-II-2i0f.cg \
	StandardShaders/VWGPU_ImageMath/subtract-II-2i0f.glsl \
	StandardShaders/VWGPU_ImageMath/tan-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/tan-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/tanh-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/tanh-1i0f.glsl \
	StandardShaders/VWGPU_ImageMath/trunc-1i0f.cg \
	StandardShaders/VWGPU_ImageMath/trunc-1i0f.glsl \
	StandardShaders/VWGPU_Interp/identity.cg \
	StandardShaders/VWGPU_Interp/identity.glsl \
	StandardShaders/VWGPU_Interp/interpolation-bicubic.cg \
	StandardShaders/VWGPU_Interp/interpolation-bicubic.glsl \
	StandardShaders/VWGPU_Interp/interpolation-bilinear.cg \
	StandardShaders/VWGPU_Interp/interpolation-bilinear.glsl \
	StandardShaders/VWGPU_Interp/interpolation-nearest-pixel.cg \
	StandardShaders/VWGPU_Interp/interpolation-nearest-pixel.glsl \
	StandardShaders/VWGPU_Manipulation/PackGrayIntoRGBA.glsl \
	StandardShaders/VWGPU_Manipulation/SelectChannel-A.glsl \
	StandardShaders/VWGPU_Manipulation/SelectChannel-B.glsl \
	StandardShaders/VWGPU_Manipulation/SelectChannel-G.glsl \
	StandardShaders/VWGPU_Manipulation/SelectChannel-R.glsl \
	StandardShaders/VWGPU_Manipulation/transpose-1i0f.cg \
	StandardShaders/VWGPU_Manipulation/transpose-1i0f.glsl \
	StandardShaders/VWGPU_Manipulation/UnpackGrayFromRGBA-R.glsl \
	StandardShaders/VWGPU_Statistics/max-channels-r.cg \
	StandardShaders/VWGPU_Statistics/max-channels-rgba.cg \
	StandardShaders/VWGPU_Statistics/max-channels-rgba.glsl \
	StandardShaders/VWGPU_Statistics/max-channels-rgb.cg \
	StandardShaders/VWGPU_Statistics/max-channels-rgb.glsl \
	StandardShaders/VWGPU_Statistics/max-channels-r.glsl \
	StandardShaders/VWGPU_Statistics/max-quad.cg \
	StandardShaders/VWGPU_Statistics/max-quad.glsl \
	StandardShaders/VWGPU_Statistics/min-channels-r.cg \
	StandardShaders/VWGPU_Statistics/min-channels-rgba.cg \
	StandardShaders/VWGPU_Statistics/min-channels-rgba.glsl \
	StandardShaders/VWGPU_Statistics/min-channels-rgb.cg \
	StandardShaders/VWGPU_Statistics/min-channels-rgb.glsl \
	StandardShaders/VWGPU_Statistics/min-channels-r.glsl \
	StandardShaders/VWGPU_Statistics/min-quad.cg \
	StandardShaders/VWGPU_Statistics/min-quad.glsl \
	StandardShaders/VWGPU_Statistics/sum-channels-r.cg \
	StandardShaders/VWGPU_Statistics/sum-channels-rgba.cg \
	StandardShaders/VWGPU_Statistics/sum-channels-rgba.glsl \
	StandardShaders/VWGPU_Statistics/sum-channels-rgb.cg \
	StandardShaders/VWGPU_Statistics/sum-channels-rgb.glsl \
	StandardShaders/VWGPU_Statistics/sum-channels-r.glsl \
	StandardShaders/VWGPU_Statistics/sum-quad.cg \
	StandardShaders/VWGPU_Statistics/sum-quad.glsl

SUBDIRS = . tests

includedir = $(prefix)/include/vw/GPU

include $(top_srcdir)/config/rules.mak
