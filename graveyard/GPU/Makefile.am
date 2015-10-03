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

Shaders.cc:
	perl "${top_srcdir}/${subdir}/Shaders.pl"

CLEANFILES = Shaders.cc

if MAKE_MODULE_GPU

include_HEADERS = Setup.h Utilities.h Shaders.h GPUProgram.h	\
  TexAlloc.h TexObj.h GPUImage.h Filter.h Manipulation.h Transform.h	\
  GenericShaders.h ImageMath.h Statistics.h Algorithms.h		\
  Interpolation.h EdgeExtension.h Expressions.h

libvwGPU_la_SOURCES = Setup.cc Utilities.cc Shaders.cc	\
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
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)
EXTRA_DIST= \
	Shaders.pl \
	Shaders/Algorithms/clamp-1i2f.cg \
	Shaders/Algorithms/clamp-1i2f.glsl \
	Shaders/Algorithms/copy-1i0f.cg \
	Shaders/Algorithms/copy-1i0f.glsl \
	Shaders/Algorithms/fill-1i4f.cg \
	Shaders/Algorithms/fill-1i4f.glsl \
	Shaders/Algorithms/normalize-1i3f.cg \
	Shaders/Algorithms/normalize-1i3f.glsl \
	Shaders/Algorithms/threshold-1i3f-r.cg \
	Shaders/Algorithms/threshold-1i3f-rgba.cg \
	Shaders/Algorithms/threshold-1i3f-rgba.glsl \
	Shaders/Algorithms/threshold-1i3f-rgb.cg \
	Shaders/Algorithms/threshold-1i3f-rgb.glsl \
	Shaders/Algorithms/threshold-1i3f-r.glsl \
	Shaders/Filter/convolution.cg \
	Shaders/Filter/convolution-columns.cg \
	Shaders/Filter/convolution-columns.glsl \
	Shaders/Filter/convolution.glsl \
	Shaders/Filter/convolution-rows.cg \
	Shaders/Filter/convolution-rows.glsl \
	Shaders/ImageMath/abs-1i0f.cg \
	Shaders/ImageMath/abs-1i0f.glsl \
	Shaders/ImageMath/acos-1i0f.cg \
	Shaders/ImageMath/acos-1i0f.glsl \
	Shaders/ImageMath/acosh-1i0f.cg \
	Shaders/ImageMath/acosh-1i0f.glsl \
	Shaders/ImageMath/add-IF-1i1f.cg \
	Shaders/ImageMath/add-IF-1i1f.glsl \
	Shaders/ImageMath/add-II-2i0f.cg \
	Shaders/ImageMath/add-II-2i0f.glsl \
	Shaders/ImageMath/asin-1i0f.cg \
	Shaders/ImageMath/asin-1i0f.glsl \
	Shaders/ImageMath/asinh-1i0f.cg \
	Shaders/ImageMath/asinh-1i0f.glsl \
	Shaders/ImageMath/atan-1i0f.cg \
	Shaders/ImageMath/atan-1i0f.glsl \
	Shaders/ImageMath/atan2-FI-1i1f.cg \
	Shaders/ImageMath/atan2-FI-1i1f.glsl \
	Shaders/ImageMath/atan2-IF-1i1f.cg \
	Shaders/ImageMath/atan2-IF-1i1f.glsl \
	Shaders/ImageMath/atan2-II-2i0f.cg \
	Shaders/ImageMath/atan2-II-2i0f.glsl \
	Shaders/ImageMath/atanh-1i0f.cg \
	Shaders/ImageMath/atanh-1i0f.glsl \
	Shaders/ImageMath/cbrt-1i0f.cg \
	Shaders/ImageMath/cbrt-1i0f.glsl \
	Shaders/ImageMath/ceil-1i0f.cg \
	Shaders/ImageMath/ceil-1i0f.glsl \
	Shaders/ImageMath/copysign-FI-1i1f.cg \
	Shaders/ImageMath/copysign-FI-1i1f.glsl \
	Shaders/ImageMath/copysign-IF-1i1f.cg \
	Shaders/ImageMath/copysign-IF-1i1f.glsl \
	Shaders/ImageMath/copysign-II-2i0f.cg \
	Shaders/ImageMath/copysign-II-2i0f.glsl \
	Shaders/ImageMath/cos-1i0f.cg \
	Shaders/ImageMath/cos-1i0f.glsl \
	Shaders/ImageMath/cosh-1i0f.cg \
	Shaders/ImageMath/cosh-1i0f.glsl \
	Shaders/ImageMath/divide-FI-1i1f.cg \
	Shaders/ImageMath/divide-FI-1i1f.glsl \
	Shaders/ImageMath/divide-IF-1i1f.cg \
	Shaders/ImageMath/divide-IF-1i1f.glsl \
	Shaders/ImageMath/divide-II-2i0f.cg \
	Shaders/ImageMath/divide-II-2i0f.glsl \
	Shaders/ImageMath/exp-1i0f.cg \
	Shaders/ImageMath/exp-1i0f.glsl \
	Shaders/ImageMath/exp2-1i0f.cg \
	Shaders/ImageMath/exp2-1i0f.glsl \
	Shaders/ImageMath/expm1-1i0f.cg \
	Shaders/ImageMath/expm1-1i0f.glsl \
	Shaders/ImageMath/fdim-2i0f.glsl \
	Shaders/ImageMath/fdim-FI-1i1f.cg \
	Shaders/ImageMath/fdim-FI-1i1f.glsl \
	Shaders/ImageMath/fdim-IF-1i1f.cg \
	Shaders/ImageMath/fdim-IF-1i1f.glsl \
	Shaders/ImageMath/fdim-II-2i0f.cg \
	Shaders/ImageMath/fdim-II-2i0f.glsl \
	Shaders/ImageMath/floor-1i0f.cg \
	Shaders/ImageMath/floor-1i0f.glsl \
	Shaders/ImageMath/hypot-FI-1i1f.cg \
	Shaders/ImageMath/hypot-FI-1i1f.glsl \
	Shaders/ImageMath/hypot-IF-1i1f.cg \
	Shaders/ImageMath/hypot-IF-1i1f.glsl \
	Shaders/ImageMath/hypot-II-2i0f.cg \
	Shaders/ImageMath/hypot-II-2i0f.glsl \
	Shaders/ImageMath/log10-1i0f.cg \
	Shaders/ImageMath/log10-1i0f.glsl \
	Shaders/ImageMath/log-1i0f.cg \
	Shaders/ImageMath/log-1i0f.glsl \
	Shaders/ImageMath/log1p-1i0f.cg \
	Shaders/ImageMath/log1p-1i0f.glsl \
	Shaders/ImageMath/log2-1i0f.cg \
	Shaders/ImageMath/log2-1i0f.glsl \
	Shaders/ImageMath/multiply-IF-1i1f.cg \
	Shaders/ImageMath/multiply-IF-1i1f.glsl \
	Shaders/ImageMath/multiply-II-2i0f.cg \
	Shaders/ImageMath/multiply-II-2i0f.glsl \
	Shaders/ImageMath/pow-FI-1i1f.cg \
	Shaders/ImageMath/pow-FI-1i1f.glsl \
	Shaders/ImageMath/pow-IF-1i1f.cg \
	Shaders/ImageMath/pow-IF-1i1f.glsl \
	Shaders/ImageMath/pow-II-2i0f.cg \
	Shaders/ImageMath/pow-II-2i0f.glsl \
	Shaders/ImageMath/round-1i0f.cg \
	Shaders/ImageMath/round-1i0f.glsl \
	Shaders/ImageMath/sin-1i0f.cg \
	Shaders/ImageMath/sin-1i0f.glsl \
	Shaders/ImageMath/sinh-1i0f.cg \
	Shaders/ImageMath/sinh-1i0f.glsl \
	Shaders/ImageMath/sqrt-1i0f.cg \
	Shaders/ImageMath/sqrt-1i0f.glsl \
	Shaders/ImageMath/subtract-FI-1i1f.cg \
	Shaders/ImageMath/subtract-FI-1i1f.glsl \
	Shaders/ImageMath/subtract-IF-1i1f.cg \
	Shaders/ImageMath/subtract-IF-1i1f.glsl \
	Shaders/ImageMath/subtract-II-2i0f.cg \
	Shaders/ImageMath/subtract-II-2i0f.glsl \
	Shaders/ImageMath/tan-1i0f.cg \
	Shaders/ImageMath/tan-1i0f.glsl \
	Shaders/ImageMath/tanh-1i0f.cg \
	Shaders/ImageMath/tanh-1i0f.glsl \
	Shaders/ImageMath/trunc-1i0f.cg \
	Shaders/ImageMath/trunc-1i0f.glsl \
	Shaders/Interp/identity.cg \
	Shaders/Interp/identity.glsl \
	Shaders/Interp/interpolation-bicubic.cg \
	Shaders/Interp/interpolation-bicubic.glsl \
	Shaders/Interp/interpolation-bilinear.cg \
	Shaders/Interp/interpolation-bilinear.glsl \
	Shaders/Interp/interpolation-nearest-pixel.cg \
	Shaders/Interp/interpolation-nearest-pixel.glsl \
	Shaders/Manipulation/PackGrayIntoRGBA.glsl \
	Shaders/Manipulation/SelectChannel-A.glsl \
	Shaders/Manipulation/SelectChannel-B.glsl \
	Shaders/Manipulation/SelectChannel-G.glsl \
	Shaders/Manipulation/SelectChannel-R.glsl \
	Shaders/Manipulation/transpose-1i0f.cg \
	Shaders/Manipulation/transpose-1i0f.glsl \
	Shaders/Manipulation/UnpackGrayFromRGBA-R.glsl \
	Shaders/Statistics/max-channels-r.cg \
	Shaders/Statistics/max-channels-rgba.cg \
	Shaders/Statistics/max-channels-rgba.glsl \
	Shaders/Statistics/max-channels-rgb.cg \
	Shaders/Statistics/max-channels-rgb.glsl \
	Shaders/Statistics/max-channels-r.glsl \
	Shaders/Statistics/max-quad.cg \
	Shaders/Statistics/max-quad.glsl \
	Shaders/Statistics/min-channels-r.cg \
	Shaders/Statistics/min-channels-rgba.cg \
	Shaders/Statistics/min-channels-rgba.glsl \
	Shaders/Statistics/min-channels-rgb.cg \
	Shaders/Statistics/min-channels-rgb.glsl \
	Shaders/Statistics/min-channels-r.glsl \
	Shaders/Statistics/min-quad.cg \
	Shaders/Statistics/min-quad.glsl \
	Shaders/Statistics/sum-channels-r.cg \
	Shaders/Statistics/sum-channels-rgba.cg \
	Shaders/Statistics/sum-channels-rgba.glsl \
	Shaders/Statistics/sum-channels-rgb.cg \
	Shaders/Statistics/sum-channels-rgb.glsl \
	Shaders/Statistics/sum-channels-r.glsl \
	Shaders/Statistics/sum-quad.cg \
	Shaders/Statistics/sum-quad.glsl

SUBDIRS = . tests

includedir = $(prefix)/include/vw/GPU

include $(top_srcdir)/config/rules.mak
