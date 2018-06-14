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

TestAlgorithms_SOURCES            = TestAlgorithms.cxx
TestAntiAliasing_SOURCES          = TestAntiAliasing.cxx
TestBlobIndex_SOURCES             = TestBlobIndex.cxx
TestBlockRasterize_SOURCES        = TestBlockRasterize.cxx
TestCensusTransform_SOURCES       = TestCensusTransform.cxx
TestConvolution_SOURCES           = TestConvolution.cxx
TestEdgeExtension_SOURCES         = TestEdgeExtension.cxx
TestErodeView_SOURCES             = TestErodeView.cxx
TestFilter_SOURCES                = TestFilter.cxx
TestImageMath_SOURCES             = TestImageMath.cxx
TestImageResource_SOURCES         = TestImageResource.cxx
TestImageViewRef_SOURCES          = TestImageViewRef.cxx
TestImageView_SOURCES             = TestImageView.cxx
TestImageViewMemory_SOURCES       = TestImageViewMemory.cxx
TestInpaintView_SOURCES           = TestInpaintView.cxx
TestInterpolation_SOURCES         = TestInterpolation.cxx
TestManipulation_SOURCES          = TestManipulation.cxx
TestMaskedImageMath_SOURCES       = TestMaskedImageMath.cxx
TestMaskedPixelMath2_SOURCES      = TestMaskedPixelMath2.cxx
TestMaskedPixelMath_SOURCES       = TestMaskedPixelMath.cxx
TestMaskViews_SOURCES             = TestMaskViews.cxx
TestPerPixelAccessorViews_SOURCES = TestPerPixelAccessorViews.cxx
TestPerPixelViews_SOURCES         = TestPerPixelViews.cxx
TestPixelMath_SOURCES             = TestPixelMath.cxx
TestPixelTypes_SOURCES            = TestPixelTypes.cxx
TestStatistics_SOURCES            = TestStatistics.cxx
TestTransform_SOURCES             = TestTransform.cxx
TestUtilityViews_SOURCES          = TestUtilityViews.cxx

TESTS = \
  TestAlgorithms \
  TestAntiAliasing \
  TestBlobIndex \
  TestBlockRasterize \
  TestCensusTransform \
  TestConvolution \
  TestEdgeExtension \
  TestErodeView \
  TestFilter \
  TestImageMath \
  TestImageResource \
  TestImageView \
  TestImageViewMemory \
  TestImageViewRef \
  TestInpaintView \
  TestInterpolation \
  TestManipulation \
  TestMaskedImageMath \
  TestMaskedPixelMath \
  TestMaskedPixelMath2 \
  TestMaskViews \
  TestPerPixelAccessorViews \
  TestPerPixelViews \
  TestPixelMath \
  TestPixelTypes \
  TestStatistics \
  TestTransform \
  TestUtilityViews

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_IMAGE_LIBS@ # @PKG_FILEIO_LIBS@

check_PROGRAMS = $(TESTS)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
