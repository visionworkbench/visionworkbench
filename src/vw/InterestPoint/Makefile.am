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

if MAKE_MODULE_INTERESTPOINT

include_HEADERS = Detector.h Detector.tcc Descriptor.h Matcher.h Extrema.h \
                  Localize.h InterestOperator.h WeightedHistogram.h    \
                  ImageOctave.h InterestData.h ImageOctaveHistory.h    \
                  InterestTraits.h MatrixIO.h LearnPCA.h               \
		  IntegralImage.h IntegralInterestOperator.h           \
		  IntegralDetector.h BoxFilter.h IntegralDescriptor.h

libvwInterestPoint_la_SOURCES = InterestData.cc Descriptor.cc   \
	          IntegralInterestOperator.cc Matcher.cc
libvwInterestPoint_la_LIBADD = @MODULE_INTERESTPOINT_LIBS@

lib_LTLIBRARIES = libvwInterestPoint.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/InterestPoint

include $(top_srcdir)/config/rules.mak
