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

if MAKE_MODULE_GPU

#TestGPU_Profile_SOURCES       = TestGPU_Profile.oldtest
#TestGPU_ImageMath_SOURCES     = TestGPU_ImageMath.oldtest
#TestGPU_EdgeExtend_SOURCES    = TestGPU_EdgeExtend.oldtest
#TestGPU_Filter_SOURCES        = TestGPU_Filter.oldtest
#TestGPU_Interpolation_SOURCES = TestGPU_Interpolation.oldtest
#TestGPU_Manipulation_SOURCES  = TestGPU_Manipulation.oldtest
#TestGPU_Statistics_SOURCES    = TestGPU_Statistics.oldtest
#TestGPU_Transform_SOURCES     = TestGPU_Transform.oldtest
#TestGPU_Algorithms_SOURCES    = TestGPU_Algorithms.oldtest
#
#TESTS = TestGPU_Profile TestGPU_ImageMath TestGPU_EdgeExtend TestGPU_Filter \
#TestGPU_Interpolation TestGPU_Manipulation TestGPU_Statistics TestGPU_Transform \
#TestGPU_Algorithms

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_GPU_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = test_images/lighthouse-1000x1000.png

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
