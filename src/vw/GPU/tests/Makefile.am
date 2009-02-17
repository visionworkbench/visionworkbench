# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_GPU

nodist_TestGPU_Profile_SOURCES = TestGPU_Profile.cxx 
nodist_TestGPU_ImageMath_SOURCES = TestGPU_ImageMath.cxx
nodist_TestGPU_EdgeExtend_SOURCES = TestGPU_EdgeExtend.cxx
nodist_TestGPU_Filter_SOURCES = TestGPU_Filter.cxx
nodist_TestGPU_Interpolation_SOURCES = TestGPU_Interpolation.cxx
nodist_TestGPU_Manipulation_SOURCES = TestGPU_Manipulation.cxx
nodist_TestGPU_Statistics_SOURCES = TestGPU_Statistics.cxx
nodist_TestGPU_Transform_SOURCES = TestGPU_Transform.cxx
nodist_TestGPU_Algorithms_SOURCES = TestGPU_Algorithms.cxx

TESTS = TestGPU_Profile TestGPU_ImageMath TestGPU_EdgeExtend TestGPU_Filter \
TestGPU_Interpolation TestGPU_Manipulation TestGPU_Statistics TestGPU_Transform \
TestGPU_Algorithms

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ -I$(top_srcdir)/thirdparty/cxxtest
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_GPU_LIBS@ @PKG_VW_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(TESTS:%=%.h)
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
