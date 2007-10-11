# __BEGIN_LICENSE__
# 
# Copyright (C) 2006 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration
# (NASA).  All Rights Reserved.
# 
# Copyright 2006 Carnegie Mellon University. All rights reserved.
# 
# This software is distributed under the NASA Open Source Agreement
# (NOSA), version 1.3.  The NOSA has been approved by the Open Source
# Initiative.  See the file COPYING at the top of the distribution
# directory tree for the complete NOSA document.
# 
# THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
# KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
# LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
# SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
# THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
# DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
# 
# __END_LICENSE__

########################################################################
# sources
########################################################################

if MAKE_MODULE_GPU

TESTS = TestGPU_Profile TestGPU_ImageMath TestGPU_EdgeExtend TestGPU_Filter \
TestGPU_Interpolation TestGPU_Manipulation TestGPU_Statistics TestGPU_Transform \
TestGPU_Algorithms

nodist_TestGPU_Profile_SOURCES = TestGPU_Profile.cxx 
nodist_TestGPU_ImageMath_SOURCES = TestGPU_ImageMath.cxx
nodist_TestGPU_EdgeExtend_SOURCES = TestGPU_EdgeExtend.cxx
nodist_TestGPU_Filter_SOURCES = TestGPU_Filter.cxx
nodist_TestGPU_Interpolation_SOURCES = TestGPU_Interpolation.cxx
nodist_TestGPU_Manipulation_SOURCES = TestGPU_Manipulation.cxx
nodist_TestGPU_Statistics_SOURCES = TestGPU_Statistics.cxx
nodist_TestGPU_Transform_SOURCES = TestGPU_Transform.cxx
nodist_TestGPU_Algorithms_SOURCES = TestGPU_Algorithms.cxx

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ -I$(top_srcdir)/thirdparty/cxxtest
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_GPU_LIBS@ @PKG_VW_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(patsubst %,%.h,$(TESTS))
CLEANFILES = $(patsubst %,%.cxx,$(TESTS))

include $(top_srcdir)/config/rules.mak
