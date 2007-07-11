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

if MAKE_MODULE_MATH

nodist_TestVector_SOURCES = TestVector.cxx
nodist_TestMatrix_SOURCES = TestMatrix.cxx
nodist_TestQuaternion_SOURCES = TestQuaternion.cxx
nodist_TestBBox_SOURCES = TestBBox.cxx
nodist_TestFunctions_SOURCES = TestFunctions.cxx
nodist_TestFunctors_SOURCES = TestFunctors.cxx
nodist_TestNelderMead_SOURCES = TestNelderMead.cxx
#nodist_TestConjugateGradient_SOURCES = TestConjugateGradient.cxx
nodist_TestSpatialTree_SOURCES = TestSpatialTree.cxx

if HAVE_PKG_LAPACK
nodist_TestLinearAlgebra_SOURCES = TestLinearAlgebra.cxx
nodist_TestGeometry_SOURCES = TestGeometry.cxx
nodist_TestLevenbergMarquardt_SOURCES = TestLevenbergMarquardt.cxx
TestLinearAlgebra = TestLinearAlgebra TestGeometry TestLevenbergMarquardt
endif

#if HAVE_PKG_PPL
#nodist_TestBConvex_SOURCES = TestBConvex.cxx
#TestPPL = TestBConvex
#endif

TESTS = TestVector TestMatrix TestQuaternion TestBBox TestFunctions	 \
        TestFunctors TestNelderMead TestSpatialTree $(TestLinearAlgebra) \
        $(TestPPL) #TestConjugateGradient
endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @MODULE_MATH_CPPFLAGS@ -I$(top_srcdir)/thirdparty/cxxtest
AM_LDFLAGS  = @PKG_MATH_LDFLAGS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(patsubst %,%.h,$(TESTS))
CLEANFILES = $(patsubst %,%.cxx,$(TESTS))

include $(top_srcdir)/config/rules.mak
