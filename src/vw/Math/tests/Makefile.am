# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_MATH

if HAVE_PKG_LAPACK
TestLinearAlgebra = TestLinearAlgebra TestGeometry TestLevenbergMarquardt TestPoseEstimation
endif

TESTS = TestVector TestMatrix TestQuaternion TestBBox TestFunctions	\
        TestFunctors TestNelderMead TestKDTree $(TestLinearAlgebra)	\
        TestEuler TestParticleSwarmOptimization

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_MATH_LIBS@ @PKG_CORE_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = *.cxx

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
