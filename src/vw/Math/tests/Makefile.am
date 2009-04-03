# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
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
nodist_TestKDTree_SOURCES = TestKDTree.cxx
nodist_TestEuler_SOURCES = TestEuler.cxx
nodist_TestParticleSwarmOptimization_SOURCES = TestParticleSwarmOptimization.cxx
#nodist_TestConjugateGradient_SOURCES = TestConjugateGradient.cxx

if HAVE_PKG_LAPACK
nodist_TestLinearAlgebra_SOURCES = TestLinearAlgebra.cxx
nodist_TestGeometry_SOURCES = TestGeometry.cxx
nodist_TestLevenbergMarquardt_SOURCES = TestLevenbergMarquardt.cxx
nodist_TestPoseEstimation_SOURCES = TestPoseEstimation.cxx
TestLinearAlgebra = TestLinearAlgebra TestGeometry	\
	TestLevenbergMarquardt TestPoseEstimation
endif

TESTS = TestVector TestMatrix TestQuaternion TestBBox TestFunctions	\
        TestFunctors TestNelderMead TestKDTree $(TestLinearAlgebra)	\
        TestEuler TestParticleSwarmOptimization #TestConjugateGradient

include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ $(TEST_CPPFLAGS)
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_MATH_LIBS@ @PKG_CORE_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(TESTS:%=%.h)
CLEANFILES = $(TESTS:%=%.cxx)

include $(top_srcdir)/config/rules.mak
