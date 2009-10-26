# __BEGIN_LICENSE__
# Copyright (C) 2006-2009 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_MATH

TestVector_SOURCES                    = TestVector.cxx
TestMatrix_SOURCES                    = TestMatrix.cxx
TestQuaternion_SOURCES                = TestQuaternion.cxx
TestBBox_SOURCES                      = TestBBox.cxx
TestFunctions_SOURCES                 = TestFunctions.cxx
TestFunctors_SOURCES                  = TestFunctors.cxx
TestNelderMead_SOURCES                = TestNelderMead.cxx
TestKDTree_SOURCES                    = TestKDTree.cxx
TestEuler_SOURCES                     = TestEuler.cxx
TestParticleSwarmOptimization_SOURCES = TestParticleSwarmOptimization.cxx

if HAVE_PKG_LAPACK

TestLinearAlgebra_SOURCES      = TestLinearAlgebra.cxx
TestGeometry_SOURCES           = TestGeometry.cxx
TestLevenbergMarquardt_SOURCES = TestLevenbergMarquardt.cxx
TestPoseEstimation_SOURCES     = TestPoseEstimation.cxx

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
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_MATH_LIBS@

check_PROGRAMS = $(TESTS)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
