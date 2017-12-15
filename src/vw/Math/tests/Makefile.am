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
TestStatistics_SOURCES                = TestStatistics.cxx
TestMatrixSparseSkyline_SOURCES       = TestMatrixSparseSkyline.cxx
TestConjugateGradient_SOURCES         = TestConjugateGradient.cxx
TestFLANNTree_SOURCES                 = TestFLANNTree.cxx
TestGaussianClustering_SOURCES        = TestGaussianClustering.cxx

if HAVE_PKG_LAPACK

TestLinearAlgebra_SOURCES      = TestLinearAlgebra.cxx
TestGeometry_SOURCES           = TestGeometry.cxx
TestLevenbergMarquardt_SOURCES = TestLevenbergMarquardt.cxx
TestPoseEstimation_SOURCES     = TestPoseEstimation.cxx

TestLinearAlgebra = TestLinearAlgebra TestGeometry TestLevenbergMarquardt TestPoseEstimation
endif

TESTS = TestVector TestMatrix TestQuaternion TestBBox TestFunctions     \
        TestFunctors TestNelderMead TestKDTree $(TestLinearAlgebra)     \
        TestEuler TestParticleSwarmOptimization TestStatistics          \
        TestMatrixSparseSkyline TestConjugateGradient TestFLANNTree     \
        TestGaussianClustering

#include $(top_srcdir)/config/instantiate.am

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_MATH_LIBS@

check_PROGRAMS = $(TESTS)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
