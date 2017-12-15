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

# TODO: This should be required
if HAVE_PKG_LAPACK
lapack_headers = LinearAlgebra.h LevenbergMarquardt.h PoseEstimation.h LapackExports.h
lapack_sources = LapackExports.cc
endif

# TODO: This should be required
if HAVE_PKG_FLANN
flann_headers = FLANNTree.h
flann_sources = FLANNTree.cc
endif

include_HEADERS = Geometry.h Vector.h Matrix.h BBox.h BBox.tcc Functions.h Functors.h	\
		  Quaternion.h EulerAngles.h ConjugateGradient.h	\
		  NelderMead.h Statistics.h Statistics.tcc DisjointSet.h		\
		  MinimumSpanningTree.h KDTree.h ParticleSwarmOptimization.h \
		  BresenhamLine.h GaussianClustering.h \
		  RANSAC.h MatrixSparseSkyline.h $(lapack_headers) $(flann_headers)

libvwMath_la_SOURCES = Geometry.cc Quaternion.cc MinimumSpanningTree.cc $(lapack_sources) $(flann_sources)
libvwMath_la_LIBADD = @MODULE_MATH_LIBS@

lib_LTLIBRARIES = libvwMath.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Math

include $(top_srcdir)/config/rules.mak
