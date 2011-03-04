# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_MATH

if HAVE_PKG_LAPACK
lapack_headers = LinearAlgebra.h LevenbergMarquardt.h Geometry.h PoseEstimation.h LapackExports.h
lapack_sources = LapackExports.cc
endif

if HAVE_PKG_FLANN
flann_headers = FLANNTree.h
endif

include_HEADERS = Vector.h Matrix.h BBox.h Functions.h Functors.h	\
                  Quaternion.h EulerAngles.h ConjugateGradient.h	\
                  NelderMead.h Statistics.h DisjointSet.h		\
                  MinimumSpanningTree.h KDTree.h ParticleSwarmOptimization.h \
                  RANSAC.h MatrixSparseSkyline.h $(lapack_headers) $(flann_headers)

libvwMath_la_SOURCES = MinimumSpanningTree.cc $(lapack_sources)
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
