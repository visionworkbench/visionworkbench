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

if HAVE_PKG_LAPACK
lapack_headers = LinearAlgebra.h LevenbergMarquardt.h Geometry.h
endif

if HAVE_PKG_PPL
if HAVE_PKG_GMP
ppl_headers = BConvex.h
ppl_sources = BConvex.cc
endif
endif

include_HEADERS = Vector.h Matrix.h BShape.h BBox.h BBall.h	\
                  Functions.h Functors.h Quaternion.h EulerAngles.h	\
                  ConjugateGradient.h NelderMead.h SpatialTree.h	\
                  Statistics.h DisjointSet.h MinimumSpanningTree.h	\
                  $(lapack_headers) $(ppl_headers)

libvwMath_la_SOURCES = SpatialTree.cc MinimumSpanningTree.cc $(ppl_sources)
libvwMath_la_LDFLAGS= -version-info @LIBTOOL_VERSION@
libvwMath_la_LIBADD = @MODULE_MATH_LDFLAGS@

lib_LTLIBRARIES = libvwMath.la

endif

########################################################################
# general
########################################################################

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Math

AM_CPPFLAGS = @MODULE_MATH_CPPFLAGS@
AM_LDFLAGS  = 

include $(top_srcdir)/config/rules.mak
