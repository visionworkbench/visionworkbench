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

MODULES= cartography composite core fileio filter image imagealgo imagemanip imagemath pixel pixelcast qtree transform vwmath

BUILT_SOURCES = $(MODULES:%=%_wrap.cpp)
CLEANFILES = $(BUILT_SOURCES) $(MODULES:%=%_wrap.h) $(MODULES:%=%.py)

_cartography_la_SOURCES = cartography_wrap.cpp _cartography.i
_cartography_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_cartography_la_LDFLAGS  = -module $(VW_LDFLAGS)
_cartography_la_LIBADD   = $(PKG_CARTOGRAPHY_LIBS)

_composite_la_SOURCES = composite_wrap.cpp _composite.i
_composite_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_composite_la_LDFLAGS = -module $(VW_LDFLAGS)
_composite_la_LIBADD = $(PKG_MOSAIC_LIBS)

_core_la_SOURCES = core_wrap.cpp _core.i
_core_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS)
_core_la_LDFLAGS = -module $(VW_LDFLAGS)
_core_la_LIBADD = $(PKG_VW_LIBS)

_fileio_la_SOURCES = fileio_wrap.cpp _fileio.i
_fileio_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_fileio_la_LDFLAGS = -module $(VW_LDFLAGS)
_fileio_la_LIBADD = $(PKG_VW_LIBS)

_filter_la_SOURCES = filter_wrap.cpp _filter.i
_filter_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_filter_la_LDFLAGS = -module $(VW_LDFLAGS)
_filter_la_LIBADD = $(PKG_VW_LIBS)

_image_la_SOURCES = image_wrap.cpp _image.i
_image_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_image_la_LDFLAGS = -module $(VW_LDFLAGS)
_image_la_LIBADD = $(PKG_VW_LIBS)

_imagealgo_la_SOURCES = imagealgo_wrap.cpp _imagealgo.i
_imagealgo_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_imagealgo_la_LDFLAGS = -module $(VW_LDFLAGS)
_imagealgo_la_LIBADD = $(PKG_VW_LIBS)

_imagemanip_la_SOURCES = imagemanip_wrap.cpp _imagemanip.i
_imagemanip_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_imagemanip_la_LDFLAGS = -module $(VW_LDFLAGS)
_imagemanip_la_LIBADD = $(PKG_VW_LIBS)

_imagemath_la_SOURCES = imagemath_wrap.cpp _imagemath.i
_imagemath_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_imagemath_la_LDFLAGS = -module $(VW_LDFLAGS)
_imagemath_la_LIBADD = $(PKG_VW_LIBS)

_pixel_la_SOURCES = pixel_wrap.cpp _pixel.i
_pixel_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_pixel_la_LDFLAGS = -module $(VW_LDFLAGS)
_pixel_la_LIBADD = $(PKG_VW_LIBS)

_pixelcast_la_SOURCES = pixelcast_wrap.cpp _pixelcast.i
_pixelcast_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_pixelcast_la_LDFLAGS = -module $(VW_LDFLAGS)
_pixelcast_la_LIBADD = $(PKG_VW_LIBS)

_qtree_la_SOURCES = qtree_wrap.cpp _qtree.i
_qtree_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_qtree_la_LDFLAGS = -module $(VW_LDFLAGS)
_qtree_la_LIBADD = $(PKG_MOSAIC_LIBS)

_transform_la_SOURCES = transform_wrap.cpp _transform.i
_transform_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_transform_la_LDFLAGS = -module $(VW_LDFLAGS)
_transform_la_LIBADD = $(PKG_VW_LIBS)

_vwmath_la_SOURCES = vwmath_wrap.cpp _vwmath.i
_vwmath_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_vwmath_la_LDFLAGS = -module $(VW_LDFLAGS)
_vwmath_la_LIBADD = $(PKG_VW_LIBS)

python_vwdir = ${pythondir}/vw
pyexec_vwdir = ${pyexecdir}/vw

python_vw_PYTHON = __init__.py cartography.py composite.py core.py fileio.py filter.py image.py imagealgo.py imagemanip.py imagemath.py pixel.py pixelcast.py qtree.py transform.py vwmath.py
pyexec_vw_LTLIBRARIES = _cartography.la _composite.la _core.la _fileio.la _filter.la _image.la _imagealgo.la _imagemanip.la _imagemath.la _pixel.la _pixelcast.la _qtree.la _transform.la _vwmath.la

%_wrap.cpp : _%.i
	$(SWIG) $(SWIG_PYTHON_OPT) $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) -o $@ $<
