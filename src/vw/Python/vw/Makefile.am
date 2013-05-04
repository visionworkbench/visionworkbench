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


# cartography composite core fileio filter image imagealgo imagemanip imagemath pixel pixelcast qtree transform vwmath

EXTRA_DIST = \
    _cartography.i _composite.i _core.i _fileio.i _filter.i _image.i           \
    _imagealgo.i _imagemanip.i _imagemath.i _pixel.i _pixelcast.i              \
    _qtree.i _qtree.h _transform.i _vwmath.i numpy.i vwutil.i __init__.py

CLEANFILES = \
    cartography_wrap.h composite_wrap.h core_wrap.h fileio_wrap.h              \
        filter_wrap.h image_wrap.h imagealgo_wrap.h imagemanip_wrap.h          \
        imagemath_wrap.h pixel_wrap.h pixelcast_wrap.h qtree_wrap.h            \
        transform_wrap.h vwmath_wrap.h                                         \
    cartography.py composite.py core.py fileio.py filter.py image.py           \
        imagealgo.py imagemanip.py imagemath.py pixel.py pixelcast.py qtree.py \
        transform.py vwmath.py                                                 \
    _cartography.cpp _composite.cpp _core.cpp _fileio.cpp _filter.cpp          \
        _image.cpp _imagealgo.cpp _imagemanip.cpp _imagemath.cpp _pixel.cpp    \
        _pixelcast.cpp _qtree.cpp _transform.cpp _vwmath.cpp                   \
     _cartography.h _transform.h

if MAKE_MODULE_PYTHON

nodist__cartography_la_SOURCES = _cartography.cpp
_cartography_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_cartography_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_cartography_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_CARTOGRAPHY_LIBS)

nodist__composite_la_SOURCES = _composite.cpp
_composite_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_composite_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_composite_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_MOSAIC_LIBS)

nodist__core_la_SOURCES = _core.cpp
_core_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS)
_core_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_core_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_VW_LIBS)

nodist__fileio_la_SOURCES = _fileio.cpp
_fileio_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_fileio_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_fileio_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_VW_LIBS)

nodist__filter_la_SOURCES = _filter.cpp
_filter_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_filter_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_filter_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_VW_LIBS)

nodist__image_la_SOURCES = _image.cpp
_image_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_image_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_image_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_VW_LIBS)

nodist__imagealgo_la_SOURCES = _imagealgo.cpp
_imagealgo_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_imagealgo_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_imagealgo_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_VW_LIBS)

nodist__imagemanip_la_SOURCES = _imagemanip.cpp
_imagemanip_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_imagemanip_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_imagemanip_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_VW_LIBS)

nodist__imagemath_la_SOURCES = _imagemath.cpp
_imagemath_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_imagemath_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_imagemath_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_VW_LIBS)

nodist__pixel_la_SOURCES = _pixel.cpp
_pixel_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_pixel_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_pixel_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_VW_LIBS)

nodist__pixelcast_la_SOURCES = _pixelcast.cpp
_pixelcast_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_pixelcast_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_pixelcast_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_VW_LIBS)

nodist__qtree_la_SOURCES = _qtree.cpp
_qtree_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_qtree_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_qtree_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_MOSAIC_LIBS)

nodist__transform_la_SOURCES = _transform.cpp
_transform_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_transform_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_transform_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_VW_LIBS)

nodist__vwmath_la_SOURCES = _vwmath.cpp
_vwmath_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(VW_CPPFLAGS) $(NUMPY_CPPFLAGS)
_vwmath_la_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module $(VW_LDFLAGS)
_vwmath_la_LIBADD = $(PYTHON_EXTRA_LIBS) $(PKG_VW_LIBS)

python_vwdir = ${pythondir}/vw
pyexec_vwdir = ${pyexecdir}/vw

python_vw_PYTHON = __init__.py
nodist_python_vw_PYTHON = cartography.py composite.py core.py fileio.py filter.py image.py imagealgo.py imagemanip.py imagemath.py pixel.py pixelcast.py qtree.py transform.py vwmath.py
pyexec_vw_LTLIBRARIES = _cartography.la _composite.la _core.la _fileio.la _filter.la _image.la _imagealgo.la _imagemanip.la _imagemath.la _pixel.la _pixelcast.la _qtree.la _transform.la _vwmath.la

.i.cpp:
	$(SWIG) $(SWIG_PYTHON_OPT) $(SWIG_PYTHON_CPPFLAGS) $(VW_CPPFLAGS) -o $@ $<

endif
