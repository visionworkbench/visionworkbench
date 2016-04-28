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

if MAKE_MODULE_CARTOGRAPHY

if HAVE_PKG_GDAL
gdal_headers = GeoReferenceResourceGDAL.h
gdal_sources = GeoReferenceResourceGDAL.cc
endif

if MAKE_MODULE_CAMERA
camerabbox_headers = CameraBBox.h
camerabbox_sources = CameraBBox.cc
endif


include_HEADERS = GeoReference.h  GeoReferenceUtils.h                  \
                  GeoTransform.h Datum.h SimplePointImageManipulation.h   \
                  PointImageManipulation.h Map2CamTrans.h                 \
                  OrthoImageView.h GeoReferenceResourcePDS.h              \
                  Projection.h ToastTransform.h Chipper.h $(gdal_headers) \
                  $(camerabbox_headers)


#detail_HEADERS =  detail/BresenhamLine.h

libvwCartography_la_SOURCES = Datum.cc GeoReference.cc GeoTransform.cc  \
                  GeoReferenceResourcePDS.cc ToastTransform.cc          \
                  PointImageManipulation.cc GeoReferenceUtils.cc        \
                  Map2CamTrans.cc Chipper.cc $(gdal_sources)            \
                  $(camerabbox_sources)

nodist_libvwCartography_la_SOURCES = 

libvwCartography_la_LIBADD = @MODULE_CARTOGRAPHY_LIBS@

lib_LTLIBRARIES = libvwCartography.la

endif

########################################################################
# general
########################################################################

BUILT_SOURCES = 

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests
CLEANFILES = 
EXTRA_DIST = 

includedir = $(prefix)/include/vw/Cartography
#detaildir = $(includedir)/detail

include $(top_srcdir)/config/rules.mak
