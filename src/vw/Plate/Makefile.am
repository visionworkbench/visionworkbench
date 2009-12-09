# __BEGIN_LICENSE__
# Copyright (C) 2006-2009 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

lib_LTLIBRARIES =

if MAKE_MODULE_PLATE

# This is necessary because otherwise automake doesn't see vwPlate as a
# dependency for compiler-order calculations
PLATE_LOCAL_LIBS = libvwPlate.la @MODULE_PLATE_LIBS@

protocol_headers = ProtoBuffers.pb.h
protocol_sources = ProtoBuffers.pb.cc

BUILT_SOURCES = $(protocol_sources)

include_HEADERS = Exception.h Blob.h BlobManager.h Index.h		\
		  LocalIndex.h RemoteIndex.h Tree.h PlateFile.h Amqp.h	\
		  PlateManager.h ToastPlateManager.h KmlPlateManager.h	\
		  PlateView.h RpcServices.h IndexService.h
nodist_include_HEADERS = $(protocol_headers)

libvwPlate_la_SOURCES = Blob.cc BlobManager.cc Index.cc LocalIndex.cc	\
			RemoteIndex.cc PlateFile.cc RpcServices.cc	\
			IndexService.cc Amqp.cc PlateManager.cc		\
			ToastPlateManager.cc
nodist_libvwPlate_la_SOURCES = $(protocol_sources)

libvwPlate_la_LIBADD  = @MODULE_PLATE_LIBS@
libvwPlate_la_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

lib_LTLIBRARIES += libvwPlate.la

endif

# ----------------
#  mod_plate
# ----------------
if MAKE_MODPLATE

mod_plate_la_SOURCES = mod_plate.c mod_plate_io.cc
mod_plate_la_CPPFLAGS = @APXS_CFLAGS@ @VW_CPPFLAGS@
mod_plate_la_LIBADD = $(PLATE_LOCAL_LIBS)
mod_plate_la_LDFLAGS = @APXS_LDFLAGS@ -module -rpath $(APXS_INSTALLDIR) -avoid-version
noinst_LTLIBRARIES = mod_plate.la

apxs_install:
	$(APXS) -i mod_plate.la

endif

if ENABLE_EXCEPTIONS

index_server_SOURCES = index_server.cc
index_server_LDADD   = $(PLATE_LOCAL_LIBS)

image2plate_SOURCES = image2plate.cc
image2plate_LDADD   = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ $(PLATE_LOCAL_LIBS)

tiles2plate_SOURCES = tiles2plate.cc
tiles2plate_LDADD = $(PLATE_LOCAL_LIBS)

plate2tiles_SOURCES = plate2tiles.cc
plate2tiles_LDADD    = $(PLATE_LOCAL_LIBS)

plate2dem_SOURCES = plate2dem.cc
plate2dem_LDADD    = @PKG_CARTOGRAPHY_LIBS@ $(PLATE_LOCAL_LIBS)

index_client_SOURCES = index_client.cc
index_client_LDADD = $(PLATE_LOCAL_LIBS)

amqp_perftest_SOURCES = amqp_perftest.cc
amqp_perftest_LDADD = $(PLATE_LOCAL_LIBS)

index_perftest_SOURCES = index_perftest.cc
index_perftest_LDADD = $(PLATE_LOCAL_LIBS)

bin_PROGRAMS = plate2dem plate2tiles tiles2plate image2plate	\
	       index_server index_client amqp_perftest index_perftest

endif


########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@

SUBDIRS = . tests
CLEANFILES = $(protocol_sources) $(protocol_headers)

includedir = $(prefix)/include/vw/Plate

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/thirdparty/protobuf.mak
