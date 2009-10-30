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

protocol_headers = ProtoBuffers.pb.h 
protocol_sources = ProtoBuffers.pb.cc

include_HEADERS = Blob.h BlobManager.h Index.h RemoteIndex.h Tree.h	\
		  PlateFile.h Amqp.h PlateManager.h			\
		  ToastPlateManager.h KmlPlateManager.h			\
		  PlateView.h
nodist_include_HEADERS = $(protocol_headers)

libvwPlate_la_SOURCES = Blob.cc BlobManager.cc	\
			Index.cc RemoteIndex.cc Amqp.cc			\
			ToastPlateManager.cc KmlPlateManager.cc
nodist_libvwPlate_la_SOURCES = $(protocol_sources)

libvwPlate_la_LIBADD = @MODULE_PLATE_LIBS@

lib_LTLIBRARIES += libvwPlate.la 

endif 

# ----------------
#  mod_plate
# ----------------
#mod_plate_io.o:
#	$(CXX) -I/opt/local/include/apr-1 -I$(APXS_INCLUDEDIR) @VW_CPPFLAGS@ -c mod_plate_io.cc


if MAKE_MODPLATE

mod_plate_la_SOURCES = mod_plate.c mod_plate_io.cc
mod_plate_la_CPPFLAGS = -I/opt/local/include/apr-1 -I$(APXS_INCLUDEDIR) @VW_CPPFLAGS@
mod_plate_la_LDFLAGS = @MODULE_PLATE_LIBS@ -module
lib_LTLIBRARIES += mod_plate.la

apxs_install:
	sudo $(APXS) -i mod_plate.la

endif

if ENABLE_EXCEPTIONS

mosaic_worker_SOURCES = mosaic_worker.cc
mosaic_worker_LDADD    = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ @PKG_PLATE_LIBS@ 

plateindex_server_SOURCES = plateindex_server.cc
plateindex_server_LDADD    = @PKG_PLATE_LIBS@ libvwPlate.la

image2plate_SOURCES = image2plate.cc
image2plate_LDADD    = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ @PKG_PLATE_LIBS@ 

plate2tiles_SOURCES = plate2tiles.cc
plate2tiles_LDADD    = @PKG_PLATE_LIBS@ 

plate2dem_SOURCES = plate2dem.cc
plate2dem_LDADD    = @PKG_CARTOGRAPHY_LIBS@ @PKG_PLATE_LIBS@ 

bin_PROGRAMS = plate2dem plate2tiles image2plate  mosaic_worker plateindex_server
noinst_PROGRAMS = 

endif


########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ -version-info @LIBTOOL_VERSION@

SUBDIRS = . tests
CLEANFILES = $(protocol_sources) $(protocol_headers)

includedir = $(prefix)/include/vw/Plate

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/thirdparty/protobuf.mak
