# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

lib_LTLIBRARIES =

if MAKE_MODULE_PLATE

# This is necessary because otherwise automake doesn't see vwPlate as a
# dependency for compile-order calculations
PLATE_LOCAL_LIBS = libvwPlate.la @MODULE_PLATE_LIBS@

protocol_headers = IndexService.pb.h Rpc.pb.h
protocol_sources = IndexService.pb.cc Rpc.pb.cc

BUILT_SOURCES = $(protocol_sources)

include_HEADERS =           \
  AmqpConnection.h          \
  Blob.h                    \
  BlobManager.h             \
  common.h                  \
  Exception.h               \
  HTTPUtils.h               \
  Index.h                   \
  IndexPage.h               \
  IndexService.h            \
  LocalIndex.h              \
  PagedIndex.h              \
  PlateCarreePlateManager.h \
  PlateFile.h               \
  PlateManager.h            \
  PlateView.h               \
  PolarStereoPlateManager.h \
  RemoteIndex.h             \
  RpcServices.h             \
  SnapshotManager.h         \
  TileManipulation.h        \
  ToastDem.h                \
  ToastPlateManager.h
#  WMSService.h

nodist_include_HEADERS = $(protocol_headers)
noinst_HEADERS = mod_plate.h mod_plate_utils.h mod_plate_core.h mod_plate_handlers.h

libvwPlate_la_SOURCES =      \
  AmqpConnection.cc          \
  Blob.cc                    \
  BlobManager.cc             \
  HTTPUtils.cc               \
  Index.cc                   \
  IndexPage.cc               \
  IndexService.cc            \
  LocalIndex.cc              \
  PagedIndex.cc              \
  PlateCarreePlateManager.cc \
  PlateFile.cc               \
  PlateManager.cc            \
  PolarStereoPlateManager.cc \
  RemoteIndex.cc             \
  RpcServices.cc             \
  SnapshotManager.cc         \
  TileManipulation.cc        \
  ToastDem.cc                \
  ToastPlateManager.cc
#  WMSService.cc
nodist_libvwPlate_la_SOURCES = $(protocol_sources)

libvwPlate_la_LIBADD  = @MODULE_PLATE_LIBS@
libvwPlate_la_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

lib_LTLIBRARIES += libvwPlate.la

endif

# ----------------
#  mod_plate
# ----------------
if MAKE_MODPLATE

mod_plate_la_SOURCES = mod_plate.c mod_plate_handlers.cc mod_plate_utils.cc mod_plate_core.cc
mod_plate_la_CPPFLAGS = @APXS_CFLAGS@ @VW_CPPFLAGS@
mod_plate_la_LIBADD = $(PLATE_LOCAL_LIBS)
mod_plate_la_LDFLAGS = @APXS_LDFLAGS@ -module -rpath $(APXS_INSTALLDIR) -avoid-version
noinst_LTLIBRARIES = mod_plate.la

install-exec-local:
	@if ! $(LIBTOOL) --mode=install cp mod_plate.la $(APXS_INSTALLDIR); then\
		echo "********************************************************************************";\
		echo "*                                                                              *";\
		echo "*                                                                              *";\
		echo "*                    WARNING: Failed to install mod_plate!                     *";\
		echo "*                                                                              *";\
		echo "*                                                                              *";\
		echo "********************************************************************************";\
	fi
	rm -f $(APXS_INSTALLDIR)/mod_plate.la

apxs_install: mod_plate.la
	$(LIBTOOL) --mode=install cp mod_plate.la $(APXS_INSTALLDIR)
	rm -f $(APXS_INSTALLDIR)/mod_plate.la

endif

if ENABLE_EXCEPTIONS

index_server_SOURCES = index_server.cc
index_server_LDADD   = $(PLATE_LOCAL_LIBS)

#wms_server_SOURCES = wms_server.cc
#wms_server_LDADD = $(PLATE_LOCAL_LIBS)
#
#wms_client_SOURCES = wms_client.cc
#wms_client_LDADD = $(PLATE_LOCAL_LIBS)

image2plate_SOURCES = image2plate.cc
image2plate_LDADD   = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ $(PLATE_LOCAL_LIBS)

snapshot_SOURCES = snapshot.cc
snapshot_LDADD   = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ $(PLATE_LOCAL_LIBS)

rebuild_index_SOURCES = rebuild_index.cc
rebuild_index_LDADD   = @PKG_CARTOGRAPHY_LIBS@ @PKG_MOSAIC_LIBS@ $(PLATE_LOCAL_LIBS)

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

rpc_tool_SOURCES = rpc_tool.cc
rpc_tool_LDADD = $(PLATE_LOCAL_LIBS)

lustre_torture_SOURCES = lustre_torture.cc
lustre_torture_LDADD = $(PLATE_LOCAL_LIBS)

hirise2tif_SOURCES = hirise2tif.cc
hirise2tif_LDADD   = @PKG_CARTOGRAPHY_LIBS@ $(PLATE_LOCAL_LIBS)

plate2plate_SOURCES = plate2plate.cc
plate2plate_LDADD   = @PKG_CARTOGRAPHY_LIBS@ $(PLATE_LOCAL_LIBS)

platereduce_SOURCES = platereduce.cc
platereduce_LDADD   = $(PLATE_LOCAL_LIBS)

bin_PROGRAMS =   \
  amqp_perftest  \
  hirise2tif     \
  image2plate    \
  index_client   \
  index_perftest \
  index_server   \
  lustre_torture \
  plate2dem      \
  plate2plate    \
  plate2tiles    \
  platereduce    \
  rebuild_index  \
  rpc_tool       \
  snapshot       \
  tiles2plate
  #wms_client
  #wms_server

endif


########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@

SUBDIRS = . tests google
CLEANFILES = $(protocol_sources) $(protocol_headers)

includedir = $(prefix)/include/vw/Plate

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/thirdparty/protobuf.mak
