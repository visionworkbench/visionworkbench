########################################################################
# rules for compiling Google Protocol Buffers (using protoc)
########################################################################

SUFFIXES += .proto .pb.cc

PROTOC_ARGS =

# This bit of pwd uglyness is because autoconf 2.59 & automake 1.9.6 don't
# canonicalize abs_top_srcdir... instead it ends up being relative to
# $abs_srcdir (e.g. $abs_srcdir/../../). Protobuf can't handle relative paths.

.proto.pb.cc :
	$(AM_V_GEN)( \
		SRC=`cd $(abs_top_srcdir)/src && pwd` ;\
		OBJ=`cd $(abs_top_builddir)/src && pwd` ;\
		FILEDIR=`cd $$(dirname $<) && pwd` ;\
		FILENAME="$$FILEDIR/$$(basename $<)" ;\
		cd "$$SRC" && \
		$(PROTOC) -I"$$SRC" -I"$$OBJ" --cpp_out="$$OBJ" $(PROTOC_ARGS) "$$FILENAME" \
	)
