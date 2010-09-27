########################################################################
# rules for compiling Google Protocol Buffers (using protoc)
########################################################################

SUFFIXES += .proto .pb.cc

PROTOC = protoc
PROTOC_ARGS =

.proto.pb.cc :
	$(AM_V_GEN)oldpwd=`pwd` && ( \
		cd $(abs_top_srcdir)/src && \
		$(PROTOC) -I$(abs_top_srcdir)/src -I$(abs_top_builddir)/src --cpp_out="$(abs_top_builddir)/src" $(PROTOC_ARGS) "$(abs_srcdir)/"$< \
	)
