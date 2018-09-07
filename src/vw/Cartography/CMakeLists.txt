

# Build all protobuf files and list names in this variable
#generate_protobuf_files(PROTO_GEN)

# Use wrapper function at this level to avoid code duplication
add_library_wrapper(VwCartography "${VW_CARTOGRAPHY_SRC_FILES}" "${VW_CARTOGRAPHY_TEST_FILES}" "${VW_CARTOGRAPHY_LIB_DEPENDENCIES}")


