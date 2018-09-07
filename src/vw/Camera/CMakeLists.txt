 
## Build all protobuf files and list names in this variable
#generate_protobuf_files(PROTO_GEN)

# Use wrapper function at this level to avoid code duplication
add_library_wrapper(VwCamera "${VW_CAMERA_SRC_FILES}" "${VW_CAMERA_TEST_FILES}" "${VW_CAMERA_LIB_DEPENDENCIES}")


