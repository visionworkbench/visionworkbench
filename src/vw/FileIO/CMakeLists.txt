# Use wrapper function at this level to avoid code duplication
add_library_wrapper(VwFileIO "${VW_FILEIO_SRC_FILES}" "${VW_FILEIO_TEST_FILES}" "${VW_FILEIO_LIB_DEPENDENCIES}")
