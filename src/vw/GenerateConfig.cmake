# Generate vw_config.h from vw_config.h.in, and vw_date_config.h, which
# changes often, from vw_date_config.h.in.

message("Generating config file: ${CMAKE_BINARY_DIR}/src/vw/vw_config.h")
configure_file(${CMAKE_SOURCE_DIR}/src/vw/vw_config.h.in ${CMAKE_BINARY_DIR}/src/vw/vw_config.h)

message("Generating config file: ${CMAKE_BINARY_DIR}/src/vw/vw_date_config.h")
configure_file(${CMAKE_SOURCE_DIR}/src/vw/vw_date_config.h.in ${CMAKE_BINARY_DIR}/src/vw/vw_date_config.h)
