cmake_minimum_required (VERSION 2.8)

# Tell cmake to look in the /cmake folder.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

message("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")
add_subdirectory(src)

# TODO: Precompiled headers!

# TODO: Flood support!

# TODO: Better conditional building!

# TODO: Streamline library dependencies!

# TODO: Should probably do more at this top level

# TODO: Compare ninja speed to make
