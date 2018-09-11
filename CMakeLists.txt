cmake_minimum_required (VERSION 3.2)

# Tell cmake to look in the /cmake folder.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

## Add precompiled header tool
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/CMakePCHCompiler)

add_subdirectory(src)

# TODO: Precompiled headers!
#       - Wait for response to issue

# TODO: Build on the Mac

# TODO: Should probably do more at this top level

# TODO: Compare ninja speed to make

