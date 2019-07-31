cmake_minimum_required (VERSION 3.2)

# Make it possible to append to these from the command line
set(CMAKE_CXX_FLAGS "" CACHE STRING "")
set(CMAKE_C_FLAGS "" CACHE STRING "")

# Tell cmake to look in the /cmake folder.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

## Add precompiled header tool
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/CMakePCHCompiler)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 -lm")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 -lm")

add_subdirectory(src)

# TODO: Precompiled headers!
#       - Wait for response to issue

# TODO: Build on the Mac

# TODO: Should probably do more at this top level

# TODO: Compare ninja speed to make

