cmake_minimum_required (VERSION 3.2)

# Make it possible to append to these from the command line
set(CMAKE_CXX_FLAGS "" CACHE STRING "")
set(CMAKE_C_FLAGS "" CACHE STRING "")

# Tell cmake to look in the /cmake folder.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

## Add precompiled header tool
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/CMakePCHCompiler)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
set(CMAKE_EXE_LINKER_FLAGS "-L${BINARYBUILDER_INSTALL_DIR}/lib -lm -pthread -lpthread")

set(CMAKE_BUILD_TYPE "Release")
message("Setting CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

add_subdirectory(src)

