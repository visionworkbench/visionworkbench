if(APPLE)
    # Needed for conda build on OSX with conda
    cmake_policy(SET CMP0025 NEW) 
endif(APPLE)

cmake_minimum_required (VERSION 3.15)

# Make it possible to append to these from the command line
set(CMAKE_CXX_FLAGS "" CACHE STRING "")
set(CMAKE_C_FLAGS "" CACHE STRING "")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_EXE_LINKER_FLAGS "-pthread -lpthread")

if (APPLE)
    # A workaround for the clang included with conda build
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mlinker-version=305")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mlinker-version=305")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mlinker-version=305")
else()
    # On Linux need to link to additional libraries
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -lrt -lgcc_s")
endif()

# Tell cmake to look in the /cmake folder.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Make sure we do a release type
set(CMAKE_BUILD_TYPE "Release")
message("Setting CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

project(visionworkbench)

add_subdirectory(src)
