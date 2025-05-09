cmake_minimum_required(VERSION 3.15.5)

project(visionworkbench)

if(APPLE)
    cmake_policy(SET CMP0025 NEW)
    set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13)
endif(APPLE)

# Make it possible to append to these from the command line
set(CMAKE_CXX_FLAGS "" CACHE STRING "")
set(CMAKE_C_FLAGS "" CACHE STRING "")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++14")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_EXE_LINKER_FLAGS "-pthread -lpthread")

set(USE_OPENEXR ON CACHE BOOL "Use EXR (default ON)")

if (APPLE)
    # A workaround for the clang included with conda build
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mlinker-version=305")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mlinker-version=305")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mlinker-version=305")
    # Prevent a pedantic error in recent clang
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-enum-constexpr-conversion")
else()
    # On Linux need to link to additional libraries
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -lrt -lgcc_s")
endif()

# Make sure we do a release type
set(CMAKE_BUILD_TYPE "Release")
message("Setting CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

# Tell cmake to look in the cmake folder
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

add_subdirectory(src)
