# Add path for custom modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

include(Utilities) # Contains functions

# Use CCache if it is available
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  #message("Using ccache tool...")
  #set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  #set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

if (UNIX AND NOT APPLE)
  execute_process(COMMAND ${CMAKE_C_COMPILER} -fuse-ld=gold -Wl,--version ERROR_QUIET OUTPUT_VARIABLE ld_version)
  if ("${ld_version}" MATCHES "GNU gold")
    message("Using gold linker...")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=gold -Wl,--disable-new-dtags")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=gold -Wl,--disable-new-dtags")
  endif()
endif()


#===============================================================================
# Configuration options

set(VW_ENABLE_HDR "0" CACHE BOOL "Build the HDR module.")
set(VW_NUM_THREADS "4" CACHE STRING "The default number of threads to use.")

if("${CUSTOM_BUILD_DIR}" STREQUAL "")
    set(CUSTOM_BUILD_DIR "${ASP_DEPS_DIR}")
endif()
    
if("${ASP_DEPS_DIR}" STREQUAL "")
  message(FATAL_ERROR "You need to set ASP_DEPS_DIR.")
endif()

if("${CUSTOM_BUILD_DIR}" STREQUAL "")
    set(CUSTOM_BUILD_DIR "${ASP_DEPS_DIR}")
endif()

# Where to search for dependencies
set(CMAKE_PREFIX_PATH "${ASP_DEPS_DIR}")

# Generated protocol buffers go here
include_directories(${CMAKE_BINARY_DIR}/src) 

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
# foreach(dir ${dirs})
#   message(STATUS "dir='${dir}'")
# endforeach()

# Default this to empty
SET_PROPERTY(GLOBAL PROPERTY storedPchFirstLibrary "NA")

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Make sure installed files can find installed libraries
if (APPLE)
  SET(CMAKE_INSTALL_RPATH "@loader_path/../lib" ${CMAKE_INSTALL_RPATH})
elseif (UNIX)
  SET(CMAKE_INSTALL_RPATH "$ORIGIN/../lib" ${CMAKE_INSTALL_RPATH})
endif()
SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} ${ASP_DEPS_DIR}/lib)

# --- Custom options ---
option(BUILD_SHARED_LIBS "Produce shared libraries." TRUE)

# --- Fixed options ---
set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON )

# Enable/disable SSE optimizations
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
  set(VW_ENABLE_SSE 0)
else()
  set(VW_ENABLE_SSE 1)
endif()
if (VW_ENABLE_SSE)
  set(CXXFLAGS "${CXXFLAGS} -msse4.1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
  message(STATUS "Enabling SSE")
else()
  if (CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
    message(STATUS "SSE optimizations are not available on Mac Arm.")
  else()
    # Older Intel CPUs do not support SSE4.1
    message(STATUS, "Disabling SSE.")
    set(CXXFLAGS "${CXXFLAGS} -mno-sse4.1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-sse4.1")
  endif()
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#==================================================================================
# Find dependencies

# assume built-in pthreads on MacOS
# IF(APPLE)
#     set(CMAKE_THREAD_LIBS_INIT "-lpthread")
#     set(CMAKE_HAVE_THREADS_LIBRARY 1)
#     set(CMAKE_USE_WIN32_THREADS_INIT 0)
#     set(CMAKE_USE_PTHREADS_INIT 1)
#     set(THREADS_PREFER_PTHREAD_FLAG ON)
# ENDIF()
    
find_package(Threads REQUIRED)

message("--------------------------------------------------------------")
message("Searching for boost...")

set(REQUIRED_BOOST_LIBS program_options system filesystem regex date_time thread iostreams)
#set(Boost_DEBUG 1) # Uncomment for useful finding boost debug output
if(ASP_DEPS_DIR) # Look in the directory of ASP dependencies
  set(Boost_ROOT       "${ASP_DEPS_DIR}")
  set(Boost_LIBRARYDIR "${ASP_DEPS_DIR}/lib")
  set(Boost_INCLUDEDIR "${ASP_DEPS_DIR}/include")
  set(Boost_NO_SYSTEM_PATHS ON) # Do not search anywhere else
  set(Boost_DEBUG ON)
  set(Boost_DETAILED_FAILURE_MSG ON)
  set(CMAKE_CXX_COMPILER_ARCHITECTURE_ID "x64") # needed for pfe
endif()  
# Now search for Boost using the available information
find_package(Boost REQUIRED COMPONENTS ${REQUIRED_BOOST_LIBS})

if(Boost_FOUND)
    message("Successfully found Boost")
    message("Boost include dir = ${Boost_INCLUDE_DIR}")
    message("Boost library dir = ${Boost_LIBRARY_DIRS}")
    message("Boost libraries   = ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIRS})
endif(Boost_FOUND)
message("--------------------------------------------------------------")

set(REQUIRED     True )
set(NOT_REQUIRED False)

find_external_library(OPENBLAS ${ASP_DEPS_DIR} "" "openblas" "" REQUIRED) # For arm
find_external_library(LAPACK ${ASP_DEPS_DIR} "" "lapack;blas" "" REQUIRED)
set(VW_HAVE_PKG_LAPACK 1)

find_external_library(FLANN   ${ASP_DEPS_DIR} "" "flann_cpp" REQUIRED)
find_external_library(LZ4     ${ASP_DEPS_DIR} "" "lz4"       REQUIRED)
find_external_library(JPEG    ${ASP_DEPS_DIR} "" "jpeg"      REQUIRED)
find_external_library(PNG     ${ASP_DEPS_DIR} "" "png;png16" REQUIRED)
find_external_library(TIFF    ${ASP_DEPS_DIR} "" "tiff"      REQUIRED)
#find_external_library(GEOTIFF ${ASP_DEPS_DIR} "" "geotiff"   REQUIRED)
find_external_library(Z       ${ASP_DEPS_DIR} "" "z"         REQUIRED)
find_external_library(PROJ4   ${ASP_DEPS_DIR} "" "proj"      REQUIRED)
find_external_library(GDAL    ${CUSTOM_BUILD_DIR} "" "gdal"  REQUIRED)

# TODO: Make a function?
set(OPENCV_LIB_NAMES opencv_calib3d     opencv_reg
                     opencv_core        opencv_shape
                     opencv_features2d  opencv_stitching
                     opencv_flann       opencv_superres
                     opencv_hdf         opencv_surface_matching
                     opencv_highgui     opencv_videoio
                     opencv_imgcodecs   opencv_video
                     opencv_imgproc     opencv_xfeatures2d
                     opencv_ml          opencv_ximgproc
                     opencv_objdetect   opencv_xobjdetect
                     opencv_photo       opencv_xphoto)
find_external_library(OPENCV ${ASP_DEPS_DIR} "opencv4" "${OPENCV_LIB_NAMES}" REQUIRED)

if (USE_OPENEXR)
  set(ILMBASE_LIB_NAMES Half Iex Imath IlmThread)
  find_external_library(ILMBASE ${ASP_DEPS_DIR} "" "${ILMBASE_LIB_NAMES}" REQUIRED)
  find_external_library(OPENEXR ${ASP_DEPS_DIR} "OpenEXR" "IlmImf;IlmThread" REQUIRED)
  set(VW_HAVE_PKG_OPENEXR 1)
else()
  message(STATUS "OpenEXR will not be used.")
  set(VW_HAVE_PKG_OPENEXR 0)
endif()

# CMake will try to find these packages if needed. They are not part
# of the usual dependencies.
if (VW_ENABLE_PYTHON)
  find_package(PythonLibs)
  if(PYTHONLIBS_FOUND)
    set(VW_HAVE_PKG_PYTHONLIBS 1)
    message("-- Found PythonLibs " ${PYTHON_LIBRARIES})
  else()
    message("Failed to find PythonLibs")
  endif()
endif()

###########################################################################
# These are some options that the user might like to set manually

set(VW_DEFAULT_CACHE_SIZE_MB 768)

# enable image bounds checking (SLOW!) 
set(VW_ENABLE_BOUNDS_CHECK 0)

# ~/.vwrc support
set(VW_ENABLE_CONFIG_FILE 1) 

# TODO: Remove this as an option!
# enable the C++ exception mechanism 
set(VW_ENABLE_EXCEPTIONS 1) 

# Until it is needed we don't support this in CMake.
set(VW_ENABLE_PYTHON 0)


# Define to `int' if <sys/types.h> does not define. 
#define VW_ssize_t

# TODO: What is this for?
set(ENABLE_EXCEPTIONS True)

###########################################################################

# Libraries use paths like vw/Core/FundamentalTypes.h, so we need to
# include the top of the tree.  These must take precedence over
# anything in the conda or system directories. That is why this is
# towards the bottom.
include_directories(BEFORE ../)

# Now that we have found all our external dependencies, generate a config.h file
include("GenerateConfig.cmake")

message("--------------------------------------------------------------")
#message("Listing all include directories...")
#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()
#message("--------------------------------------------------------------")

# Add all the header files at the top level to the install command
get_all_source_files("." MAIN_HEADER_FILES)
foreach(f ${MAIN_HEADER_FILES})
  INSTALL(FILES ${f} DESTINATION include/vw)
endforeach()

# Install a small tif file for testing when conda-forge builds the package
install(FILES Cartography/tests/tinyDemAN.tif DESTINATION share/vw/tests)

#=========================================================================
# Define each of the libaries

# --- VW_BUNDLEADJUSTMENT ------------------------------------------------------------
get_all_source_files( "BundleAdjustment"       VW_BUNDLEADJUSTMENT_SRC_FILES)
get_all_source_files( "BundleAdjustment/tests" VW_BUNDLEADJUSTMENT_TEST_FILES)
set(VW_BUNDLEADJUSTMENT_LIB_DEPENDENCIES VwMosaic VwFileIO VwInterestPoint VwStereo VwCartography)

# --- VW_CAMERA ------------------------------------------------------------
get_all_source_files( "Camera"       VW_CAMERA_SRC_FILES)
get_all_source_files( "Camera/tests" VW_CAMERA_TEST_FILES)
set(VW_CAMERA_LIB_DEPENDENCIES VwMath VwImage VwFileIO)

# --- VW_CORE ------------------------------------------------------------
get_all_source_files( "Core"       VW_CORE_SRC_FILES)
get_all_source_files( "Core/tests" VW_CORE_TEST_FILES)
set(VW_CORE_LIB_DEPENDENCIES ${Boost_LIBRARIES} pthread)

# --- VW_FILEIO ------------------------------------------------------------
get_all_source_files("FileIO" VW_FILEIO_SRC_FILES)
get_all_source_files("FileIO/tests" VW_FILEIO_TEST_FILES)
set(VW_FILEIO_LIB_DEPENDENCIES VwMath VwImage ${GDAL_LIBRARIES} ${PROJ4_LIBRARIES}
    ${JPEG_LIBRARIES} ${TIFF_LIBRARIES} ${PNG_LIBRARIES} 
    ${ILMBASE_LIBRARIES} ${Z_LIBRARIES})
if (USE_OPENEXR)
  list(APPEND VW_FILEIO_LIB_DEPENDENCIES ${OPENEXR_LIBRARIES})
endif()

# --- VW_GEOMETRY ------------------------------------------------------------
get_all_source_files( "Geometry"       VW_GEOMETRY_SRC_FILES)
get_all_source_files( "Geometry/tests" VW_GEOMETRY_TEST_FILES)
set(VW_GEOMETRY_LIB_DEPENDENCIES VwCore VwMath)

# --- VW_CARTOGRAPHY ------------------------------------------------------------
get_all_source_files( "Cartography"       VW_CARTOGRAPHY_SRC_FILES)
get_all_source_files( "Cartography/tests" VW_CARTOGRAPHY_TEST_FILES)
set(VW_CARTOGRAPHY_LIB_DEPENDENCIES VwFileIO VwCamera VwGeometry
    ${GDAL_LIBRARIES} ${PROJ4_LIBRARIES}
    ${JPEG_LIBRARIES} ${TIFF_LIBRARIES} ${PNG_LIBRARIES})

# --- VW_HDR ------------------------------------------------------------
get_all_source_files( "HDR"       VW_HDR_SRC_FILES)
get_all_source_files( "HDR/tests" VW_HDR_TEST_FILES)
set(VW_HDR_LIB_DEPENDENCIES VwFileIO VwCamera )

# --- VW_IMAGE ------------------------------------------------------------
get_all_source_files( "Image"       VW_IMAGE_SRC_FILES)
get_all_source_files( "Image/tests" VW_IMAGE_TEST_FILES)
set(VW_IMAGE_LIB_DEPENDENCIES VwMath ${OPENCV_LIBRARIES})

# --- VW_INTERESTPOINT ------------------------------------------------------------
get_all_source_files( "InterestPoint"       VW_INTERESTPOINT_SRC_FILES)
get_all_source_files( "InterestPoint/tests" VW_INTERESTPOINT_TEST_FILES)
set(VW_INTERESTPOINT_LIB_DEPENDENCIES VwCamera VwFileIO)

# --- VW_MATH ------------------------------------------------------------
get_all_source_files( "Math"       VW_MATH_SRC_FILES)
get_all_source_files( "Math/tests" VW_MATH_TEST_FILES)
set(VW_MATH_LIB_DEPENDENCIES VwCore ${FLANN_LIBRARIES} ${LAPACK_LIBRARIES} ${LZ4_LIBRARIES})

# --- VW_MOSAIC ------------------------------------------------------------
get_all_source_files( "Mosaic"       VW_MOSAIC_SRC_FILES)
get_all_source_files( "Mosaic/tests" VW_MOSAIC_TEST_FILES)
set(VW_MOSAIC_LIB_DEPENDENCIES VwFileIO VwCamera VwCartography)

# --- VW_STEREO ------------------------------------------------------------
get_all_source_files( "Stereo"       VW_STEREO_SRC_FILES)
get_all_source_files( "Stereo/tests" VW_STEREO_TEST_FILES)
set(VW_STEREO_LIB_DEPENDENCIES VwMath VwImage VwCamera VwGeometry)

#=========================================================================
# Add all of the library subdirectories

add_subdirectory(Core)
add_subdirectory(Math)
add_subdirectory(Geometry)
add_subdirectory(FileIO)
add_subdirectory(Camera)
add_subdirectory(Cartography)
add_subdirectory(Image)
add_subdirectory(InterestPoint)
add_subdirectory(Mosaic)
add_subdirectory(Stereo)
add_subdirectory(BundleAdjustment)
if (${VW_ENABLE_HDR})
  add_subdirectory(HDR)
endif()
if (${VW_ENABLE_PYTHON})
  add_subdirectory(Python)
endif()

# Add the executable subdirectory
add_subdirectory(tools)

