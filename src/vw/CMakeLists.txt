# Add path for custom modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

include(Utilities) # Contains functions

# Use CCache if it is available.
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

# TODO: Change some of these to options?
set(BINARYBUILDER_INSTALL_DIR "" CACHE FILEPATH "Path to Binary Builder installation directory")
set(VW_ENABLE_SSE "1" CACHE BOOL "Allow SSE optimizations.")
set(VW_ENABLE_HDR "0" CACHE BOOL "Build the HDR module.")
set(VW_NUM_THREADS "4" CACHE STRING "The default number of threads to use.")

# --- Directory paths ---
if(BINARYBUILDER_INSTALL_DIR)
  message("Looking for dependencies in BinaryBuilder directory: ${BINARYBUILDER_INSTALL_DIR}") 
  # Add BINARYBUILDER folders to CMAKE paths
  set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "${BINARYBUILDER_INSTALL_DIR}")
  set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "${BINARYBUILDER_INSTALL_DIR}")
endif()

if("${ISIS_DEPS_DIR}" STREQUAL "")
  message(FATAL_ERROR "You need to set ISIS_DEPS_DIR: ${ISIS_DEPS_DIR}")
endif()

# Libraries use paths like "vw/blah/blah.h" so we need to include the top of the tree
include_directories( ../ )
include_directories( ${CMAKE_BINARY_DIR}/src/ ) # Generated protocol buffers go here!

# Look for libraries in both of these BB locations.
#link_directories(${BINARYBUILDER_INSTALL_DIR}/lib)
#link_directories(${BINARYBUILDER_INSTALL_DIR}/lib64)

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
SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} ${CMAKE_INSTALL_PREFIX}/lib ${ISIS_DEPS_DIR}/lib ${ISIS_DEPS_DIR}/lib64)


# --- Custom options ---
option(BUILD_SHARED_LIBS "Produce shared libraries." TRUE)

# --- Fixed options ---
set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#==================================================================================
# Find dependencies

# assume built-in pthreads on MacOS
IF(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
ENDIF()
    
find_package(Threads REQUIRED)

message("--------------------------------------------------------------")
message("Searching for boost...")

# TODO: Make into a function!
# Define boost requirements
set(BOOST_MIN_VERSION "1.65")
set(REQUIRED_BOOST_LIBS program_options system filesystem regex date_time thread iostreams)
#set(Boost_DEBUG 1) # Uncomment for useful finding boost debug output
#if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
if(ISIS_DEPS_DIR) # Look in BinaryBuilder directory
  set(BOOST_ROOT       "${ISIS_DEPS_DIR}")
  set(BOOST_LIBRARYDIR "${ISIS_DEPS_DIR}/lib")
  set(BOOST_INCLUDEDIR "${ISIS_DEPS_DIR}/include")

  set(Boost_DEBUG ON)
  set(Boost_DETAILED_FAILURE_MSG ON)
  set(CMAKE_CXX_COMPILER_ARCHITECTURE_ID "x64") # needed for pfe
  #set(Boost_NO_SYSTEM_PATHS ON) #This seems to confuse boost
endif()  
# Now search for Boost using the available information
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS ${REQUIRED_BOOST_LIBS})

if(Boost_FOUND)
    #set(Boost_LIBRARY_DIR ${BINARYBUILDER_INSTALL_DIR}/lib) # Override this value
    message("Successfully found Boost version ${Boost_VERSION}.")
    message("Boost include dir = ${Boost_INCLUDE_DIR}")
    message("Boost library dir = ${Boost_LIBRARY_DIRS}")
    message("Boost libraries   = ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIRS})
endif(Boost_FOUND)
message("--------------------------------------------------------------")

# TODO: More use of find_package to find things not in BB?

set(REQUIRED     True )
set(NOT_REQUIRED False)


# if (APPLE) # Use the system version
#   find_package(LAPACK)
#   find_package(BLAS)
# else() # Use the conda version
find_external_library(LAPACK ${ISIS_DEPS_DIR} "" "lapack;blas" "" REQUIRED)
set(LAPACK_FOUND 1)
# endif()
    
if(LAPACK_FOUND)
  set(VW_HAVE_PKG_LAPACK 1)
  #if (APPLE)
  #  set(VW_HAVE_PKG_APPLE_LAPACK 1)
  #else()
    set(VW_HAVE_PKG_FLAPACK 1) # --> This specifies that it is a Fortran derived version
                               #      How to check this if the user provides Fortran?
  #endif()
endif()

find_external_library(FLANN   ${ISIS_DEPS_DIR} "" "flann_cpp" REQUIRED)
find_external_library(JPEG    ${ISIS_DEPS_DIR} "" "jpeg"      REQUIRED)
find_external_library(PNG     ${ISIS_DEPS_DIR} "" "png;png16" REQUIRED)
find_external_library(TIFF    ${ISIS_DEPS_DIR} "" "tiff"      REQUIRED)
find_external_library(GEOTIFF ${ISIS_DEPS_DIR} "" "geotiff"   REQUIRED)
find_external_library(Z       ${ISIS_DEPS_DIR} "" "z"         REQUIRED)
find_external_library(PROJ4   ${ISIS_DEPS_DIR} "" "proj"      REQUIRED)
find_external_library(GDAL    ${ISIS_DEPS_DIR} "" "gdal"      REQUIRED)

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
find_external_library(OPENCV ${ISIS_DEPS_DIR} "" "${OPENCV_LIB_NAMES}" REQUIRED)


set(ILMBASE_LIB_NAMES Half Iex Imath IlmThread)
find_external_library(ILMBASE ${ISIS_DEPS_DIR} "" "${ILMBASE_LIB_NAMES}" REQUIRED)
find_external_library(OPENEXR ${ISIS_DEPS_DIR} "OpenEXR" "IlmImf;IlmThread" REQUIRED)

# CMake will try to find these packages if needed,
#  they are not included with BinaryBuilder.
if (VW_ENABLE_PYTHON)
  find_package(PythonLibs)
  if(PYTHONLIBS_FOUND)
    set(VW_HAVE_PKG_PYTHONLIBS 1)
    message("-- Found PythonLibs " ${PYTHON_LIBRARIES})
  else()
    message("Failed to find PythonLibs")
  endif()
endif()

# These packages are only used for ancient VW tools which are not really supported

#find_package(SigC++)
#if(SigC++_FOUND)
#  #message("-- Found SigC++ "  ${SigC++_LIBRARIES})
#  #message("Include dirs: " ${SigC++_INCLUDE_DIRS})
#  include_directories("${SigC++_INCLUDE_DIRS}")
#else()
#  message("Failed to find SIC++")
#endif()

#find_package(FreeType)
#if(FREETYPE_FOUND)
#  message("-- Found FreeType " ${FREETYPE_LIBRARIES})
#  #message("Include dirs: "  ${FREETYPE_INCLUDE_DIRS})
#  include_directories("${FREETYPE_INCLUDE_DIRS}")
#else()
#  message("Failed to find FreeType")
#endif()

#find_package(Cairo)
#if(CAIRO_FOUND)
#  set(VW_HAVE_PKG_CAIRO 1)
#  message("-- Found Cairo " ${CAIRO_LIBRARIES})
#  include_directories("${CAIRO_INCLUDE_DIR}")
#else()
#  message("Failed to find Cairo")
#endif()

#find_package(Cairomm)
#if(CAIROMM_FOUND)
#  set(VW_HAVE_PKG_CAIROMM 1)
#  message("-- Found CairoMM " ${CAIROMM_LIBRARIES})
#  include_directories("${CAIROMM_INCLUDE_DIR}")
#else()
#  message("Failed to find CairoMM")
#endif()


#TODO: Look for each of these installations!
#set(VW_HAVE_PKG_BIGTIFF 1) 
#set(VW_HAVE_PKG_HDF 0)


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

# enable SSE optimizations in some places (development) 
if (VW_ENABLE_SSE)
  set(CXXFLAGS "${CXXFLAGS} -msse4.1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
   message(STATUS, "Enabling SSE")
else() # Compatibility with older systems.
  message(STATUS, "Disabling SSE")
  set(CXXFLAGS "${CXXFLAGS} -mno-sse4.1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-sse4.1")
endif()


# Define to `int' if <sys/types.h> does not define. 
#define VW_ssize_t

# TODO: What is this for?
set(ENABLE_EXCEPTIONS True)

###########################################################################

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

# --- VW_CARTOGRAPHY ------------------------------------------------------------
get_all_source_files( "Cartography"       VW_CARTOGRAPHY_SRC_FILES)
get_all_source_files( "Cartography/tests" VW_CARTOGRAPHY_TEST_FILES)
set(VW_CARTOGRAPHY_LIB_DEPENDENCIES VwFileIO VwCamera ${GDAL_LIBRARIES} 
    ${PROJ4_LIBRARIES} ${GEOTIFF_LIBRARIES}
    ${JPEG_LIBRARIES} ${TIFF_LIBRARIES} ${PNG_LIBRARIES})

# --- VW_CORE ------------------------------------------------------------
get_all_source_files( "Core"       VW_CORE_SRC_FILES)
get_all_source_files( "Core/tests" VW_CORE_TEST_FILES)
set(VW_CORE_LIB_DEPENDENCIES ${Boost_LIBRARIES} pthread)

# --- VW_FILEIO ------------------------------------------------------------
# This is a complicated folder so we need to select files manually
#get_all_source_files( "FileIO" VW_FILEIO_SRC_FILES)
#get_all_source_files( "FileIO/tests" VW_FILEIO_TEST_FILES)
set(VW_FILEIO_LIB_DEPENDENCIES VwMath VwImage ${GDAL_LIBRARIES} ${PROJ4_LIBRARIES}
    ${GEOTIFF_LIBRARIES} ${JPEG_LIBRARIES} ${TIFF_LIBRARIES} ${PNG_LIBRARIES} 
    ${ILMBASE_LIBRARIES} ${OPENEXR_LIBRARIES} ${Z_LIBRARIES})

# --- VW_GEOMETRY ------------------------------------------------------------
get_all_source_files( "Geometry"       VW_GEOMETRY_SRC_FILES)
get_all_source_files( "Geometry/tests" VW_GEOMETRY_TEST_FILES)
set(VW_GEOMETRY_LIB_DEPENDENCIES VwCore VwMath)

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
set(VW_MATH_LIB_DEPENDENCIES VwCore ${FLANN_LIBRARIES} ${LAPACK_LIBRARIES})

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


