


# Add path for custom modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}" )
message("Module path = ${CMAKE_MODULE_PATH}")


#===============================================================================
# Configuration options

# --- Directory paths ---
set(BINARYBUILDER_INSTALL_DIR "" CACHE FILEPATH "Path to Binary Builder installation directory")
if(BINARYBUILDER_INSTALL_DIR)
  message("Looking for dependencies in BinaryBuilder directory: ${BINARYBUILDER_INSTALL_DIR}") 
  # Add BINARYBUILDER folders to CMAKE paths
  set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "${BINARYBUILDER_INSTALL_DIR}")
  set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "${BINARYBUILDER_INSTALL_DIR}")
endif()

# Libraries use paths like "vw/blah/blah.h" so we need to include the top of the tree
include_directories( ../ )
include_directories( ${CMAKE_BINARY_DIR}/src/ ) # Generated protocol buffers go here!


get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

# --- Custom options ---
option(BUILD_SHARED_LIBS "Produce shared libraries." TRUE)

# --- Fixed options ---
set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON )
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#==================================================================================
# Find dependencies

# TODO: Create a cmake folder to contain the various find files.

message("--------------------------------------------------------------")
message("Searching for boost...")

# Define boost requirements
set(BOOST_MIN_VERSION "1.59.0")
set(REQUIRED_BOOST_LIBS program_options system filesystem regex date_time thread iostreams)
#set(Boost_DEBUG 1) # Uncomment for useful finding boost debug output
if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(BOOST_ROOT       "${BINARYBUILDER_INSTALL_DIR}/include/boost-1_59")
  set(BOOST_LIBRARYDIR "${BINARYBUILDER_INSTALL_DIR}/lib")
  set(BOOST_INCLUDEDIR "${BINARYBUILDER_INSTALL_DIR}/include/boost-1_59/boost")
  set(Boost_NO_SYSTEM_PATHS ON) #This does not seem to work
endif()  
# Now search for Boost using the available information
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS ${REQUIRED_BOOST_LIBS})

if(Boost_FOUND)
    message("Successfully found Boost version ${Boost_VERSION}!")
    message("---Boost Include Dir = ${Boost_INCLUDE_DIR}")
    #set(Boost_LIBRARY_DIR ${BINARYBUILDER_INSTALL_DIR}/lib) # Override this value
    message("---Boost Library Dir = ${Boost_LIBRARY_DIRS}")
    message("---Boost libraries   = ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    #include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIR})
    #link_directories(${Boost_LIBRARY_DIRS})
endif(Boost_FOUND)
message("--------------------------------------------------------------")



# TODO: More use of find_package to find things not in BB?

set(REQUIRED     True )
set(NOT_REQUIRED False)

# TODO: Only needed for ASP!!!
## Get NAIF includes (Needed for reading/writing SPICE data)
#find_external_library(NAIF "naif" "" REQUIRED)

find_external_library(LAPACK "" "liblapack.so;libblas.so" REQUIRED)
if(LAPACK_FOUND)
  set(VW_HAVE_PKG_FLAPACK 1) # --> This specifies that it is a Fortran derived version
                             #      How to check this if the user provides Fortran?
endif()


find_external_library(FLANN "" "libflann_cpp.so" NOT_REQUIRED)

find_external_library(JPEG "" "libjpeg.so" NOT_REQUIRED)
find_external_library(PNG "" "libpng.so;libpng16.so" NOT_REQUIRED)
find_external_library(TIFF "" "libtiff.so" NOT_REQUIRED)

find_external_library(Z "" "libz.so" NOT_REQUIRED)
find_external_library(GDAL "" "libgdal.so" NOT_REQUIRED)
find_external_library(PROJ4 "" "libproj.so" NOT_REQUIRED)

set(OPENCV_LIB_NAMES libopencv_calib3d.so     libopencv_reg.so
                     libopencv_core.so        libopencv_shape.so
                     libopencv_features2d.so  libopencv_stitching.so
                     libopencv_flann.so       libopencv_superres.so
                     libopencv_hdf.so         libopencv_surface_matching.so
                     libopencv_highgui.so     libopencv_videoio.so
                     libopencv_imgcodecs.so   libopencv_video.so
                     libopencv_imgproc.so     libopencv_xfeatures2d.so
                     libopencv_ml.so          libopencv_ximgproc.so
                     libopencv_objdetect.so   libopencv_xobjdetect.so
                     libopencv_photo.so       libopencv_xphoto.so)
find_external_library(OPENCV "" "${OPENCV_LIB_NAMES}" NOT_REQUIRED)


set(ILMBASE_LIB_NAMES libHalf.so libIex.so libImath.so libIlmThread.so libIexMath.so)
find_external_library(ILMBASE "" "${ILMBASE_LIB_NAMES}" NOT_REQUIRED)

find_external_library(OPENEXR "OpenEXR" "libIlmImf.so;libIlmThread.so" NOT_REQUIRED)


# CMake will try to find these packages,
#  they are not included with BinaryBuilder.

find_package(PythonLibs)
if(PYTHONLIBS_FOUND)
  set(VW_HAVE_PKG_PYTHONLIBS 1)
  message("-- Found PYTHONLIBS " ${PYTHON_LIBRARIES})
else()
  message("Failed to find PythonLibs")
endif()


find_package(SigC++)
if(SigC++_FOUND)
  #message("-- Found SigC++ "  ${SigC++_LIBRARIES})
  #message("Include dirs: " ${SigC++_INCLUDE_DIRS})
  include_directories("${SigC++_INCLUDE_DIRS}")
else()
  message("Failed to find SIC++")
endif()

find_package(FreeType)
if(FREETYPE_FOUND)
  message("-- Found FreeType " ${FREETYPE_LIBRARIES})
  #message("Include dirs: "  ${FREETYPE_INCLUDE_DIRS})
  include_directories("${FREETYPE_INCLUDE_DIRS}")
else()
  message("Failed to find FreeType")
endif()

find_package(Cairo)
if(CAIRO_FOUND)
  set(VW_HAVE_PKG_CAIRO 1)
  message("-- Found Cairo " ${CAIRO_LIBRARIES})
  include_directories("${CAIRO_INCLUDE_DIR}")
else()
  message("Failed to find Cairo")
endif()

find_package(Cairomm)
if(CAIROMM_FOUND)
  set(VW_HAVE_PKG_CAIROMM 1)
  message("-- Found CairoMM " ${CAIROMM_LIBRARIES})
  include_directories("${CAIROMM_INCLUDE_DIR}")
else()
  message("Failed to find CairoMM")
endif()


#TODO: Look for each of these installations!
#set(VW_HAVE_PKG_BIGTIFF 1) 
#set(VW_HAVE_PKG_HDF 0)


###########################################################################
# These are some options that the user might like to set manually

set(VW_DEFAULT_CACHE_SIZE_MB 768)

# set the default number of processing threads for multi-threaded operations
set(VW_NUM_THREADS 4)

# enable image bounds checking (SLOW!) 
set(VW_ENABLE_BOUNDS_CHECK 0)

# ~/.vwrc support
set(VW_ENABLE_CONFIG_FILE 1) 

# enable the C++ exception mechanism 
set(VW_ENABLE_EXCEPTIONS 1) 

# enable SSE optimizations in some places (development) 
set(VW_ENABLE_SSE 1) 
if (VW_ENABLE_SSE)
  set(CXXFLAGS "${CXXFLAGS} -msse4.1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
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
set(VW_BUNDLEADJUSTMENT_LIB_DEPENDENCIES VwMosaic VwFileIO VwInterestPoint VwStereo)

# --- VW_CAMERA ------------------------------------------------------------
get_all_source_files( "Camera"       VW_CAMERA_SRC_FILES)
get_all_source_files( "Camera/tests" VW_CAMERA_TEST_FILES)
set(VW_CAMERA_LIB_DEPENDENCIES VwMath VwImage)

# --- VW_CARTOGRAPHY ------------------------------------------------------------
get_all_source_files( "Cartography"       VW_CARTOGRAPHY_SRC_FILES)
get_all_source_files( "Cartography/tests" VW_CARTOGRAPHY_TEST_FILES)
set(VW_CARTOGRAPHY_LIB_DEPENDENCIES VwFileIO VwCamera ${GDAL_LIBRARIES} ${PROJ4_LIBRARIES})

# --- VW_CORE ------------------------------------------------------------
get_all_source_files( "Core"       VW_CORE_SRC_FILES)
get_all_source_files( "Core/tests" VW_CORE_TEST_FILES)
set(VW_CORE_LIB_DEPENDENCIES ${Boost_LIBRARIES})

# --- VW_FILEIO ------------------------------------------------------------
# This is a complicated folder so we need to select files manually
#get_all_source_files( "FileIO" VW_FILEIO_SRC_FILES)
#get_all_source_files( "FileIO/tests" VW_FILEIO_TEST_FILES)
set(VW_FILEIO_LIB_DEPENDENCIES VwMath VwImage
    ${JPEG_LIBRARIES} ${TIFF_LIBRARIES} ${PNG_LIBRARIES} ${OPENEXR_LIBRARIES} ${Z_LIBRARIES})

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


# PYTHON directory????

# --- VW_STEREO ------------------------------------------------------------
get_all_source_files( "Stereo"       VW_STEREO_SRC_FILES)
get_all_source_files( "Stereo/tests" VW_STEREO_TEST_FILES)
set(VW_STEREO_LIB_DEPENDENCIES VwMath VwImage VwCamera VwGeometry)

#=========================================================================
# Add all of the library subdirectories

add_subdirectory(BundleAdjustment)
add_subdirectory(Camera)
add_subdirectory(Cartography)
add_subdirectory(Core)
add_subdirectory(FileIO)
add_subdirectory(Geometry)
#add_subdirectory(GPU) # DEPRECATED
#add_subdirectory(gui) # DEPRECATED
add_subdirectory(HDR)
add_subdirectory(Image)
add_subdirectory(InterestPoint)
add_subdirectory(Math)
add_subdirectory(Mosaic)
#add_subdirectory(Plate) # DEPRECATED
add_subdirectory(Stereo)
add_subdirectory(tools)



#=========================================================================
# Targets are defined in /tools/CMakeLists.txt




