


# Add path for custom modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}" )
message("Module path = ${CMAKE_MODULE_PATH}")


#===============================================================================
# Configuration options

# --- Directory paths ---
#set(BASESYSTEM_INSTALL_DIR "" CACHE FILEPATH "Path to AutoTools installation directory")
set(BASESYSTEM_INSTALL_DIR "/home/smcmich1/repo/BaseSystem")
if(NOT BASESYSTEM_INSTALL_DIR)
  message(SEND_ERROR "You need to set BASESYSTEM_INSTALL_DIR: ${BASESYSTEM_INSTALL_DIR}")
endif()

message("Base System path = ${BASESYSTEM_INSTALL_DIR}")


# Add BaseSystem folders to CMAKE paths
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "${BASESYSTEM_INSTALL_DIR}")
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "${BASESYSTEM_INSTALL_DIR}")



# Libraries use paths like "asp/blah/blah.h" so we just need to add this one include path
include_directories( ../ )


# --- Custom options ---
option(BUILD_SHARED_LIBS "Produce shared libraries." TRUE)

# --- Fixed options ---
set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON )
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#==================================================================================
# Find dependencies

message("--------------------------------------------------------------")
message("Searching for boost...")

# TODO: Set this to make sure the version matches Binary Builder
# - Alternately, ONLY accept the boost version in Binary Builder
# BINARY_BUILDER_BOOST_VERSION

set(BOOST_ROOT "${BASESYSTEM_INSTALL_DIR}/include/boost-1_55")
set(BOOST_LIBRARYDIR  "${BASESYSTEM_INSTALL_DIR}/lib")
set(BOOST_INCLUDEDIR "${BASESYSTEM_INSTALL_DIR}/include/boost-1_55/boost")
set(Boost_NO_SYSTEM_PATHS ON) #This does not seem to work
#set(Boost_DEBUG 1) # Uncomment for useful finding boost debug output
set(BOOST_MIN_VERSION "1.55.0")
set(REQUIRED_BOOST_LIBS program_options system filesystem regex date_time thread iostreams)
#TODO Set additional boost libraries

find_package(Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS ${REQUIRED_BOOST_LIBS})

if(Boost_FOUND)
    message("Successfully found Boost version ${Boost_VERSION}!")
    message("---Boost Include Dir = ${Boost_INCLUDE_DIR}")
    #set(Boost_LIBRARY_DIR ${BASESYSTEM_INSTALL_DIR}/lib) # Override this value
    message("---Boost Library Dir = ${Boost_LIBRARY_DIRS}")
    message("---Boost libraries = ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    #include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIR})
    #link_directories(${Boost_LIBRARY_DIRS})
endif(Boost_FOUND)
message("--------------------------------------------------------------")


message("--------------------------------------------------------------")
message("Searching for QT...")
# TODO: Make sure this finds the correct install!
set (QT_QMAKE_EXECUTABLE  ${BASESYSTEM_INSTALL_DIR}/bin/qmake)
set (QT_MOC_EXECUTABLE    ${BASESYSTEM_INSTALL_DIR}/bin/moc)
set (QT_RCC_EXECUTABLE    ${BASESYSTEM_INSTALL_DIR}/bin/rcc)
set (QT_UIC_EXECUTABLE    ${BASESYSTEM_INSTALL_DIR}/bin/uic)
#SET(CMAKE_PREFIX_PATH ${BASESYSTEM_INSTALL_DIR} ${CMAKE_PREFIX_PATH}) 
set(REQUIRED_QT_LIBS QtCore QtXml QtGui ) #TODO: Figure this out
find_package(Qt4 4.8.5 REQUIRED ${REQUIRED_QT_LIBS})


if(Qt4_FOUND)
    include(${QT_USE_FILE}) # Required after using QT find command
    message("QT include dirs  " ${QT_INCLUDES})

    # By default only QtCore and QtGui modules are enabled
    # other modules must be enabled like this:
    SET( QT_USE_QT3SUPPORT TRUE )   
    SET( QT_USE_QTXML TRUE )

    # We need this to be able to include headers produced by uic in our
    # code (CMAKE_BINARY_DIR holds a path to the build directory, while
    # INCLUDE_DIRECTORIES() works just like INCLUDEPATH from qmake)
    INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )

endif(Qt4_FOUND)
message("--------------------------------------------------------------")

#message("--------------------------------------------------------------")
#message("Searching for Google Protocol Buffers...")

#set(BOOST_ROOT "${BASESYSTEM_INSTALL_DIR}/boost")
#set(BOOST_MIN_VERSION "1.54.0")
#FIND_PACKAGE(Protobuf ${BOOST_MIN_VERSION} REQUIRED)

#if(PROTOBUF_FOUND)
#    message("Successfully found Google Protocol Buffers")
#    message("---Protobuf library file = ${PROTOBUF_LIBRARY}")
#    message("---Protobuf Include Dir = ${PROTOBUF_INCLUDE_DIR}")
#    include_directories(${PROTOBUF_INCLUDE_DIRS})
#    target_link_libraries(foo ${PROTOBUF_LIBRARIES})
#endif()
#message("--------------------------------------------------------------")
set(PROTOBUF_FOUND 1)
set(VW_HAVE_PKG_PROTOBUF 1)
set(PROTOBUF_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include   )
set(PROTOBUF_LIBRARY     ${BASESYSTEM_INSTALL_DIR}/lib/libprotobuf.so)
#include_directories("${PROTOBUF_INCLUDE_DIR}")
include_directories(${PROTOBUF_INCLUDE_DIR})

# Many of these are currently found in BaseSystem.
# - We should look for them elsewhere if BaseSystem does not exist.


# Get NAIF includes (Needed for reading/writing SPICE data)
set(NAIF_FOUND 1)
if (NAIF_FOUND)
  set(NAIF_INCLUDE_DIRS ${BASESYSTEM_INSTALL_DIR}/include/naif/)
  include_directories("${NAIF_INCLUDE_DIRS}")
  message("-- Found NAIF at " ${NAIF_INCLUDE_DIRS})
else()
  message("Failed to find NAIF")
endif()

set(LAPACK_FOUND 1)
if(LAPACK_FOUND)
  set(VW_HAVE_PKG_LAPACK 1)
  set(LAPACK_LIBRARIES   ${BASESYSTEM_INSTALL_DIR}/lib/liblapack.so
                       ${BASESYSTEM_INSTALL_DIR}/lib/libblas.so)
  set(LAPACK_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include/)
  message("-- Found LAPACK at " ${LAPACK_LIBRARIES})
else()
  message("Failed to find LAPACK")
endif()

set(FLANN_FOUND 1)
if(FLANN_FOUND)
  set(VW_HAVE_PKG_FLANN 1)
  set(FLANN_LIBRARIES   ${BASESYSTEM_INSTALL_DIR}/lib/libflann_cpp.so)
  set(FLANN_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include/)
  message("-- Found FLANN at " ${FLANN_LIBRARIES})
else()
  message("Failed to find FLANN")
endif()

set(JPEG_FOUND 1)
if(JPEG_FOUND)
  set(VW_HAVE_PKG_JPEG 1)
  set(JPEG_LIBRARIES   ${BASESYSTEM_INSTALL_DIR}/lib/libjpeg.so)
  set(JPEG_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include/)  
  message("-- Found JPEG " ${JPEG_LIBRARY})
else()
  message("Failed to find JPEG")
endif()

find_package(PNG REQUIRED)
if(PNG_FOUND)
  set(VW_HAVE_PKG_PNG 1)
  #message("Found PNG " ${PNG_LIBRARY})
else()
  message("Failed to find PNG")
endif()
#set(PNG_FOUND) # PNG may be available in BaseSystem too!



find_package(TIFF REQUIRED)
if(TIFF_FOUND)
  set(VW_HAVE_PKG_TIFF 1)
  #message("Found TIFF " ${TIFF_LIBRARY})
else()
  message("Failed to find TIFF")
endif()

find_package(GDAL REQUIRED)
if(GDAL_FOUND)
  set(VW_HAVE_PKG_GDAL 1)
  #message("Found GDAL at " ${GDAL_LIBRARY})
  include_directories(${GDAL_INCLUDE_DIR})
else()
  message("Failed to find GDAL")
endif()

find_package(Protobuf REQUIRED)
if(PROTOBUF_FOUND)
  set(VW_HAVE_PKG_PROTOBUF 1)
  #message("Found PROTOBUF " ${PROTOBUF_LIBRARY})
else()
  message("Failed to find PROTOBUF")
endif()

find_package(PythonLibs REQUIRED)
if(PYTHONLIBS_FOUND)
  set(VW_HAVE_PKG_PYTHONLIBS 1)
  #message("Found PYTHONLIBS " ${PYTHON_LIBRARIES})
else()
  message("Failed to find PYTHONLIBS")
endif()


#TODO: Use FindProj cmake file to do this!
set(VW_HAVE_PKG_PROJ4 1)
set(PROJ4_LIBRARY     ${BASESYSTEM_INSTALL_DIR}/lib/libproj.so)
set(PROJ4_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include)


set(VW_HAVE_PKG_OPENEXR 1)
set(OPENEXR_LIBRARIES   ${BASESYSTEM_INSTALL_DIR}/lib/libIlmImf.so
                        ${BASESYSTEM_INSTALL_DIR}/lib/libIlmThread.so )
set(OPENEXR_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include/OpenEXR)
include_directories("${OPENEXR_INCLUDE_DIR}")

# TODO: Where is this in BinaryBuilder?
set(VW_HAVE_PKG_ILMBASE 1)
set(VW_HAVE_PKG_FLAPACK 1)


find_package(OpenCV)
if(OPENCV_FOUND)
  set(VW_HAVE_PKG_OPENCV 1)
  #message("Found OPENCV " ${OPENCV_LIBRARIES})
  include_directories("${OPENCV_INCLUDE_DIR}")
else()
  message("Failed to find OPENCV")
endif()

find_package(SigC++)
if(SigC++_FOUND)
  #message("Found SigC++ " ${SigC++_LIBRARIES})
  #message("Include dirs: " ${SigC++_INCLUDE_DIRS})
  include_directories("${SigC++_INCLUDE_DIRS}")
else()
  message("Failed to find SIC++")
endif()

find_package(FreeType)
if(FREETYPE_FOUND)
  #message("Found FreeType " ${FREETYPE_LIBRARIES})
  #message("Include dirs: " ${FREETYPE_INCLUDE_DIRS})
  include_directories("${FREETYPE_INCLUDE_DIRS}")
else()
  message("Failed to find FREETYPE")
endif()

find_package(Cairo)
if(CAIRO_FOUND)
  set(VW_HAVE_PKG_CAIRO 1)
  #message("Found Cairo " ${CAIRO_LIBRARIES})
  include_directories("${CAIRO_INCLUDE_DIR}")
else()
  message("Failed to find CAIRO")
endif()

find_package(Cairomm)
if(CAIROMM_FOUND)
  set(VW_HAVE_PKG_CAIROMM 1)
  #message("Found CairoMM " ${CAIROMM_LIBRARIES})
  include_directories("${CAIROMM_INCLUDE_DIR}")
else()
  message("Failed to find CAIROMM")
endif()




#TODO: Look for each of these installations!
#set(VW_HAVE_PKG_BIGTIFF 1) TODO: Do we need this one?
set(VW_HAVE_PKG_HDF 0)
#set(VW_HAVE_PKG_RABBITMQ_C 0) # Only used for plate!
set(VW_HAVE_PKG_ZEROMQ 1)
#set(VW_HAVE_PKG_LIBKML 0) # Only used for plate!
set(VW_HAVE_PKG_ILMBASE 1)



# TODO: What is this for?
set(ENABLE_EXCEPTIONS True)


message("--------------------------------------------------------------")
#message("Listing all include directories...")
#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()
#message("--------------------------------------------------------------")


# Now that we have found all our external dependencies, generate a config.h file
include("GenerateConfig.cmake")
#GenerateConfig()


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
set(VW_BUNDLEADJUSTMENT_LIB_DEPENDENCIES VwCore VwMath)

# --- VW_CAMERA ------------------------------------------------------------
get_all_source_files( "Camera"       VW_CAMERA_SRC_FILES)
get_all_source_files( "Camera/tests" VW_CAMERA_TEST_FILES)
set(VW_CAMERA_LIB_DEPENDENCIES VwCore VwMath)

# --- VW_CARTOGRAPHY ------------------------------------------------------------
get_all_source_files( "Cartography"       VW_CARTOGRAPHY_SRC_FILES)
get_all_source_files( "Cartography/tests" VW_CARTOGRAPHY_TEST_FILES)
set(VW_CARTOGRAPHY_LIB_DEPENDENCIES VwCore VwMath)

# --- VW_CORE ------------------------------------------------------------
get_all_source_files( "Core"       VW_CORE_SRC_FILES)
get_all_source_files( "Core/tests" VW_CORE_TEST_FILES)
set(VW_CORE_LIB_DEPENDENCIES )

# --- VW_FILEIO ------------------------------------------------------------
# This is a complicated folder so we need to select files manually
#get_all_source_files( "FileIO" VW_FILEIO_SRC_FILES)
#get_all_source_files( "FileIO/tests" VW_FILEIO_TEST_FILES)
set(VW_FILEIO_LIB_DEPENDENCIES VwCore)

# --- VW_GEOMETRY ------------------------------------------------------------
get_all_source_files( "Geometry"       VW_GEOMETRY_SRC_FILES)
get_all_source_files( "Geometry/tests" VW_GEOMETRY_TEST_FILES)
set(VW_GEOMETRY_LIB_DEPENDENCIES VwCore VwMath)

# DEPRECATED
## --- VW_GPU ------------------------------------------------------------
#get_all_source_files( "GPU"       VW_GPU_SRC_FILES)
#get_all_source_files( "GPU/tests" VW_GPU_TEST_FILES)
#set(VW_GPU_LIB_DEPENDENCIES VwCore VwMath)

# DEPRECATED
## --- VW_GUI ------------------------------------------------------------
#get_all_source_files( "gui"       VW_GUI_SRC_FILES)
#get_all_source_files( "gui/tests" VW_GUI_TEST_FILES)
#set(VW_GUI_LIB_DEPENDENCIES VwCore VwMath)

# --- VW_HDR ------------------------------------------------------------
get_all_source_files( "HDR"       VW_HDR_SRC_FILES)
get_all_source_files( "HDR/tests" VW_HDR_TEST_FILES)
set(VW_HDR_LIB_DEPENDENCIES VwCore VwMath)

# --- VW_IMAGE ------------------------------------------------------------
# This is a complicated folder so we need to select files manually
#get_all_source_files( "Image" VW_IMAGE_SRC_FILES)
#get_all_source_files( "Image/tests" VW_IMAGE_TEST_FILES)
set(VW_IMAGE_LIB_DEPENDENCIES VwCore)

# --- VW_INTERESTPOINT ------------------------------------------------------------
get_all_source_files( "InterestPoint"       VW_INTERESTPOINT_SRC_FILES)
get_all_source_files( "InterestPoint/tests" VW_INTERESTPOINT_TEST_FILES)
set(VW_INTERESTPOINT_LIB_DEPENDENCIES VwCore VwMath)

# --- VW_MATH ------------------------------------------------------------
get_all_source_files( "Math"       VW_MATH_SRC_FILES)
get_all_source_files( "Math/tests" VW_MATH_TEST_FILES)
set(VW_IMAGE_LIB_DEPENDENCIES VwMath)

# --- VW_MOSAIC ------------------------------------------------------------
get_all_source_files( "Mosaic"       VW_MOSAIC_SRC_FILES)
get_all_source_files( "Mosaic/tests" VW_MOSAIC_TEST_FILES)
set(VW_MOSAIC_LIB_DEPENDENCIES VwMath)

# DEPRECATED
## --- VW_PLATE ------------------------------------------------------------
#get_all_source_files( "Plate"       VW_PLATE_SRC_FILES)
#get_all_source_files( "Plate/tests" VW_PLATE_TEST_FILES)
#set(VW_PLATE_LIB_DEPENDENCIES VwMath)

# PYTHON directory????

# --- VW_STEREO ------------------------------------------------------------
get_all_source_files( "Stereo"       VW_STEREO_SRC_FILES)
get_all_source_files( "Stereo/tests" VW_STEREO_TEST_FILES)
set(VW_STEREO_LIB_DEPENDENCIES VwMath)

#=========================================================================


# TODO: Don't just use this for all tools
set(FULL_LIBRARY_LIST
	${LAPACK_LIBRARIES} 
	${JPEG_LIBRARIES} 
	${GDAL_LIBRARY}
	${PROJ4_LIBRARY}
	${OPENEXR_LIBRARIES}
	${OPENCV_LIBRARIES}
	${PROTOBUF_LIBRARY}
  ${SigC++_LIBRARIES}
	${CAIRO_LIBRARIES}
	${CAIROMM_LIBRARIES}
	VwBundleAdjustment
    VwCamera
    VwCartography
    VwCore   
    VwFileIO
    VwGeometry
    #VwGpu
    #VwGui
    VwHdr
    VwImage
    VwInterestPoint
    VwMath
    VwMosaic
    #VwPlate
    VwStereo
	${Boost_LIBRARIES} 
)



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







