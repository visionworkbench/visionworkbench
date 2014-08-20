


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

# Libraries use paths like "vw/blah/blah.h" so we just need to add this one include path
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

# Define boost requirements
set(BOOST_MIN_VERSION "1.55.0")
set(REQUIRED_BOOST_LIBS program_options system filesystem regex date_time thread iostreams)
#set(Boost_DEBUG 1) # Uncomment for useful finding boost debug output
if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(BOOST_ROOT       "${BINARYBUILDER_INSTALL_DIR}/include/boost-1_55")
  set(BOOST_LIBRARYDIR "${BINARYBUILDER_INSTALL_DIR}/lib")
  set(BOOST_INCLUDEDIR "${BINARYBUILDER_INSTALL_DIR}/include/boost-1_55/boost")
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


message("--------------------------------------------------------------")
message("Searching for QT...")
set(REQUIRED_QT_LIBS QtCore QtXml QtGui )
if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set (QT_QMAKE_EXECUTABLE  ${BINARYBUILDER_INSTALL_DIR}/bin/qmake)
  set (QT_MOC_EXECUTABLE    ${BINARYBUILDER_INSTALL_DIR}/bin/moc)
  set (QT_RCC_EXECUTABLE    ${BINARYBUILDER_INSTALL_DIR}/bin/rcc)
  set (QT_UIC_EXECUTABLE    ${BINARYBUILDER_INSTALL_DIR}/bin/uic)
  SET(CMAKE_PREFIX_PATH ${BINARYBUILDER_INSTALL_DIR} ${CMAKE_PREFIX_PATH}) 
endif()

find_package(Qt4 4.8.5 REQUIRED ${REQUIRED_QT_LIBS})

if(Qt4_FOUND)
    include(${QT_USE_FILE}) # Required after using QT find command
    message("QT include dirs  " ${QT_INCLUDES})

    # By default only QtCore and QtGui modules are enabled
    # other modules must be enabled like this:
    SET( QT_USE_QT3SUPPORT TRUE )   
    SET( QT_USE_QTXML      TRUE )

    # We need this to be able to include headers produced by uic in our
    # code (CMAKE_BINARY_DIR holds a path to the build directory, while
    # INCLUDE_DIRECTORIES() works just like INCLUDEPATH from qmake)
    INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )

endif(Qt4_FOUND)
message("--------------------------------------------------------------")

message("--------------------------------------------------------------")
message("Searching for Google Protocol Buffers...")

find_package(Protobuf REQUIRED)
if(PROTOBUF_FOUND)
    set(VW_HAVE_PKG_PROTOBUF 1)
    message("Successfully found Google Protocol Buffers")
    message("---Protobuf library file = ${PROTOBUF_LIBRARIES}")
    message("---Protobuf Include Dir  = ${PROTOBUF_INCLUDE_DIR}")
    include_directories(${PROTOBUF_INCLUDE_DIR})
endif()

#set(PROTOBUF_INCLUDE_DIR ${BINARYBUILDER_INSTALL_DIR}/include   )
#set(PROTOBUF_LIBRARY     ${BINARYBUILDER_INSTALL_DIR}/lib/libprotobuf.so)
##include_directories("${PROTOBUF_INCLUDE_DIR}")
#include_directories(${PROTOBUF_INCLUDE_DIR})

message("--------------------------------------------------------------")


# These packages are included in BinaryBuilder so look for them there.
# - If there is no BinaryBuilder directory the user may need to specify a path.

# Get NAIF includes (Needed for reading/writing SPICE data)
if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(NAIF_FOUND 1)
  set(NAIF_INCLUDE_DIR ${BINARYBUILDER_INSTALL_DIR}/include/naif/)
else() # The user needs to specify this!
  set(NAIF_INCLUDE_DIR "" CACHE FILEPATH "Path to NAIF include files")
endif()
if (NAIF_FOUND)
  message("-- Found NAIF at " ${NAIF_INCLUDE_DIR})
  include_directories("${NAIF_INCLUDE_DIR}")
else()
  message("Failed to find Naif")
endif()


if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(LAPACK_FOUND 1)
  set(LAPACK_LIBRARIES   ${BINARYBUILDER_INSTALL_DIR}/lib/liblapack.so
                         ${BINARYBUILDER_INSTALL_DIR}/lib/libblas.so)
#  set(LAPACK_INCLUDE_DIR ${BINARYBUILDER_INSTALL_DIR}/include/)
else() # The user needs to specify this!
  #set(LAPACK_INCLUDE_DIR "" CACHE FILEPATH "Path to LAPACK include files")
    find_package(LAPACK REQUIRED)
endif()
if(LAPACK_FOUND)
  set(VW_HAVE_PKG_LAPACK 1)
  message("-- Found LAPACK at " ${LAPACK_LIBRARIES})
#  include_directories("${LAPACK_INCLUDE_DIR}")
  
  set(VW_HAVE_PKG_FLAPACK 1) # --> This specifies that it is a Fortran derived version
                             #      How to check this if the user provides Fortran?
else()
  message( FATAL_ERROR "Failed to find REQUIRED library Lapack." )
endif()


if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(FLANN_FOUND 1)
  set(FLANN_LIBRARIES   ${BINARYBUILDER_INSTALL_DIR}/lib/libflann_cpp.so)
  set(FLANN_INCLUDE_DIR ${BINARYBUILDER_INSTALL_DIR}/include/)
else() # The user needs to specify this!
  set(FLANN_INCLUDE_DIR "" CACHE FILEPATH "Path to FLANN include files")
  set(FLANN_LIBRARIES   "" CACHE FILEPATH "Path to FLANN library files")
endif()
if(FLANN_FOUND)
  set(VW_HAVE_PKG_FLANN 1)
  message("-- Found FLANN at " ${FLANN_LIBRARIES})
  include_directories("${FLANN_INCLUDE_DIR}")
else()
  message("Failed to find FLANN")
endif()


if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(JPEG_FOUND 1)
  set(JPEG_LIBRARIES   ${BINARYBUILDER_INSTALL_DIR}/lib/libjpeg.so)
  #set(JPEG_LIBRARY     ${BINARYBUILDER_INSTALL_DIR}/lib/)
  set(JPEG_INCLUDE_DIR ${BINARYBUILDER_INSTALL_DIR}/include/)  
else() # The user needs to specify this!
  #set(JPEG_INCLUDE_DIR "" CACHE FILEPATH "Path to JPEG include files")
  find_package(JPEG REQUIRED)
endif()
if(JPEG_FOUND)
  set(VW_HAVE_PKG_JPEG 1)
  message("-- Found JPEG " ${JPEG_LIBRARIES})
  include_directories("${JPEG_INCLUDE_DIR}")
else()
  message("Failed to find JPEG")
endif()


if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(PNG_FOUND 1)
  set(PNG_LIBRARIES   ${BINARYBUILDER_INSTALL_DIR}/lib/libpng16.so
                      ${BINARYBUILDER_INSTALL_DIR}/lib/libpng.so)
  set(PNG_INCLUDE_DIR ${BINARYBUILDER_INSTALL_DIR}/include/)  
else()
  find_package(PNG REQUIRED)
endif()
if(PNG_FOUND)
  set(VW_HAVE_PKG_PNG 1)
  message("-- Found PNG " ${PNG_LIBRARIES})
  include_directories("${PNG_INCLUDE_DIR}")
else()
  message("Failed to find PNG")
endif()


if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(TIFF_FOUND 1)
  set(TIFF_LIBRARIES   ${BINARYBUILDER_INSTALL_DIR}/lib/libtiff.so)
  set(TIFF_INCLUDE_DIR ${BINARYBUILDER_INSTALL_DIR}/include/)  
else()
  find_package(TIFF REQUIRED)
endif()
if(TIFF_FOUND)
  set(VW_HAVE_PKG_TIFF 1)
  message("-- Found TIFF " ${TIFF_LIBRARIES})
  include_directories("${TIFF_INCLUDE_DIR}")
else()
  message("Failed to find TIFF")
endif()


if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(GDAL_FOUND 1)
  set(GDAL_LIBRARIES   ${BINARYBUILDER_INSTALL_DIR}/lib/libgdal.so)
  set(GDAL_INCLUDE_DIR ${BINARYBUILDER_INSTALL_DIR}/include/)  
else()
  find_package(GDAL REQUIRED)
endif()
if(GDAL_FOUND)
  set(VW_HAVE_PKG_GDAL 1)
  message("-- Found GDAL at " ${GDAL_LIBRARIES})
  include_directories(${GDAL_INCLUDE_DIR})
else()
  message("Failed to find GDAL")
endif()


if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(PROJ4_FOUND 1)  
  set(PROJ4_LIBRARIES   ${BINARYBUILDER_INSTALL_DIR}/lib/libproj.so)
  set(PROJ4_INCLUDE_DIR ${BINARYBUILDER_INSTALL_DIR}/include)
else()
  set(PROJ4_INCLUDE_DIR "" CACHE FILEPATH "Path to OpenEXR include files")
  set(PROJ4_LIBRARIES   "" CACHE FILEPATH "Path to OpenEXR libraries")  
endif()
if(PROJ4_FOUND)
  set(VW_HAVE_PKG_PROJ4 1)
  message("-- Found Proj4 at " ${PROJ4_LIBRARIES})
  include_directories("${PROJ4_INCLUDE_DIR}")
else()
  message("Failed to find Proj4")
endif()  


if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(ILMBASE_FOUND 1)
  set(ILMBASE_LIBRARIES   ${BINARYBUILDER_INSTALL_DIR}/lib/libHalf.so
                          ${BINARYBUILDER_INSTALL_DIR}/lib/libIex.so
                          ${BINARYBUILDER_INSTALL_DIR}/lib/libImath.so
                          ${BINARYBUILDER_INSTALL_DIR}/lib/libIlmThread.so 
                          ${BINARYBUILDER_INSTALL_DIR}/lib/libIexMath.so)
  set(ILMBASE_INCLUDE_DIR ${BINARYBUILDER_INSTALL_DIR}/include) 
else()
  set(ILMBASE_INCLUDE_DIR "" CACHE FILEPATH "Path to ILMBase include files")
  set(ILMBASE_LIB_DIR     "" CACHE FILEPATH "Path to ILMBase library files")
  set(ILMBASE_LIBRARIES   ${ILMBASE_LIB_DIR}/lib/libHalf.so
                          ${ILMBASE_LIB_DIR}/lib/libIex.so
                          ${ILMBASE_LIB_DIR}/lib/libImath.so
                          ${ILMBASE_LIB_DIR}/lib/libIlmThread.so 
                          ${ILMBASE_LIB_DIR}/lib/libIexMath.so)                         
endif()
if(ILMBASE_FOUND)
  set(VW_HAVE_PKG_ILMBASE 1)
  message("-- Found ILMBase at " ${ILMBASE_LIBRARIES})
  include_directories("${ILMBASE_INCLUDE_DIR}")
else()
  message("Failed to find ILMBase")
endif()


if(BINARYBUILDER_INSTALL_DIR) # Look in BinaryBuilder directory
  set(OPENEXR_FOUND 1)
  set(OPENEXR_LIBRARIES   ${BINARYBUILDER_INSTALL_DIR}/lib/libIlmImf.so
                          ${BINARYBUILDER_INSTALL_DIR}/lib/libIlmThread.so )
  set(OPENEXR_INCLUDE_DIR ${BINARYBUILDER_INSTALL_DIR}/include/OpenEXR) 
else()
  set(OPENEXR_INCLUDE_DIR "" CACHE FILEPATH "Path to OpenEXR include files")
  set(OPENEXR_LIB_DIR     "" CACHE FILEPATH "Path to OpenEXR library files")
  set(OPENEXR_LIBRARIES   ${OPENEXR_LIB_DIR}/libIlmImf.so
                          ${OPENEXR_LIB_DIR}/libIlmThread.so )
endif()
if(OPENEXR_FOUND)
  set(VW_HAVE_PKG_OPENEXR 1)
  message("Found OpenEXR at " ${OPENEXR_LIBRARIES})
  include_directories("${OPENEXR_INCLUDE_DIR}")
else()
  message("Failed to find OpenEXR")
endif()  
  



# CMake will try to find these packages,
#  they are not included with BinaryBuilder.


find_package(PythonLibs REQUIRED)
if(PYTHONLIBS_FOUND)
  set(VW_HAVE_PKG_PYTHONLIBS 1)
  message("-- Found PYTHONLIBS " ${PYTHON_LIBRARIES})
else()
  message("Failed to find PythonLibs")
endif()


find_package(OpenCV)
if(OPENCV_FOUND)
  set(VW_HAVE_PKG_OPENCV 1)
  #message("Found OPENCV " ${OPENCV_LIBRARIES})
  include_directories("${OPENCV_INCLUDE_DIR}")
else()
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
#set(VW_HAVE_PKG_BIGTIFF 1) TODO: Do we need this one?
#set(VW_HAVE_PKG_HDF 0)
#set(VW_HAVE_PKG_RABBITMQ_C 0) # Only used for plate!
#set(VW_HAVE_PKG_ZEROMQ 1)
#set(VW_HAVE_PKG_LIBKML 0) # Only used for plate!



# TODO: What is this for?
set(ENABLE_EXCEPTIONS True)


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
set(VW_BUNDLEADJUSTMENT_LIB_DEPENDENCIES VwCore VwMosaic VwCartography)

# --- VW_CAMERA ------------------------------------------------------------
get_all_source_files( "Camera"       VW_CAMERA_SRC_FILES)
get_all_source_files( "Camera/tests" VW_CAMERA_TEST_FILES)
set(VW_CAMERA_LIB_DEPENDENCIES VwMath VwImage)

# --- VW_CARTOGRAPHY ------------------------------------------------------------
get_all_source_files( "Cartography"       VW_CARTOGRAPHY_SRC_FILES)
get_all_source_files( "Cartography/tests" VW_CARTOGRAPHY_TEST_FILES)
set(VW_CARTOGRAPHY_LIB_DEPENDENCIES VwFileIO VwCamera)

# --- VW_CORE ------------------------------------------------------------
get_all_source_files( "Core"       VW_CORE_SRC_FILES)
get_all_source_files( "Core/tests" VW_CORE_TEST_FILES)
set(VW_CORE_LIB_DEPENDENCIES )

# --- VW_FILEIO ------------------------------------------------------------
# This is a complicated folder so we need to select files manually
#get_all_source_files( "FileIO" VW_FILEIO_SRC_FILES)
#get_all_source_files( "FileIO/tests" VW_FILEIO_TEST_FILES)
set(VW_FILEIO_LIB_DEPENDENCIES VwMath VwImage)

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
set(VW_HDR_LIB_DEPENDENCIES VwFileIO VwCamera)

# --- VW_IMAGE ------------------------------------------------------------
# This is a complicated folder so we need to select files manually
#get_all_source_files( "Image" VW_IMAGE_SRC_FILES)
#get_all_source_files( "Image/tests" VW_IMAGE_TEST_FILES)
set(VW_IMAGE_LIB_DEPENDENCIES VwMath)

# --- VW_INTERESTPOINT ------------------------------------------------------------
get_all_source_files( "InterestPoint"       VW_INTERESTPOINT_SRC_FILES)
get_all_source_files( "InterestPoint/tests" VW_INTERESTPOINT_TEST_FILES)
set(VW_INTERESTPOINT_LIB_DEPENDENCIES VwCamera VwFileIO)

# --- VW_MATH ------------------------------------------------------------
get_all_source_files( "Math"       VW_MATH_SRC_FILES)
get_all_source_files( "Math/tests" VW_MATH_TEST_FILES)
set(VW_IMAGE_LIB_DEPENDENCIES VwCore)

# --- VW_MOSAIC ------------------------------------------------------------
get_all_source_files( "Mosaic"       VW_MOSAIC_SRC_FILES)
get_all_source_files( "Mosaic/tests" VW_MOSAIC_TEST_FILES)
set(VW_MOSAIC_LIB_DEPENDENCIES VwFileIO VwCamera)

# DEPRECATED
## --- VW_PLATE ------------------------------------------------------------
#get_all_source_files( "Plate"       VW_PLATE_SRC_FILES)
#get_all_source_files( "Plate/tests" VW_PLATE_TEST_FILES)
#set(VW_PLATE_LIB_DEPENDENCIES VwMath)

# PYTHON directory????

# --- VW_STEREO ------------------------------------------------------------
get_all_source_files( "Stereo"       VW_STEREO_SRC_FILES)
get_all_source_files( "Stereo/tests" VW_STEREO_TEST_FILES)
set(VW_STEREO_LIB_DEPENDENCIES VwMath VwImage VwCamera VwGeometry)

#=========================================================================


# TODO: Don't just use this for all tools
set(FULL_LIBRARY_LIST
      ${LAPACK_LIBRARIES} 
      ${JPEG_LIBRARIES} 
      ${GDAL_LIBRARIES}
      ${TIFF_LIBRARIES}
      ${PNG_LIBRARIES}
      ${PROJ4_LIBRARIES}
      ${OPENEXR_LIBRARIES}
      ${OPENCV_LIBRARIES}
      ${PROTOBUF_LIBRARIES}
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







