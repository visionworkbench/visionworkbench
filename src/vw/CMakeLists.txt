


# Add path for custom modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}" )



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

#find_package(Boost 1.50.0
#  COMPONENTS
#  program_options
#  system
#  filesystem
#  thread
#  )

message("--------------------------------------------------------------")
message("Searching for boost...")

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

# Get NAIF includes (Needed for reading/writing SPICE data)
set(NAIF_FOUND True)
set(NAIF_INCLUDE_DIRS ${BASESYSTEM_INSTALL_DIR}/include/naif/)
include_directories("${NAIF_INCLUDE_DIRS}")

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
set(PROTOBUF_FOUND)
set(PROTOBUF_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include   )
set(PROTOBUF_LIBRARY     ${BASESYSTEM_INSTALL_DIR}/lib/libprotobuf.so)
#include_directories("${PROTOBUF_INCLUDE_DIR}")
include_directories(${PROTOBUF_INCLUDE_DIR})





# Look for more required packages
#find_package(LAPACK REQUIRED)  This should be in BaseSystem
set(HAVE_PKG_LAPACK LAPACK_FOUND)
if(LAPACK_FOUND)
  message("Found LAPACK")
endif()

#TODO: Find in correct place!  Do we want to use find or force to BaseSystem?
#find_package(JPEG REQUIRED)
set(JPEG_FOUND True)
set(JPEG_LIBRARIES   ${BASESYSTEM_INSTALL_DIR}/lib/libjpeg.so)
set(JPEG_LIBRARY     ${BASESYSTEM_INSTALL_DIR}/lib/)
set(JPEG_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include/)

set(HAVE_PKG_JPEG JPEG_FOUND)
if(JPEG_FOUND)
  message("Found JPEG " ${JPEG_LIBRARY})
endif()

find_package(PNG REQUIRED)
set(HAVE_PKG_PNG PNG_FOUND)
if(PNGFOUND)
  message("Found PNG " ${PNG_LIBRARY})
endif()

find_package(TIFF REQUIRED)
set(HAVE_PKG_TIFF TIFF_FOUND)
if(TIFF_FOUND)
  message("Found TIFF " ${TIFF_LIBRARY})
endif()

find_package(GDAL REQUIRED)
set(HAVE_PKG_GDAL GDAL_FOUND)
if(GDAL_FOUND)
  message("Found GDAL at " ${GDAL_LIBRARY})
  include_directories(${GDAL_INCLUDE_DIR})
endif()

find_package(Protobuf REQUIRED)
set(HAVE_PKG_PROTOBUF PROTOBUF_FOUND)
if(PROTOBUF_FOUND)
  message("Found PROTOBUF " ${PROTOBUF_LIBRARY})
endif()

find_package(PythonLibs REQUIRED)
set(HAVE_PKG_PYTHONLIBS PYTHONLIBS_FOUND)
if(PYTHONLIBS_FOUND)
  message("Found PYTHONLIBS " ${PYTHON_LIBRARIES})
endif()


#TODO: Use FindProj cmake file to do this!
set(HAVE_PKG_PROJ4 true)
set(PROJ4_LIBRARY     ${BASESYSTEM_INSTALL_DIR}/lib/libproj.so)
set(PROJ4_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include)


set(HAVE_PKG_OPENEXR true)
set(OPENEXR_LIBRARIES   ${BASESYSTEM_INSTALL_DIR}/lib/libIlmImf.so
                        ${BASESYSTEM_INSTALL_DIR}/lib/libIlmThread.so )
set(OPENEXR_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include/OpenEXR)
include_directories("${OPENEXR_INCLUDE_DIR}")


#TODO: Look for each of these installations!
set(HAVE_PKG_HDF)
set(HAVE_PKG_CAIROMM)
set(HAVE_PKG_RABBITMQ_C)
set(HAVE_PKG_ZEROMQ)
set(HAVE_PKG_LIBKML)
set(HAVE_PKG_OPENCV)
set(HAVE_PKG_FLANN)

# TODO: What is this for?
set(ENABLE_EXCEPTIONS True)




message("--------------------------------------------------------------")
message("Listing all include directories...")
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()
message("--------------------------------------------------------------")



#=========================================================================
# Define each of the libaries
#TODO: Handle .lo files!


# --- VW_CAMERA ------------------------------------------------------------
get_all_source_files( "Camera" VW_CAMERA_SRC_FILES)
get_all_source_files( "Camera/tests" VW_CAMERA_TEST_FILES)
set(VW_CAMERA_LIB_DEPENDENCIES VwCore)

#message("${VW_CAMERA_SRC_FILES}")

# --- VW_CORE ------------------------------------------------------------
get_all_source_files( "Core" VW_CORE_SRC_FILES)
get_all_source_files( "Core/tests" VW_CORE_TEST_FILES)
set(VW_CORE_LIB_DEPENDENCIES )

# --- VW_FILEIO ------------------------------------------------------------
# This is a complicated folder so we need to select files manually
#get_all_source_files( "FileIO" VW_FILEIO_SRC_FILES)
#get_all_source_files( "FileIO/tests" VW_FILEIO_TEST_FILES)
set(VW_FILEIO_LIB_DEPENDENCIES VwCore)

# --- VW_IMAGE ------------------------------------------------------------
get_all_source_files( "Image" VW_IMAGE_SRC_FILES)
get_all_source_files( "Image/tests" VW_IMAGE_TEST_FILES)
set(VW_IMAGE_LIB_DEPENDENCIES VwCore)



#=========================================================================
# Add all of the library subdirectories

add_subdirectory(Camera)
add_subdirectory(Core)
add_subdirectory(FileIO)
add_subdirectory(Image)
add_subdirectory(tools)



#=========================================================================
# Define targets

# Targets are defined in /tools/CMakeLists.txt







