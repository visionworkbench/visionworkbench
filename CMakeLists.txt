cmake_minimum_required(VERSION 3.8)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_DEBUG_POSTFIX -gd)
set(CMAKE_CXX_STANDARD 11)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

project(VisionWorkbench VERSION 2.6.0)

option(BUILD_SHARED_LIBS "build shared libs" ON)
if(BUILD_SHARED_LIBS)
    add_definitions(-DVW_SHARED_LIB -DVW_EXPORT)
endif(BUILD_SHARED_LIBS)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /MP /wd4305 /wd4018 /wd4715 /wd4800 /wd4819 /wd4101 /wd4251 /wd4996 /wd4275 /wd4244 /wd4267")
endif(MSVC)

option(VW_ENABLE_EXCEPTIONS "" ON)
option(VW_ENABLE_SSE "" ON)
option(VW_HAVE_PKG_CORE "" ON)
option(VW_HAVE_PKG_MATH "" ON)
option(VW_HAVE_PKG_IMAGE "" ON)
option(VW_HAVE_PKG_FILEIO "" ON)
option(VW_HAVE_PKG_CAMERA "" ON)
option(VW_HAVE_PKG_MOSAIC "" ON)
option(VW_HAVE_PKG_HDR "" ON)
option(VW_HAVE_PKG_CARTOGRAPHY "" ON)
option(VW_HAVE_PKG_GEOMETRY "" ON)
option(VW_HAVE_PKG_GPU "" OFF)
option(VW_HAVE_PKG_INTERESTPOINT "" ON)
option(VW_HAVE_PKG_STEREO "" ON)
option(VW_HAVE_PKG_BUNDLEADJUSTMENT "" ON)
option(VW_HAVE_PKG_TOOL "" ON)
option(VW_HAVE_PKG_FLOOD_DETECT "" ON)
# obliged libs
if(MSVC)
    set(Boost_USE_STATIC_LIBS ON)
endif(MSVC)
find_package(Boost REQUIRED COMPONENTS program_options system thread filesystem iostreams)
set(VW_HAVE_PKG_BOOST_IOSTREAMS 1)

# image io
find_package(JPEG REQUIRED)
find_package(PNG REQUIRED)
find_package(GDAL REQUIRED)
find_package(TIFF REQUIRED)
set(VW_HAVE_PKG_JPEG 1)
set(VW_HAVE_PKG_PNG 1)
set(VW_HAVE_PKG_GDAL 1)
set(VW_HAVE_PKG_TIFF 1)

# opencv
option(VW_HAVE_PKG_OPENCV "" ON)
find_package(OpenCV REQUIRED) # OpenCV_LIBS OpenCV_INCLUDE_DIRS 
set_target_properties(${OpenCV_LIBS} PROPERTIES MAP_IMPORTED_CONFIG_RELWITHDEBINFO Release)

# blas and lapack
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
set(VW_HAVE_PKG_STANDALONE_LAPACK_AND_BLAS 1)

set(VW_3RDPARTY_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR})
include_directories(${VW_3RDPARTY_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})

# flann
find_package(Flann REQUIRED) #
set(VW_HAVE_PKG_FLANN 1)

# proj4
find_package(Proj REQUIRED)
set(VW_HAVE_PKG_PROJ4 1)

add_subdirectory(src)

# configure file
set(VW_LIBRARIES "")

if(VW_ENABLE_EXCEPTIONS)
    set(VW_ENABLE_EXCEPTIONS 1)
endif(VW_ENABLE_EXCEPTIONS)
if(VW_HAVE_PKG_CORE)
    set(VW_HAVE_PKG_CORE 1)
    list(APPEND VW_LIBRARIES vw::vw_core)
endif()
if(VW_HAVE_PKG_MATH)
    set(VW_HAVE_PKG_MATH 1)
    list(APPEND VW_LIBRARIES vw::vw_math)    
endif()
if(VW_HAVE_PKG_IMAGE)
    set(VW_HAVE_PKG_IMAGE 1)
    list(APPEND VW_LIBRARIES vw::vw_image)    
endif()
if(VW_HAVE_PKG_FILEIO)
    set(VW_HAVE_PKG_FILEIO 1)
    list(APPEND VW_LIBRARIES vw::vw_fileio)    
endif()
if(VW_HAVE_PKG_CAMERA)
    set(VW_HAVE_PKG_CAMERA 1)
    list(APPEND VW_LIBRARIES vw::vw_camera)    
endif()
if(VW_HAVE_PKG_MOSAIC)
    set(VW_HAVE_PKG_MOSAIC 1)
    list(APPEND VW_LIBRARIES vw::vw_mosaic)    
endif()
if(VW_HAVE_PKG_CARTOGRAPHY)
    set(VW_HAVE_PKG_CARTOGRAPHY 1)
    list(APPEND VW_LIBRARIES vw::vw_cartography)    
endif()
if(VW_HAVE_PKG_GEOMETRY)
    set(VW_HAVE_PKG_GEOMETRY 1)
    list(APPEND VW_LIBRARIES vw::vw_geometry)    
endif()
if(VW_HAVE_PKG_GPU)
    set(VW_HAVE_PKG_GPU 1)
endif()
if(VW_HAVE_PKG_INTERESTPOINT)
    set(VW_HAVE_PKG_INTERESTPOINT 1)
    list(APPEND VW_LIBRARIES vw::vw_interestpoint)        
endif()
if(VW_HAVE_PKG_STEREO)
    set(VW_HAVE_PKG_STEREO 1)
    list(APPEND VW_LIBRARIES vw::vw_stereo)    
endif()
if(VW_HAVE_PKG_BUNDLEADJUSTMENT)
    set(VW_HAVE_PKG_BUNDLEADJUSTMENT 1)
    list(APPEND VW_LIBRARIES vw::vw_bundleadjustment)    
endif(VW_HAVE_PKG_BUNDLEADJUSTMENT)
if(VW_HAVE_PKG_HDR)
    set(VW_HAVE_PKG_HDR 1)
    list(APPEND VW_LIBRARIES vw::vw_hdr)    
endif(VW_HAVE_PKG_HDR)

if(VW_HAVE_PKG_OPENCV)
    set(VW_HAVE_PKG_OPENCV 1)
endif(VW_HAVE_PKG_OPENCV)
if(VW_ENABLE_SSE)
    set(VW_ENABLE_SSE 1)
endif(VW_ENABLE_SSE)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/vw/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/vw/config.h @ONLY)

# settings for installation
set(config_install_dir "lib/cmake/vw")
set(include_install_dir "include")
set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(version_config "${generated_dir}/vw-config-version.cmake")
set(project_config "${generated_dir}/vw-config.cmake")
set(TARGETS_EXPORT_NAME "vw-targets")
set(namespace "vw::")

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${version_config}" COMPATIBILITY SameMajorVersion VERSION ${PROJECT_VERSION}
)

configure_package_config_file(
    "cmake/config.cmake.in"
    "${project_config}"
    INSTALL_DESTINATION "${config_install_dir}"
)

# config.h
install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/vw/config.h
    DESTINATION "${include_install_dir}/vw"
)

file(GLOB VW_ADDITION_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/src/vw/*.h)

install(
    FILES ${VW_ADDITION_HEADERS}
    DESTINATION "${include_install_dir}/vw"
)

install(
    FILES "${project_config}" "${version_config}"
    DESTINATION "${config_install_dir}"
)

install(
    EXPORT "${TARGETS_EXPORT_NAME}"
    NAMESPACE "${namespace}"
    DESTINATION "${config_install_dir}"
)