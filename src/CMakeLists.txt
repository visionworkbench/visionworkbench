



# TODO: Move everything up a directory!

include(Utilities) # Contains functions

#===============================================================================
# Project information

# Need this to make sure we get the compiler from PATH.
# - cc and c++ both need to be on the PATH.
find_program(CMAKE_C_COMPILER   NAMES $ENV{CC}  gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

#project (VisionWorkbench C CXX CXXPCH CPCH)
project (VisionWorkbench C CXX)

# Verify that we have a new enough compiler.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.5")
    message(FATAL_ERROR "GCC version >= 4.8.5 is required!")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "8.1.0")
    message(FATAL_ERROR "Clang version >= 8.1.0 is required!")
  endif()
else()
  message(FATAL_ERROR "Unknown compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()

# Short and long name of this package
set(PACKAGE            "VisionWorkbench")
set(PACKAGE_NAME       "NASA Vision Workbench")

# Version number of the package
set(VERSION            "2.2.0_post")
set(PACKAGE_VERSION    ${VERSION})

# Full name and version number
set(PACKAGE_STRING     "${PACKAGE_NAME} ${VERSION}")

# Define to the address where bug reports for this package should be sent. 
set(PACKAGE_BUGREPORT  "vision-workbench@lists.nasa.gov")

# Define to the one symbol short name of this package.
set(PACKAGE_TARNAME    "${PACKAGE_NAME}")

# Main website associated with the software package
set(PACKAGE_URL        "http://ti.arc.nasa.gov/tech/asr/intelligent-robotics/nasa-vision-workbench/")


# Fetch the Git hash and store in VW_COMMIT_ID
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE VW_COMMIT_ID
  OUTPUT_STRIP_TRAILING_WHITESPACE
)


#===============================================================================

# Set the install path if not specified by the user.
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "default install path" FORCE ) 
endif()
message("Cmake install directory set to: ${CMAKE_INSTALL_PREFIX}")

enable_testing() # Needed for unit tests to work

set(VW_STATIC_LIBRARIES_LIST "")

# Add command to build the instructions pdf file
# - First delete temporary qtree file if it exists
# - Next run the doc-generate tool (which must be built), generates files in /image directory.
#   This step is really slow!
# - Finally call pdflatex to generate the output file
add_custom_target(workbook rm -rf images/Walker.qtree 
                  COMMAND ${CMAKE_BINARY_DIR}/vw/tools/doc-generate
                  COMMAND pdflatex workbook.tex
                  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/docs/workbook")


# Add command to build doxygen files
add_custom_target(doxygen doxygen
                  COMMAND ln -s ${CMAKE_SOURCE_DIR}/docs/doxygen/index.html ${CMAKE_SOURCE_DIR}/docs/doxygen-output.html
                  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")



# Add external dependencies which still need to be built

# Build GTest
include_directories(../thirdparty/gtest/include)
#include_directories(../thirdparty/gtest/)
add_library(gtest      ../thirdparty/gtest/src/gtest-all.cc)
add_library(gtest_main ../thirdparty/gtest/src/gtest_main.cc)
target_link_libraries(gtest_main gtest) 

target_compile_definitions(gtest PRIVATE GTEST_USE_OWN_TR1_TUPLE=1)
target_compile_definitions(gtest_main PRIVATE GTEST_USE_OWN_TR1_TUPLE=1)




# Add all of our code
add_subdirectory(vw)







#========================================================================================


# Deprecated - We no longer use protobuf files!
function(generate_protobuf_files PROTO_GEN_OUT)

  # Finds all .proto files in the current dir
  file(GLOB PROTO_INPUT RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.proto")

  # Protobuf compiler
  find_program(PROTOC protoc)
  set(PROTOC_C_OUT_FLAG --cpp_out)

  # Set where generated files go to and add that directory to the include path
  set(PROTO_GEN_DIR ${CMAKE_CURRENT_BINARY_DIR})
  include_directories(${PROTO_GEN_DIR})

  # For each input protobuf file
  foreach(PROTO_FILE ${PROTO_INPUT})
    # Get the name of the file without extension
    get_filename_component(PROTO_NAME ${PROTO_FILE} NAME_WE)
    
    # Add the two generated files to PROTO_GEN variable
    set(PROTO_GEN ${PROTO_GEN}
        ${PROTO_GEN_DIR}/${PROTO_NAME}.pb.h
        ${PROTO_GEN_DIR}/${PROTO_NAME}.pb.cc)
  endforeach()

  # Add the custom command that will generate all the files
  # - The generated files will be put in the CMake build directory, not the source tree.
  #message("Output files: ${PROTO_GEN}")
  add_custom_command(OUTPUT   ${PROTO_GEN}
                     COMMAND  ${PROTOC} ${PROTO_INPUT} ${PROTOC_C_OUT_FLAG} ${CMAKE_CURRENT_BINARY_DIR}
                     DEPENDS  ${PROTO_INPUT}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     COMMENT "Generating Protocol Buffers...")

  set(${PROTO_GEN_OUT} ${PROTO_GEN} PARENT_SCOPE) # Set up output variable
endfunction(generate_protobuf_files)



