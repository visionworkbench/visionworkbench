include(Utilities) # Contains functions

# Need this to make sure we get the compiler from PATH.
# - cc and c++ both need to be on the PATH.
find_program(CMAKE_C_COMPILER   NAMES $ENV{CC}  gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project (VisionWorkbench C CXX)

# Verify that we have a new enough compiler.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.0.0")
    message(FATAL_ERROR "GCC version >= 5.0.0 is required!")
  endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "8.1.0")
    message(FATAL_ERROR "Clang version >= 8.1.0 is required!")
  endif()
else()
  message(FATAL_ERROR "Unknown compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()

# Short and long name of this package
set(PACKAGE            "VisionWorkbench")
set(PACKAGE_NAME       "NASA Vision Workbench")

# The line below is used by the docs/ so don't change it
# without testing the result.  Versioning is based on https://semver.org
#set(PACKAGE_VERSION "3.5.0")
set(PACKAGE_VERSION "3.6.0-alpha")

# Full name and version number
set(PACKAGE_STRING     "${PACKAGE_NAME} ${PACKAGE_VERSION}")

# Define to the address where bug reports for this package should be sent. 
set(PACKAGE_BUGREPORT  "ames-stereo-pipeline-support@googlegroups.com")

# Define to the one symbol short name of this package.
set(PACKAGE_TARNAME "${PACKAGE_NAME}")

# Main website associated with the software package
set(PACKAGE_URL "https://github.com/NeoGeographyToolkit/StereoPipeline")

# Fetch the Git hash and store in VW_COMMIT_ID
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE VW_COMMIT_ID
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

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
