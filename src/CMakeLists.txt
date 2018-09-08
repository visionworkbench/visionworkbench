



# TODO: Move everything up a directory!

# Define functions --> Move these to another file


# Obtains a file list with all the files in a directory properly formatted
function( get_all_source_files relativePath outputFileList)

  # Load all matching files into TEMP
  file(GLOB TEMP
      "${CMAKE_CURRENT_SOURCE_DIR}/${relativePath}/*.h"
      "${CMAKE_CURRENT_SOURCE_DIR}/${relativePath}/*.cc"
      "${CMAKE_CURRENT_SOURCE_DIR}/${relativePath}/*.cxx"      
      "${CMAKE_CURRENT_SOURCE_DIR}/${relativePath}/*.tcc"
  )
  set(fileList) # Empty list
  foreach(f ${TEMP}) # Iterate through TEMP
    get_filename_component(FILENAME ${f} NAME) # Extract just the file name
    set(fileList ${fileList} ${FILENAME}) # Append to the list
  endforeach(f)
  set(${outputFileList} ${fileList} PARENT_SCOPE) 
endfunction(get_all_source_files)

# Look for a library dependency, starting with the BinaryBuilder folder.
function(find_external_library name bbIncludeFolder libNameList required)

  # Define the variable names we will create
  set(FOUND_NAME "${name}_FOUND")
  set(LIB_NAME   "${name}_LIBRARIES")
  set(INC_NAME   "${name}_INCLUDE_DIR")
  set(VW_NAME    "VW_HAVE_PKG_${name}")

  # Look in the BB directory if it was provided, otherwise
  #  make halfhearted attempt to find the dependency.
  if(BINARYBUILDER_INSTALL_DIR)
    set(${FOUND_NAME} 1)

    # Add each lib file that was provided.
    set(${${LIB_NAME}} "")
    foreach(lib ${libNameList})
      set(${LIB_NAME}   ${${LIB_NAME}} ${BINARYBUILDER_INSTALL_DIR}/lib/${lib})
    endforeach()
    
    set(${INC_NAME}   ${BINARYBUILDER_INSTALL_DIR}/include/${bbIncludeFolder})
  else()
    # TODO: Provide effective findX.cmake files to handle these.
    find_package(${name} REQUIRED)
  endif()
  # Check and display our results
  if(${FOUND_NAME})
    set(${VW_NAME} 1)
    message("-- Found ${name} at " ${${INC_NAME}})
    include_directories("${${INC_NAME}}")
  else()
    if (${required})
      message( FATAL_ERROR "Failed to find REQUIRED library ${name}." )
    else()
      message("Failed to find ${name}")
    endif()
  endif()

  # Pass the results back up to the parent function
  set(${FOUND_NAME} ${${FOUND_NAME}} PARENT_SCOPE)
  set(${LIB_NAME}   ${${LIB_NAME}}   PARENT_SCOPE)
  set(${INC_NAME}   ${${INC_NAME}}   PARENT_SCOPE)
  set(${VW_NAME}    ${${VW_NAME}}    PARENT_SCOPE)

endfunction(find_external_library)

# Define a custom make target that will run all tests with normal gtest output.
# - Normally you can run 'make test' to run all tests but the output is brief.
# - With this you can run 'make gtest_all' to run all tests with more output.
add_custom_target(gtest_all)
# Call this function once for each gtest target.
macro(add_to_custom_test_target test_target)
  add_custom_target(${test_target}_runtest
                    COMMAND ${test_target} #cmake 2.6 required
                    DEPENDS ${test_target}
                    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  add_dependencies(gtest_all ${test_target}_runtest)
endmacro()



# Function to add a library to the project.
# - This is called in each library folder directory.
function(add_library_wrapper libName fileList testFileList dependencyList)

  # Set up the library
  add_library(${libName} ${fileList})

  set_target_properties(${libName} PROPERTIES LINKER_LANGUAGE CXX)   
  #message("For ${libName}, linking DEPS: ${dependencyList}")
  target_link_libraries(${libName} "${dependencyList}")

  install(TARGETS ${libName} DESTINATION lib)

  # Set all the header files to be installed to the include directory
  foreach(f ${fileList})
    get_filename_component(extension ${f} EXT) # Get file extension  
    string( TOLOWER "${extension}" extensionLower )
    if( extensionLower STREQUAL ".h" OR extensionLower STREQUAL ".hpp" OR extensionLower STREQUAL ".tcc")
      set(fullPath "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
      STRING(REGEX MATCH "vw/.*/" dir ${fullPath})
      INSTALL(FILES ${f} DESTINATION include/${dir})
    endif()
  endforeach(f)


  # Add unit test for each test file given
  set(TEST_MAIN_PATH "${CMAKE_SOURCE_DIR}/src/test/test_main.cc")
  foreach(f ${testFileList})

    get_filename_component(filename ${f} NAME_WE) # Get file name without extension
    set(executableName "${libName}_${filename}")   # Generate a name for the executable   

    #message("Adding test target ${executableName}")

    # Add executable with shared main file and this file
    # - This executeable should not be built unless running tests.
    add_executable( ${executableName} EXCLUDE_FROM_ALL  ${TEST_MAIN_PATH} ./tests/${f} )      

    # Link test executable against current library, gtest, and gtest_main
    #target_link_libraries(${executableName} gtest "${libName}" ${GTEST_BOTH_LIBRARIES})
    target_link_libraries(${executableName} gtest gtest_main ${FULL_LIBRARY_LIST})


    # These variables need to be set for each test directory
    set_property (TARGET ${executableName} APPEND PROPERTY COMPILE_DEFINITIONS "TEST_OBJDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/tests\"")
    set_property (TARGET ${executableName} APPEND PROPERTY COMPILE_DEFINITIONS "TEST_SRCDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/tests\"")

    add_test(${executableName} ${executableName}) 
    add_to_custom_test_target(${executableName})  # Add to the verbose test make target.
  endforeach(f)

endfunction( add_library_wrapper )




#===============================================================================
# Project information

# Need this to make sure we get the compiler from PATH.
# - cc and c++ both need to be on the PATH.
find_program(CMAKE_C_COMPILER   NAMES $ENV{CC}  gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project (VisionWorkbench C CXX)

# Verify that we have a new enough compiler.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9.3")
    message(FATAL_ERROR "GCC version >= 4.9.3 is required!")
  endif()
else()
  message(FATAL_ERROR "TODO: Get MacOS compiler version!")
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


#===============================================================================

# Set the install path
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/build/install" CACHE PATH "default install path" FORCE )
endif()
#set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/build/install" CACHE INTERNAL "Prefix prepended to install directories")
message("Cmake install directory: ${CMAKE_INSTALL_PREFIX}")

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



