
# TODO: Need to move this up a level to handle other release stuff like documentation


cmake_minimum_required (VERSION 2.8)


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
function( add_library_wrapper libName fileList testFileList)

  #message("${CMAKE_CURRENT_SOURCE_DIR}")

  # Set up the library
  add_library(${libName} ${fileList})

  set_target_properties("${libName}" PROPERTIES LINKER_LANGUAGE CXX)   
  install(TARGETS "${libName}" DESTINATION lib)
  
  # Set all the header files to be installed to the include directory
  foreach(f ${fileList})
    get_filename_component(extension ${f} EXT) # Get file extension  
    string( TOLOWER "${extension}" extensionLower )
    if( extensionLower STREQUAL ".h" OR extensionLower STREQUAL ".hpp" OR extensionLower STREQUAL ".tcc")
      #install(FILES "${f}" DESTINATION include)
      #INSTALL_HEADERS_WITH_DIRECTORY("${f}")
      
      set(fullPath "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
      STRING(REGEX MATCH "vw/.*/" dir ${fullPath})
      #message("${dir}")
      INSTALL(FILES ${f} DESTINATION include/${dir})
      
    endif()
  endforeach(f)
  


  # Add unit test for each test file given
  set(TEST_MAIN_PATH "${CMAKE_SOURCE_DIR}/test/test_main.cc")
  foreach(f ${testFileList})
  
    get_filename_component(filename ${f} NAME_WE) # Get file name without extension
    set(executableName "${libName}_${filename}")   # Generate a name for the executable   
    
    #message("Adding test target ${executableName}")
    
    # Add executable with shared main file and this file
    add_executable( ${executableName}   ${TEST_MAIN_PATH} ./tests/${f} )      

    # Link test executable against current library, gtest, and gtest_main
    #target_link_libraries(${executableName} gtest "${libName}" ${GTEST_BOTH_LIBRARIES})
    target_link_libraries(${executableName} gtest gtest_main ${FULL_LIBRARY_LIST})

    
    # These variables need to be set for each test directory: TODO: What should they be?
    #  TEST_OBJDIR = "\"$(abs_top_builddir)/$(subdir)\""
    #  TEST_SRCDIR = "\"$(abs_top_srcdir)/$(subdir)\""
    set_property (TARGET ${executableName} APPEND PROPERTY COMPILE_DEFINITIONS "TEST_OBJDIR=\"${CMAKE_CURRENT_SOURCE_DIR}\"")
    set_property (TARGET ${executableName} APPEND PROPERTY COMPILE_DEFINITIONS "TEST_SRCDIR=\"${CMAKE_CURRENT_SOURCE_DIR}\"")
    
    add_test(${executableName} ${executableName}) 
    add_to_custom_test_target(${executableName})  # Add to the verbose test make target.
  endforeach(f)
  
  
endfunction( add_library_wrapper )




#===============================================================================
# Project information

project (VisionWorkbench)

set(PACKAGE_VERSION "1.0") # TODO: Automatically set this!
set(PACKAGE_NAME "Vision Workbench")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_BUGREPORT "scott.t.mcmichael@nasa.gov")

#===============================================================================

# Set the install path
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/../cmakeBuild/install" CACHE INTERNAL "Prefix
prepended to install directories")
message("Cmake install directory: ${CMAKE_INSTALL_PREFIX}")

enable_testing() # Needed for unit tests to work

set(SP_STATIC_LIBRARIES_LIST "")

# Add command to build the instructions pdf file
# - First delete temporary qtree file if it exists
# - Next run the doc-generate tool (which must be built), generates files in /image directory.
#   This step is really slow!
# - Finally call pdflatex to generate the output file
add_custom_target(workbook rm -rf images/Walker.qtree 
                  COMMAND ${CMAKE_BINARY_DIR}/vw/tools/doc-generate
                  COMMAND pdflatex workbook.tex
                  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/../docs/workbook")


# Add command to build doxygen files
add_custom_target(doxygen doxygen
                  COMMAND ln -s ${CMAKE_SOURCE_DIR}/../docs/doxygen/index.html ${CMAKE_SOURCE_DIR}/../docs/doxygen-output.html
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






