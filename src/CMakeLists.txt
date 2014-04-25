
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
  )
  set(fileList) # Empty list
  foreach(f ${TEMP}) # Iterate through TEMP
    get_filename_component(FILENAME ${f} NAME) # Extract just the file name
    set(fileList ${fileList} ${FILENAME}) # Append to the list
  endforeach(f)
  set(${outputFileList} ${fileList} PARENT_SCOPE) 
endfunction(get_all_source_files)


# Function to add a library to the project.
# - This is called in each library folder directory.
function( add_library_wrapper libName fileList testFileList)

  # Set up the library
  add_library(${libName} ${fileList})

  set_target_properties("${libName}" PROPERTIES LINKER_LANGUAGE CXX)   
 
  
endfunction( add_library_wrapper )



#===============================================================================
# Project information

project (VisionWorkbench)

set(PACKAGE_VERSION "1.0") # TODO: Automatically set this!
set(PACKAGE_NAME "Vision Workbench")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_BUGREPORT "scott.t.mcmichael@nasa.gov")

#===============================================================================

enable_testing() # Needed for unit tests to work

set(SP_STATIC_LIBRARIES_LIST "")


# Add external dependencies which still need to be built

# Build GTest
#include_directories(../thirdparty/gtest/include)
#include_directories(../thirdparty/gtest/)
#add_library(gtest      ../thirdparty/gtest/src/gtest-all.cc)
#add_library(gtest_main ../thirdparty/gtest/src/gtest_main.cc)
#target_link_libraries(gtest_main gtest) 



# Add all of our code
add_subdirectory(vw)






