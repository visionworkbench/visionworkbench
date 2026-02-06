# This file generates config.h from config.h.in

include(CheckIncludeFiles)
include(CheckSymbolExists)
include(CheckFunctionExists)
include(CheckTypeSize)
include(CheckCXXSymbolExists)
include(CheckCXXSourceCompiles)

###########################################################################

# Check if certain include files are present
# - Define to 1 if present, blank otherwise.
# - TODO: Probably need to set to zero if not present
CHECK_INCLUDE_FILES(ext/stdio_filebuf.h VW_HAVE_EXT_STDIO_FILEBUF_H) # TODO
CHECK_INCLUDE_FILES(fenv.h              VW_HAVE_FENV_H)
CHECK_INCLUDE_FILES(inttypes.h          VW_HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(memory.h            VW_HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(pwd.h               VW_HAVE_PWD_H)
CHECK_INCLUDE_FILES(stdint.h            VW_HAVE_STDINT_H)
CHECK_INCLUDE_FILES(stdlib.h            VW_HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(strings.h           VW_HAVE_STRINGS_H)
CHECK_INCLUDE_FILES(string.h            VW_HAVE_STRING_H)
CHECK_INCLUDE_FILES(sys/stat.h          VW_HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(sys/types.h         VW_HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(dlfcn.h             VW_HAVE_DLFCN_H)
CHECK_INCLUDE_FILES(unistd.h            VW_HAVE_UNISTD_H)

# Ignore, only used by plate.
# Define to 1 if you have the ANSI C header files. 
#define VW_STDC_HEADERS @STDC_HEADERS@

###########################################################################
# Check if certain compiler features are available

set(emptyIncludeList )

CHECK_CXX_SOURCE_COMPILES("void testFunc() __attribute__((deprecated));         void testFunc(){}   int main(){return 0;}" VW_COMPILER_HAS_ATTRIBUTE_DEPRECATED)
CHECK_CXX_SOURCE_COMPILES("void testFunc() __attribute__((noreturn));           void testFunc(){}   int main(){return 0;}" VW_COMPILER_HAS_ATTRIBUTE_NORETURN)
CHECK_CXX_SOURCE_COMPILES("void testFunc() __attribute__((warn_unused_result)); void testFunc(){}   int main(){return 0;}" VW_COMPILER_HAS_ATTRIBUTE_WARN_UNUSED_RESULT)


# Check for some supported functions (could probably streamline ssize_t check)
check_cxx_symbol_exists(exp2            cmath                  VW_HAVE_EXP2)
check_cxx_symbol_exists(fabsl           cmath                  VW_HAVE_FABSL)
check_cxx_symbol_exists(feenableexcept  "fenv.h"               VW_HAVE_FEENABLEEXCEPT)
check_cxx_symbol_exists(getpid          "unistd.h;sys/types.h" VW_HAVE_GETPID)
check_cxx_symbol_exists(getpwuid        "pwd.h;sys/types.h"    VW_HAVE_GETPWUID)
check_cxx_symbol_exists(llabs           "stdlib.h"             VW_HAVE_LLABS)
check_cxx_symbol_exists(log2            cmath                  VW_HAVE_LOG2)
check_cxx_symbol_exists(mkstemps        "stdlib.h"             VW_HAVE_MKSTEMPS)
check_cxx_symbol_exists(tgamma          cmath                  VW_HAVE_TGAMMA)
CHECK_CXX_SOURCE_COMPILES("
                          #include <sys/types.h>
                          int main(){ssize_t a=2; return a;}" VW_HAVE_SSIZET)




###########################################################################
# Determine which libraries we can build

# If we made it to here we can build these modules
set(VW_HAVE_PKG_CORE 1)
set(VW_HAVE_PKG_MATH 1)
set(VW_HAVE_PKG_IMAGE 1)
set(VW_HAVE_PKG_FILEIO 1)
set(VW_HAVE_PKG_MATH 1)
set(VW_HAVE_PKG_VW 1) # This is not actually a module
set(VW_HAVE_PKG_CAMERA 1)
set(VW_HAVE_PKG_INTERESTPOINT 1)
set(VW_HAVE_PKG_CARTOGRAPHY 1)
set(VW_HAVE_PKG_MOSAIC 1)
set(VW_HAVE_PKG_HDR ${VW_ENABLE_HDR})
set(VW_HAVE_PKG_STEREO 1)
set(VW_HAVE_PKG_GEOMETRY 1)
set(VW_HAVE_PKG_BUNDLEADJUSTMENT 1)

set(VW_HAVE_PKG_TOOLS 1)

#######################################################################
# Finished setting up variables, now call the function to paste them into files

# Generate stable config file (feature detection - rarely changes)
message("Generating config file: ${CMAKE_BINARY_DIR}/src/vw/vw_config.h")
configure_file(${CMAKE_SOURCE_DIR}/src/vw/vw_config.h.in ${CMAKE_BINARY_DIR}/src/vw/vw_config.h)

# Generate config file with elements that change often 
message("Generating config file: ${CMAKE_BINARY_DIR}/src/vw/vw_date_config.h")
configure_file(${CMAKE_SOURCE_DIR}/src/vw/vw_date_config.h.in ${CMAKE_BINARY_DIR}/src/vw/vw_date_config.h)


###
###
###   TODO: These variables are used, try to find them!
###
###


#/* Define to 1 if VW has BigTIFF support */
##define VW_HAS_BIGTIFF @HAS_BIGTIFF@



#/* Define to 1 if the CG package is available. */
##define VW_HAVE_PKG_CG @HAVE_PKG_CG@


#/* Define to 1 if the HDF package is available */
##define VW_HAVE_PKG_HDF @HAVE_PKG_HDF@

#/* Define to 1 if the INTEL_LAPACK package is available. */
##define VW_HAVE_PKG_INTEL_LAPACK @HAVE_PKG_INTEL_LAPACK@


#/* Define to 1 if the CLAPACK package is available. */
##define VW_HAVE_PKG_CLAPACK @VW_HAVE_PKG_CLAPACK@

#///* Define to 1 if the pkg package is available */
#//#define VW_HAVE_PKG_APPLE_LAPACK @PKG_APPLE_LAPACK@

#/* Define to 1 if the SLAPACK package is available. */
##define VW_HAVE_PKG_SLAPACK @VW_HAVE_PKG_SLAPACK@

#///* Define to 1 if the STANDALONE_FLAPACK package is available. */
#//#define VW_HAVE_PKG_STANDALONE_FLAPACK

#///* Define to 1 if the pkg package is available */
#//#define VW_HAVE_PKG_STANDALONE_LAPACK_AND_BLAS


#/* Define to 1 if the LIBKML package is available. */
##define VW_HAVE_PKG_LIBKML @VW_HAVE_PKG_LIBKML@



#/* Define to 1 if the PTHREADS package is available. */
##define VW_HAVE_PKG_PTHREADS @VW_HAVE_PKG_PROTOBUF@ #TODO?


#//DELETE ME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#///* Define to 1 if the THREADS package is available. */
#//#define VW_HAVE_PKG_THREADS













