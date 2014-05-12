
# This file generates config.h from config.h.in

include(CheckIncludeFiles)
include(CheckSymbolExists)


###########################################################################
# These are some options that the user might like to set manually

#TODO: Where does all this get set?
set(DEFAULT_CACHE_SIZE_MB 768)

# set the default number of processing threads for multi-threaded operations
set(NUM_THREADS 4)

# enable image bounds checking (SLOW!) 
set(ENABLE_BOUND_CHECK false)

# ~/.vwrc support
set(ENABLE_CONFIG_FILE true) 

# enable the C++ exception mechanism 
set(ENABLE_EXCEPTIONS true) 

# enable SSE optimizations in some places (development) 
set(ENABLE_SSE true) 


/* Define to `int' if <sys/types.h> does not define. */
#define VW_ssize_t





###########################################################################

# Check if certain include files are present
# - Define to 1 if present, blank otherwise.
# - TODO: Probably need to set to zero if not present
CHECK_INCLUDE_FILES(ext/stdio_filebuf.h HAVE_EXT_STDIO_FILEBUF_H)
CHECK_INCLUDE_FILES(fenv.h              HAVE_FENV_H)
CHECK_INCLUDE_FILES(inttypes.h          HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(memory.h            HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(pwd.h               HAVE_PWD_H)
CHECK_INCLUDE_FILES(stdint.h            HAVE_STDINT_H)
CHECK_INCLUDE_FILES(stdlib.h            HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(strings.h           HAVE_STRINGS_H)
CHECK_INCLUDE_FILES(string.h            HAVE_STRING_H)
CHECK_INCLUDE_FILES(sys/stat.h          HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(sys/types.h         HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(dlfcn.h             HAVE_DLFCN_H)


/* Define to 1 if you have the ANSI C header files. */
#define VW_STDC_HEADERS @STDC_HEADERS@


###########################################################################
# Check if certain compiler features are available

set(emptyIncludeList )
CHECK_SYMBOL_EXISTS("__attribute__((deprecated))"         emptyIncludeList COMPILER_HAS_ATTRIBUTE_DEPRECATED)
CHECK_SYMBOL_EXISTS("__attribute__((noreturn))"           emptyIncludeList COMPILER_HAS_ATTRIBUTE_NORETURN)
CHECK_SYMBOL_EXISTS("__attribute__((warn_unused_result))" emptyIncludeList COMPILER_HAS_ATTRIBUTE_WARN_UNUSED_RESULT)

# Check for some supported functions
CHECK_SYMBOL_EXISTS("exp2"            emptyIncludeList HAVE_EXP2)
CHECK_SYMBOL_EXISTS("fabsl"           emptyIncludeList HAVE_FABSL)
CHECK_SYMBOL_EXISTS("feenableexcept2" emptyIncludeList HAVE_FEENABLEEXCEPT)
CHECK_SYMBOL_EXISTS("getpid"          emptyIncludeList HAVE_PID)
CHECK_SYMBOL_EXISTS("getpwuid"        emptyIncludeList HAVE_PWUID)
CHECK_SYMBOL_EXISTS("llabs"           emptyIncludeList HAVE_LLABS)
CHECK_SYMBOL_EXISTS("log2"            emptyIncludeList HAVE_LOG2)
CHECK_SYMBOL_EXISTS("mkstemps"        emptyIncludeList HAVE_MKSTEMPS)
CHECK_SYMBOL_EXISTS("tgamma"          emptyIncludeList HAVE_TGAMMA)
CHECK_SYMBOL_EXISTS("ssize_t"         "sys/types.h"    HAVE_SSIZET)





###########################################################################
# Determine which libraries we can build

# If we made it to here we can build these modules
# TODO: Do we need to check any other dependencies here?
set(HAVE_PKG_CORE)
set(HAVE_PKG_MATH)
set(HAVE_PKG_IMAGE)
set(HAVE_PKG_FILEIO)
set(HAVE_PKG_MATH)
set(HAVE_PKG_VW) # This is not actually a module
set(HAVE_PKG_CAMERA)
set(HAVE_PKG_INTERESTPOINT)
set(HAVE_PKG_CARTOGRAPHY)
set(HAVE_PKG_MOSAIC)
set(HAVE_PKG_HDR)
set(HAVE_PKG_STEREO)
set(HAVE_PKG_GEOMETRY)
set(HAVE_PKG_BUNDLEADJUSTMENT)

# The next few libraries are deprecated
if (false)
#if (${HAVE_PKG_RABBITMQ_C} and ${HAVE_PKG_ZEROMQ} and ${HAVE_PKG_LIBKML})
  set(HAVE_PKG_PLATE)
endif()

if (false)
#if (${HAVE_PKG_PLATE} and ${HAVE_PKG_QT})
  set(HAVE_PKG_GUI)
endif()

if (false)
#if (${HAVE_PKG_GL} and ${HAVE_PKG_GLEW})
  set(HAVE_PKG_GPU)
endif()


set(HAVE_PKG_TOOLS)





#######################################################################
# Finished setting up variables, now call the function to paste them into a file

# Each value like "@VAR@ is replaced by the CMake variable of the same name
configure_file(${CMAKE_SOURCE_DIR}/vw/config.in ${CMAKE_SOURCE_DIR}/vw/config.h)



###
###
###   TODO: These variables are used, try to find them!
###
###


#/* Define to 1 if VW has BigTIFF support */
##define VW_HAS_BIGTIFF @HAS_BIGTIFF@



#/* Define to 1 if the CG package is available. */
##define VW_HAVE_PKG_CG @HAVE_PKG_CG@


#/* Define to 1 if the HDR module is available. */
##define VW_HAVE_PKG_HDR @HAVE_PKG_HDR@

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


## Needed for plate
#///* Define to 1 if the RABBITMQ_C package is available. */
#//#define VW_HAVE_PKG_RABBITMQ_C

#/* Define to 1 if the LIBKML package is available. */
##define VW_HAVE_PKG_LIBKML @VW_HAVE_PKG_LIBKML@



#/* Define to 1 if the PTHREADS package is available. */
##define VW_HAVE_PKG_PTHREADS @VW_HAVE_PKG_PROTOBUF@


#//DELETE ME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#///* Define to 1 if the THREADS package is available. */
#//#define VW_HAVE_PKG_THREADS


#/* Define to 1 if the ZEROMQ package is available. */
##define VW_HAVE_PKG_ZEROMQ @VW_HAVE_PKG_ZEROMQ@













