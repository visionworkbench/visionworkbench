# Cross-compile toolchain for building on Mac ARM64 targeting Linux x86_64.
#
# Both VisionWorkbench and StereoPipeline have a copy of this file. Keep them
# in sync when making changes.
#
# Prerequisites:
#   - conda-forge clang 16+ with -fopenmp support (in MAC_ASP_DEPS env)
#   - lld (LLVM linker) in MAC_ASP_DEPS env
#   - Linux deps prefix with sysroot and GCC libraries (version auto-detected)
#
# Usage (VW, from build_linux/):
#   cmake .. \
#     -DCMAKE_TOOLCHAIN_FILE=../cmake/linux_cross_toolchain.cmake \
#     -DLINUX_DEPS_PREFIX=$HOME/miniconda3/envs/asp_deps_linux \
#     -DMAC_ASP_DEPS=$HOME/anaconda3/envs/asp_deps \
#     -DASP_DEPS_DIR=$HOME/miniconda3/envs/asp_deps_linux \
#     -DCMAKE_INSTALL_PREFIX=$HOME/projects/StereoPipeline/install_linux
#
# Usage (ASP, from build_linux/):
#   cmake .. \
#     -DCMAKE_TOOLCHAIN_FILE=../cmake/linux_cross_toolchain.cmake \
#     -DLINUX_DEPS_PREFIX=$HOME/miniconda3/envs/asp_deps_linux \
#     -DMAC_ASP_DEPS=$HOME/anaconda3/envs/asp_deps \
#     -DASP_DEPS_DIR=$HOME/miniconda3/envs/asp_deps_linux \
#     -DVISIONWORKBENCH_INSTALL_DIR=$HOME/projects/StereoPipeline/install_linux \
#     -DCMAKE_INSTALL_PREFIX=$HOME/projects/StereoPipeline/install_linux \
#     -DOpenMP_C_FLAGS=-fopenmp \
#     -DOpenMP_CXX_FLAGS=-fopenmp \
#     -DOpenMP_C_LIB_NAMES=omp \
#     -DOpenMP_CXX_LIB_NAMES=omp \
#     -DOpenMP_omp_LIBRARY=${LINUX_DEPS_PREFIX}/lib/libomp.so

# Propagate these variables to try_compile() projects. Without this, CMake's
# ABI detection re-invokes this toolchain file but loses the -D variables.
list(APPEND CMAKE_TRY_COMPILE_PLATFORM_VARIABLES LINUX_DEPS_PREFIX MAC_ASP_DEPS)

# Validate required variables.
if(NOT DEFINED LINUX_DEPS_PREFIX)
  message(FATAL_ERROR "Set -DLINUX_DEPS_PREFIX=/path/to/asp_deps_linux")
endif()
if(NOT DEFINED MAC_ASP_DEPS)
  message(FATAL_ERROR "Set -DMAC_ASP_DEPS=/path/to/mac/asp_deps")
endif()

# Derived paths.
set(CROSS_SYSROOT "${LINUX_DEPS_PREFIX}/x86_64-conda-linux-gnu/sysroot")
# Auto-detect GCC version in the linux prefix.
file(GLOB GCC_VERSION_DIRS "${LINUX_DEPS_PREFIX}/lib/gcc/x86_64-conda-linux-gnu/*")
list(GET GCC_VERSION_DIRS 0 GCC_LIB)
message(STATUS "Cross-compile: GCC lib dir = ${GCC_LIB}")
set(GCC_INC "${GCC_LIB}/include/c++")

# Target platform.
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR x86_64)

# Compilers (conda-forge clang from the Mac env).
set(CMAKE_C_COMPILER "${MAC_ASP_DEPS}/bin/clang")
set(CMAKE_CXX_COMPILER "${MAC_ASP_DEPS}/bin/clang++")

# Compiler flags. Pass --sysroot directly in flags instead of using
# CMAKE_SYSROOT, which strips path prefixes and corrupts other paths.
set(CROSS_COMMON_FLAGS
  "--target=x86_64-unknown-linux-gnu \
   --sysroot=${CROSS_SYSROOT} \
   --gcc-toolchain=${LINUX_DEPS_PREFIX} \
   -B${GCC_LIB} \
   -fuse-ld=lld \
   -I${CROSS_SYSROOT}/usr/include \
   -I${LINUX_DEPS_PREFIX}/include \
   -I${LINUX_DEPS_PREFIX}/include/eigen3 \
   -L${GCC_LIB} \
   -L${LINUX_DEPS_PREFIX}/lib")

set(CMAKE_C_FLAGS_INIT "${CROSS_COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS_INIT
  "${CROSS_COMMON_FLAGS} \
   -isystem ${GCC_INC} \
   -isystem ${GCC_INC}/x86_64-conda-linux-gnu \
   -Wno-enum-constexpr-conversion")

# Linker flags. Use lld from the Mac env, link libstdc++ from GCC toolchain,
# allow undefined symbols in shared libs (resolved at runtime on target).
set(CROSS_LINKER_FLAGS
  "-fuse-ld=lld \
   -L${GCC_LIB} -L${LINUX_DEPS_PREFIX}/lib \
   -lstdc++ -Wl,--allow-shlib-undefined")
set(CMAKE_EXE_LINKER_FLAGS_INIT "${CROSS_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS_INIT "${CROSS_LINKER_FLAGS}")

# Force --allow-shlib-undefined globally. Conda's pre-built Linux .so files
# reference symbols (e.g. __cxa_call_terminate@CXXABI_1.3.15) that are only
# resolved at runtime via libstdc++. Without this flag lld rejects them.
add_link_options("-Wl,--allow-shlib-undefined")

# Threads: FindThreads try_compile fails when cross-compiling (can't run
# the test binary). Force pthreads since Linux always has them.
set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)

# Search only the cross-prefix for libraries and headers.
set(CMAKE_FIND_ROOT_PATH "${LINUX_DEPS_PREFIX}")
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

# Qt moc/rcc/uic are host tools that must run on Mac, not target Linux binaries.
set(QT_HOST_PATH "${MAC_ASP_DEPS}")
set(QT_MOC_EXECUTABLE "${MAC_ASP_DEPS}/bin/moc")
set(QT_RCC_EXECUTABLE "${MAC_ASP_DEPS}/bin/rcc")
set(QT_UIC_EXECUTABLE "${MAC_ASP_DEPS}/bin/uic")
