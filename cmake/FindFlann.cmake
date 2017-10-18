# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#.rst:
# FindFlann
# --------
#
# Find the Flann library (libFlann).
#
# Imported targets
# ^^^^^^^^^^^^^^^^
#
# This module defines the following :prop_tgt:`IMPORTED` targets:
#
# ``Flann::Flann``
#   The Flann library, if found.
#
# Result variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
# ``Flann_FOUND``
#   true if the Flann headers and libraries were found
# ``Flann_INCLUDE_DIR``
#   the directory containing the Flann headers
# ``Flann_INCLUDE_DIRS``
#   the directory containing the Flann headers
# ``Flann_LIBRARIES``
#   Flann libraries to be linked
#
# Cache variables
# ^^^^^^^^^^^^^^^
#
# The following cache variables may also be set:
#
# ``Flann_INCLUDE_DIR``
#   the directory containing the Flann headers
# ``Flann_LIBRARY``
#   the path to the Flann library

find_path(Flann_INCLUDE_DIR flann/flann.hpp)

set(Flann_NAMES ${Flann_NAMES} flann flann_cpp_s libflann_cpp_s)
foreach(name ${Flann_NAMES})
  list(APPEND Flann_NAMES_DEBUG "${name}-gd")
endforeach()

if(NOT Flann_LIBRARY)
  find_library(Flann_LIBRARY_RELEASE NAMES ${Flann_NAMES})
  find_library(Flann_LIBRARY_DEBUG NAMES ${Flann_NAMES_DEBUG})
  include(SelectLibraryConfigurations)
  select_library_configurations(Flann)
  mark_as_advanced(Flann_LIBRARY_RELEASE Flann_LIBRARY_DEBUG)
endif()
unset(Flann_NAMES)
unset(Flann_NAMES_DEBUG)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Flann
                                  REQUIRED_VARS Flann_LIBRARY Flann_INCLUDE_DIR
                                  )

if(Flann_FOUND)
  set(Flann_LIBRARIES ${Flann_LIBRARY})
  set(Flann_INCLUDE_DIRS "${Flann_INCLUDE_DIR}")

  if(NOT TARGET Flann::Flann)
    add_library(Flann::Flann UNKNOWN IMPORTED)
    if(Flann_INCLUDE_DIRS)
      set_target_properties(Flann::Flann PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${Flann_INCLUDE_DIRS}")
    endif()
    if(EXISTS "${Flann_LIBRARY}")
      set_target_properties(Flann::Flann PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${Flann_LIBRARY}")
    endif()
    if(EXISTS "${Flann_LIBRARY_RELEASE}")
      set_property(TARGET Flann::Flann APPEND PROPERTY
        IMPORTED_CONFIGURATIONS RELEASE)
      set_target_properties(Flann::Flann PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
        IMPORTED_LOCATION_RELEASE "${Flann_LIBRARY_RELEASE}")
    endif()
    if(EXISTS "${Flann_LIBRARY_DEBUG}")
      set_property(TARGET Flann::Flann APPEND PROPERTY
        IMPORTED_CONFIGURATIONS DEBUG)
      set_target_properties(Flann::Flann PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
        IMPORTED_LOCATION_DEBUG "${Flann_LIBRARY_DEBUG}")
    endif()
  endif()
endif()

mark_as_advanced(Flann_INCLUDE_DIR Flann_LIBRARY)
