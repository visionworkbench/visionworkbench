# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#.rst:
# FindProj
# --------
#
# Find the Proj library (libProj).
#
# Imported targets
# ^^^^^^^^^^^^^^^^
#
# This module defines the following :prop_tgt:`IMPORTED` targets:
#
# ``Proj::Proj``
#   The Proj library, if found.
#
# Result variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
# ``Proj_FOUND``
#   true if the Proj headers and libraries were found
# ``Proj_INCLUDE_DIR``
#   the directory containing the Proj headers
# ``Proj_INCLUDE_DIRS``
#   the directory containing the Proj headers
# ``Proj_LIBRARIES``
#   Proj libraries to be linked
#
# Cache variables
# ^^^^^^^^^^^^^^^
#
# The following cache variables may also be set:
#
# ``Proj_INCLUDE_DIR``
#   the directory containing the Proj headers
# ``Proj_LIBRARY``
#   the path to the Proj library

find_path(Proj_INCLUDE_DIR proj_api.h)

set(Proj_NAMES ${Proj_NAMES} proj)
foreach(name ${Proj_NAMES})
  list(APPEND Proj_NAMES_DEBUG "${name}d")
endforeach()

if(NOT Proj_LIBRARY)
  find_library(Proj_LIBRARY_RELEASE NAMES ${Proj_NAMES})
  find_library(Proj_LIBRARY_DEBUG NAMES ${Proj_NAMES_DEBUG})
  include(SelectLibraryConfigurations)
  select_library_configurations(Proj)
  mark_as_advanced(Proj_LIBRARY_RELEASE Proj_LIBRARY_DEBUG)
endif()
unset(Proj_NAMES)
unset(Proj_NAMES_DEBUG)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Proj
                                  REQUIRED_VARS Proj_LIBRARY Proj_INCLUDE_DIR
                                  )

if(Proj_FOUND)
  set(Proj_LIBRARIES ${Proj_LIBRARY})
  set(Proj_INCLUDE_DIRS "${Proj_INCLUDE_DIR}")

  if(NOT TARGET Proj::Proj)
    add_library(Proj::Proj UNKNOWN IMPORTED)
    if(Proj_INCLUDE_DIRS)
      set_target_properties(Proj::Proj PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${Proj_INCLUDE_DIRS}")
    endif()
    if(EXISTS "${Proj_LIBRARY}")
      set_target_properties(Proj::Proj PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${Proj_LIBRARY}")
    endif()
    if(EXISTS "${Proj_LIBRARY_RELEASE}")
      set_property(TARGET Proj::Proj APPEND PROPERTY
        IMPORTED_CONFIGURATIONS RELEASE)
      set_target_properties(Proj::Proj PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
        IMPORTED_LOCATION_RELEASE "${Proj_LIBRARY_RELEASE}")
    endif()
    if(EXISTS "${Proj_LIBRARY_DEBUG}")
      set_property(TARGET Proj::Proj APPEND PROPERTY
        IMPORTED_CONFIGURATIONS DEBUG)
      set_target_properties(Proj::Proj PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
        IMPORTED_LOCATION_DEBUG "${Proj_LIBRARY_DEBUG}")
    endif()
  endif()
endif()

mark_as_advanced(Proj_INCLUDE_DIR Proj_LIBRARY)
