# FindGarfield.cmake
# Locate Garfield++ library for gas detector simulation
#
# This module defines:
#   Garfield_FOUND        - True if Garfield++ was found
#   Garfield_INCLUDE_DIRS - Include directories for Garfield++
#   Garfield_LIBRARIES    - Libraries to link against
#   Garfield_VERSION      - Version string
#
# Hints:
#   GARFIELD_HOME         - Root directory of Garfield++ installation
#   GARFIELD_ROOT         - Alternative variable name
#
# Garfield++ requires ROOT, so this module also verifies ROOT is available.

# Check for ROOT first (required by Garfield++)
find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
    set(Garfield_FOUND FALSE)
    if(Garfield_FIND_REQUIRED)
        message(FATAL_ERROR "Garfield++ requires ROOT, but ROOT was not found")
    endif()
    return()
endif()

# Set search paths
set(_garfield_search_paths
    ${GARFIELD_HOME}
    ${GARFIELD_ROOT}
    $ENV{GARFIELD_HOME}
    $ENV{GARFIELD_INSTALL}
    /usr/local
    /usr
    /opt/garfield
)

# Find include directory
find_path(Garfield_INCLUDE_DIR
    NAMES Garfield/ComponentAnalyticField.hh
          Garfield/MediumMagboltz.hh
          Garfield/TrackHeed.hh
    PATHS ${_garfield_search_paths}
    PATH_SUFFIXES include
)

# Find library
find_library(Garfield_LIBRARY
    NAMES Garfield garfield
    PATHS ${_garfield_search_paths}
    PATH_SUFFIXES lib lib64
)

# Check if Heed is available (used for primary ionization)
find_library(Garfield_HEED_LIBRARY
    NAMES Heed heed
    PATHS ${_garfield_search_paths}
    PATH_SUFFIXES lib lib64
)

# Handle version
if(Garfield_INCLUDE_DIR)
    # Try to extract version from header if available
    if(EXISTS "${Garfield_INCLUDE_DIR}/Garfield/Version.hh")
        file(STRINGS "${Garfield_INCLUDE_DIR}/Garfield/Version.hh"
             _version_line REGEX "^#define GARFIELD_VERSION")
        if(_version_line)
            string(REGEX REPLACE "^#define GARFIELD_VERSION \"([^\"]+)\".*" "\\1"
                   Garfield_VERSION "${_version_line}")
        endif()
    endif()
endif()

# Set output variables
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Garfield
    REQUIRED_VARS Garfield_LIBRARY Garfield_INCLUDE_DIR
    VERSION_VAR Garfield_VERSION
)

if(Garfield_FOUND)
    set(Garfield_INCLUDE_DIRS ${Garfield_INCLUDE_DIR})
    set(Garfield_LIBRARIES ${Garfield_LIBRARY})

    # Add Heed if found
    if(Garfield_HEED_LIBRARY)
        list(APPEND Garfield_LIBRARIES ${Garfield_HEED_LIBRARY})
    endif()

    # Add ROOT libraries (Garfield++ depends on ROOT)
    list(APPEND Garfield_LIBRARIES ${ROOT_LIBRARIES})

    # Create imported target
    if(NOT TARGET Garfield::Garfield)
        add_library(Garfield::Garfield UNKNOWN IMPORTED)
        set_target_properties(Garfield::Garfield PROPERTIES
            IMPORTED_LOCATION "${Garfield_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${Garfield_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${ROOT_LIBRARIES}"
        )
    endif()

    message(STATUS "Found Garfield++: ${Garfield_LIBRARY}")
    message(STATUS "  Include: ${Garfield_INCLUDE_DIR}")
    if(Garfield_VERSION)
        message(STATUS "  Version: ${Garfield_VERSION}")
    endif()
endif()

mark_as_advanced(Garfield_INCLUDE_DIR Garfield_LIBRARY Garfield_HEED_LIBRARY)
