# ============================================================================
# FindCeleritas.cmake
# CMake module to find Celeritas GPU-accelerated particle transport library
# ============================================================================
#
# This module finds the Celeritas library and sets the following variables:
#
#   Celeritas_FOUND        - True if Celeritas was found
#   Celeritas_VERSION      - Version string
#   Celeritas_INCLUDE_DIRS - Include directories
#   Celeritas_LIBRARIES    - Libraries to link
#
# The following targets are created:
#   Celeritas::celeritas   - Main Celeritas library
#
# Hints:
#   Celeritas_DIR          - Path to CeleritasConfig.cmake
#   CELERITAS_ROOT         - Root installation directory
#
# ============================================================================

# First try to find Celeritas via its config file
find_package(Celeritas CONFIG QUIET
    HINTS
        ${Celeritas_DIR}
        ${CELERITAS_ROOT}
        $ENV{CELERITAS_ROOT}
        $ENV{CELERITAS_DIR}
        ${CMAKE_PREFIX_PATH}
        /usr/local
        /opt/celeritas
        ${CMAKE_SOURCE_DIR}/../GEANT4_Packages/celeritas-install
    PATH_SUFFIXES
        lib/cmake/Celeritas
        lib64/cmake/Celeritas
        share/cmake/Celeritas
)

if(Celeritas_FOUND)
    message(STATUS "Found Celeritas via config: ${Celeritas_DIR}")
    return()
endif()

# Manual search if config not found
find_path(Celeritas_INCLUDE_DIR
    NAMES celeritas/celeritas_config.h
    HINTS
        ${CELERITAS_ROOT}
        $ENV{CELERITAS_ROOT}
        ${CMAKE_SOURCE_DIR}/../GEANT4_Packages/celeritas-install
    PATH_SUFFIXES include
)

find_library(Celeritas_LIBRARY
    NAMES celeritas
    HINTS
        ${CELERITAS_ROOT}
        $ENV{CELERITAS_ROOT}
        ${CMAKE_SOURCE_DIR}/../GEANT4_Packages/celeritas-install
    PATH_SUFFIXES lib lib64
)

# Extract version if header found
if(Celeritas_INCLUDE_DIR)
    file(STRINGS "${Celeritas_INCLUDE_DIR}/celeritas/celeritas_version.h"
         _version_line REGEX "^#define CELERITAS_VERSION_STRING")
    if(_version_line)
        string(REGEX REPLACE ".*\"([0-9]+\\.[0-9]+\\.[0-9]+)\".*" "\\1"
               Celeritas_VERSION "${_version_line}")
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Celeritas
    REQUIRED_VARS Celeritas_LIBRARY Celeritas_INCLUDE_DIR
    VERSION_VAR Celeritas_VERSION
)

if(Celeritas_FOUND AND NOT TARGET Celeritas::celeritas)
    add_library(Celeritas::celeritas UNKNOWN IMPORTED)
    set_target_properties(Celeritas::celeritas PROPERTIES
        IMPORTED_LOCATION "${Celeritas_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${Celeritas_INCLUDE_DIR}"
    )
endif()

set(Celeritas_INCLUDE_DIRS ${Celeritas_INCLUDE_DIR})
set(Celeritas_LIBRARIES ${Celeritas_LIBRARY})

mark_as_advanced(Celeritas_INCLUDE_DIR Celeritas_LIBRARY)
