# ============================================================================
# FindOpticks.cmake
# CMake module to find Opticks GPU optical photon propagation library
# ============================================================================
#
# This module finds the Opticks library and sets the following variables:
#
#   Opticks_FOUND        - True if Opticks was found (full installation)
#   Opticks_PARTIAL      - True if only partial Opticks found (CSGOptiX only)
#   Opticks_INCLUDE_DIRS - Include directories
#   Opticks_LIBRARIES    - Libraries to link
#   OPTICKS_HOME         - Root directory
#
# The following targets are created:
#   Opticks::Opticks     - Main Opticks library (if full install found)
#
# Hints:
#   OPTICKS_HOME         - Root installation directory
#   OPTICKS_PREFIX       - Alternative to OPTICKS_HOME
#
# Note: Full Opticks requires G4CX library which needs a specially patched
# Geant4 with G4Tree. Without this, only CSGOptiX ray tracing is available.
#
# ============================================================================

# Look for OPTICKS_HOME
if(NOT OPTICKS_HOME)
    if(DEFINED ENV{OPTICKS_HOME})
        set(OPTICKS_HOME $ENV{OPTICKS_HOME})
    elseif(DEFINED ENV{OPTICKS_PREFIX})
        set(OPTICKS_HOME $ENV{OPTICKS_PREFIX})
    else()
        # Try common locations
        foreach(_path
            "${CMAKE_SOURCE_DIR}/../GEANT4_Packages/install/opticks"
            "${CMAKE_SOURCE_DIR}/../GEANT4_Packages/opticks"
            "/opt/opticks"
            "$ENV{HOME}/opticks"
        )
            if(EXISTS "${_path}/lib/libCSG.so" OR EXISTS "${_path}/lib/libCSGOptiX.so")
                set(OPTICKS_HOME "${_path}")
                break()
            endif()
        endforeach()
    endif()
endif()

if(NOT OPTICKS_HOME)
    set(Opticks_FOUND FALSE)
    set(Opticks_PARTIAL FALSE)
    if(Opticks_FIND_REQUIRED)
        message(FATAL_ERROR "Opticks not found. Set OPTICKS_HOME environment variable.")
    endif()
    return()
endif()

message(STATUS "Looking for Opticks in: ${OPTICKS_HOME}")

# Find core Opticks libraries (always needed)
find_library(Opticks_CSG_LIBRARY
    NAMES CSG
    HINTS ${OPTICKS_HOME}/lib ${OPTICKS_HOME}/lib64
)

find_library(Opticks_SYSRAP_LIBRARY
    NAMES SysRap
    HINTS ${OPTICKS_HOME}/lib ${OPTICKS_HOME}/lib64
)

find_library(Opticks_CSGOPTIX_LIBRARY
    NAMES CSGOptiX
    HINTS ${OPTICKS_HOME}/lib ${OPTICKS_HOME}/lib64
)

# Find G4CX (Geant4-Opticks interface) - requires special Geant4 build
find_library(Opticks_G4CX_LIBRARY
    NAMES G4CX
    HINTS ${OPTICKS_HOME}/lib ${OPTICKS_HOME}/lib64
)

find_library(Opticks_U4_LIBRARY
    NAMES U4
    HINTS ${OPTICKS_HOME}/lib ${OPTICKS_HOME}/lib64
)

find_library(Opticks_GDXML_LIBRARY
    NAMES GDXML
    HINTS ${OPTICKS_HOME}/lib ${OPTICKS_HOME}/lib64
)

# Find include directories
find_path(Opticks_CSG_INCLUDE_DIR
    NAMES CSGFoundry.h
    HINTS ${OPTICKS_HOME}/include/CSG ${OPTICKS_HOME}/include
)

find_path(Opticks_SYSRAP_INCLUDE_DIR
    NAMES SEvt.hh
    HINTS ${OPTICKS_HOME}/include/SysRap ${OPTICKS_HOME}/include
)

# Check what we found
set(Opticks_PARTIAL FALSE)
set(Opticks_FOUND FALSE)

if(Opticks_CSG_LIBRARY AND Opticks_SYSRAP_LIBRARY)
    set(Opticks_PARTIAL TRUE)
    message(STATUS "Found Opticks core libraries (CSG, SysRap)")

    if(Opticks_CSGOPTIX_LIBRARY)
        message(STATUS "Found CSGOptiX (GPU ray tracing)")
    endif()

    if(Opticks_G4CX_LIBRARY AND Opticks_U4_LIBRARY)
        set(Opticks_FOUND TRUE)
        message(STATUS "Found G4CX and U4 (full Geant4-Opticks integration)")
    else()
        message(STATUS "G4CX/U4 not found - full Opticks requires specially patched Geant4 with G4Tree")
        message(STATUS "  See: https://bitbucket.org/simoncblyth/opticks/wiki/Home")
        message(STATUS "  GPU optical photon simulation will not be available")
    endif()
else()
    message(STATUS "Opticks core libraries not found")
    if(Opticks_FIND_REQUIRED)
        message(FATAL_ERROR "Required Opticks libraries not found at ${OPTICKS_HOME}")
    endif()
    return()
endif()

# Set up include directories
set(Opticks_INCLUDE_DIRS "")
if(Opticks_CSG_INCLUDE_DIR)
    list(APPEND Opticks_INCLUDE_DIRS ${Opticks_CSG_INCLUDE_DIR})
endif()
if(Opticks_SYSRAP_INCLUDE_DIR)
    list(APPEND Opticks_INCLUDE_DIRS ${Opticks_SYSRAP_INCLUDE_DIR})
endif()
list(APPEND Opticks_INCLUDE_DIRS
    ${OPTICKS_HOME}/include
    ${OPTICKS_HOME}/include/CSG
    ${OPTICKS_HOME}/include/SysRap
    ${OPTICKS_HOME}/include/QUDARap
    ${OPTICKS_HOME}/include/CSGOptiX
)

# Set up libraries based on what's available
set(Opticks_LIBRARIES "")
if(Opticks_FOUND)
    # Full installation with Geant4 integration
    list(APPEND Opticks_LIBRARIES
        ${Opticks_G4CX_LIBRARY}
        ${Opticks_U4_LIBRARY}
    )
    if(Opticks_GDXML_LIBRARY)
        list(APPEND Opticks_LIBRARIES ${Opticks_GDXML_LIBRARY})
    endif()
endif()

# Always include available core libraries
if(Opticks_CSGOPTIX_LIBRARY)
    list(APPEND Opticks_LIBRARIES ${Opticks_CSGOPTIX_LIBRARY})
endif()
if(Opticks_CSG_LIBRARY)
    list(APPEND Opticks_LIBRARIES ${Opticks_CSG_LIBRARY})
endif()
if(Opticks_SYSRAP_LIBRARY)
    list(APPEND Opticks_LIBRARIES ${Opticks_SYSRAP_LIBRARY})
endif()

# Create imported target only if full installation found
if(Opticks_FOUND)
    if(NOT TARGET Opticks::Opticks)
        add_library(Opticks::Opticks INTERFACE IMPORTED)
        set_target_properties(Opticks::Opticks PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${Opticks_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${Opticks_LIBRARIES}"
        )
    endif()
    message(STATUS "Opticks: Full installation found - GPU optical photons enabled")
else()
    message(STATUS "Opticks: Partial installation - GPU optical photons NOT available")
endif()

# Report status
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Opticks
    REQUIRED_VARS OPTICKS_HOME
    HANDLE_COMPONENTS
)

mark_as_advanced(
    Opticks_CSG_LIBRARY
    Opticks_SYSRAP_LIBRARY
    Opticks_CSGOPTIX_LIBRARY
    Opticks_G4CX_LIBRARY
    Opticks_U4_LIBRARY
    Opticks_GDXML_LIBRARY
    Opticks_CSG_INCLUDE_DIR
    Opticks_SYSRAP_INCLUDE_DIR
)
