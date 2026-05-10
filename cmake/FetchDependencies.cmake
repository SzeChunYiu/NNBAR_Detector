# ============================================================================
# FetchDependencies.cmake
# NNBAR Detector Simulation - Automatic Dependency Management
# ============================================================================
#
# This module automatically handles ALL external dependencies:
# - Bundled dependencies (Arrow, JSON, spdlog, MCPL)
# - GPU acceleration packages (Celeritas, Opticks, Garfield++)
# - Optional packages (Qt for dashboard)
#
# Priority order for each package:
# 1. Bundled in external/ or GEANT4_Packages/
# 2. System-installed packages
# 3. FetchContent automatic download
#
# ============================================================================

include(FetchContent)
include(ExternalProject)

# Directories
set(NNBAR_EXTERNAL_DIR "${PROJECT_SOURCE_DIR}/external")
set(NNBAR_PACKAGES_DIR "${PROJECT_SOURCE_DIR}/../GEANT4_Packages")

# Create packages directory if needed
file(MAKE_DIRECTORY "${NNBAR_PACKAGES_DIR}")

# Status tracking
set(NNBAR_ARROW_SOURCE "Not Found")
set(NNBAR_JSON_SOURCE "Not Found")
set(NNBAR_SPDLOG_SOURCE "Not Found")
set(NNBAR_CELERITAS_SOURCE "Not Found")
set(NNBAR_OPTICKS_SOURCE "Not Found")
set(NNBAR_GARFIELD_SOURCE "Not Found")

# ============================================================================
# Helper: Check CUDA availability for GPU packages
# ============================================================================
macro(check_cuda_availability)
    if(NOT DEFINED CUDA_AVAILABLE)
        find_package(CUDAToolkit QUIET)
        if(CUDAToolkit_FOUND)
            set(CUDA_AVAILABLE TRUE)
            set(CUDA_VERSION "${CUDAToolkit_VERSION}")
        else()
            find_program(NVCC_EXECUTABLE nvcc)
            if(NVCC_EXECUTABLE)
                set(CUDA_AVAILABLE TRUE)
                execute_process(
                    COMMAND ${NVCC_EXECUTABLE} --version
                    OUTPUT_VARIABLE NVCC_OUTPUT
                )
                string(REGEX MATCH "release ([0-9]+\\.[0-9]+)" _ ${NVCC_OUTPUT})
                set(CUDA_VERSION "${CMAKE_MATCH_1}")
            else()
                set(CUDA_AVAILABLE FALSE)
            endif()
        endif()
    endif()
endmacro()

# ============================================================================
# 1. Apache Arrow with Parquet
# ============================================================================
# Check bundled first. The checked-in *-linux bundles are ELF/x86_64 builds and
# must not be auto-selected on macOS.
if(EXISTS "${NNBAR_EXTERNAL_DIR}/arrow-install/lib/cmake/Arrow/ArrowConfig.cmake")
    set(Arrow_DIR "${NNBAR_EXTERNAL_DIR}/arrow-install/lib/cmake/Arrow")
    set(Parquet_DIR "${NNBAR_EXTERNAL_DIR}/arrow-install/lib/cmake/Parquet")
    set(NNBAR_ARROW_SOURCE "Bundled")
elseif(NOT APPLE AND EXISTS "${NNBAR_EXTERNAL_DIR}/arrow-install-linux/lib/cmake/Arrow/ArrowConfig.cmake")
    set(Arrow_DIR "${NNBAR_EXTERNAL_DIR}/arrow-install-linux/lib/cmake/Arrow")
    set(Parquet_DIR "${NNBAR_EXTERNAL_DIR}/arrow-install-linux/lib/cmake/Parquet")
    set(NNBAR_ARROW_SOURCE "Bundled Linux")
endif()

find_package(Arrow QUIET)
find_package(Parquet QUIET)

if(Arrow_FOUND AND Parquet_FOUND)
    if(NNBAR_ARROW_SOURCE MATCHES "^Bundled")
        set(NNBAR_ARROW_SOURCE "${NNBAR_ARROW_SOURCE} (v${Arrow_VERSION})")
    else()
        set(NNBAR_ARROW_SOURCE "System (v${Arrow_VERSION})")
    endif()
else()
    message(STATUS "Arrow not found, will use FetchContent...")
    set(NNBAR_ARROW_SOURCE "FetchContent (building...)")

    FetchContent_Declare(
        arrow
        GIT_REPOSITORY https://github.com/apache/arrow.git
        GIT_TAG apache-arrow-17.0.0
        SOURCE_SUBDIR cpp
    )

    # Minimal Arrow build
    set(ARROW_BUILD_STATIC OFF CACHE BOOL "" FORCE)
    set(ARROW_BUILD_SHARED ON CACHE BOOL "" FORCE)
    set(ARROW_PARQUET ON CACHE BOOL "" FORCE)
    set(ARROW_WITH_SNAPPY ON CACHE BOOL "" FORCE)
    set(ARROW_WITH_ZLIB ON CACHE BOOL "" FORCE)
    set(ARROW_WITH_ZSTD OFF CACHE BOOL "" FORCE)
    set(ARROW_WITH_LZ4 OFF CACHE BOOL "" FORCE)
    set(ARROW_WITH_BROTLI OFF CACHE BOOL "" FORCE)
    set(ARROW_COMPUTE OFF CACHE BOOL "" FORCE)
    set(ARROW_DATASET OFF CACHE BOOL "" FORCE)
    set(ARROW_FILESYSTEM ON CACHE BOOL "" FORCE)
    set(ARROW_JSON OFF CACHE BOOL "" FORCE)
    set(ARROW_CSV OFF CACHE BOOL "" FORCE)
    set(ARROW_FLIGHT OFF CACHE BOOL "" FORCE)
    set(ARROW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(PARQUET_BUILD_EXECUTABLES OFF CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(arrow)
    set(NNBAR_ARROW_SOURCE "FetchContent (v17.0.0)")

    if(NOT TARGET Arrow::arrow_shared)
        add_library(Arrow::arrow_shared ALIAS arrow_shared)
    endif()
    if(NOT TARGET Parquet::parquet_shared)
        add_library(Parquet::parquet_shared ALIAS parquet_shared)
    endif()
endif()

# ============================================================================
# 2. nlohmann/json (Header-only)
# ============================================================================
if(EXISTS "${NNBAR_EXTERNAL_DIR}/json/include/nlohmann/json.hpp")
    add_library(nlohmann_json INTERFACE)
    target_include_directories(nlohmann_json INTERFACE "${NNBAR_EXTERNAL_DIR}/json/include")
    add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)
    set(nlohmann_json_FOUND TRUE)
    set(NNBAR_JSON_SOURCE "Bundled (v3.11.3)")
else()
    find_package(nlohmann_json QUIET)
    if(nlohmann_json_FOUND)
        set(NNBAR_JSON_SOURCE "System (v${nlohmann_json_VERSION})")
    endif()
endif()

if(NOT nlohmann_json_FOUND)
    FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
    )
    set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(nlohmann_json)
    set(NNBAR_JSON_SOURCE "FetchContent (v3.11.3)")
endif()

# ============================================================================
# 3. spdlog (Logging library)
# ============================================================================
if(EXISTS "${NNBAR_EXTERNAL_DIR}/spdlog-install/lib/cmake/spdlog/spdlogConfig.cmake")
    set(spdlog_DIR "${NNBAR_EXTERNAL_DIR}/spdlog-install/lib/cmake/spdlog")
    set(NNBAR_SPDLOG_SOURCE "Bundled")
elseif(NOT APPLE AND EXISTS "${NNBAR_EXTERNAL_DIR}/spdlog-install-linux/lib/cmake/spdlog/spdlogConfig.cmake")
    set(spdlog_DIR "${NNBAR_EXTERNAL_DIR}/spdlog-install-linux/lib/cmake/spdlog")
    set(NNBAR_SPDLOG_SOURCE "Bundled Linux")
endif()

find_package(spdlog QUIET)

if(spdlog_FOUND)
    if(NNBAR_SPDLOG_SOURCE MATCHES "^Bundled")
        set(NNBAR_SPDLOG_SOURCE "${NNBAR_SPDLOG_SOURCE} (v${spdlog_VERSION})")
    else()
        set(NNBAR_SPDLOG_SOURCE "System (v${spdlog_VERSION})")
    endif()
else()
    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.12.0
    )
    set(SPDLOG_BUILD_EXAMPLE OFF CACHE BOOL "" FORCE)
    set(SPDLOG_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(spdlog)
    set(NNBAR_SPDLOG_SOURCE "FetchContent (v1.12.0)")
endif()

# ============================================================================
# 4. Celeritas (GPU EM Physics) - Optional
# ============================================================================
if(WITH_CELERITAS)
    check_cuda_availability()

    if(NOT CUDA_AVAILABLE)
        message(WARNING "Celeritas requires CUDA. Disabling WITH_CELERITAS.")
        set(WITH_CELERITAS OFF CACHE BOOL "" FORCE)
    elseif(Geant4_VERSION VERSION_LESS "11.0")
        message(WARNING "Celeritas requires Geant4 >= 11.0. You have ${Geant4_VERSION}.")
        set(WITH_CELERITAS OFF CACHE BOOL "" FORCE)
    else()
        # Check bundled
        if(EXISTS "${NNBAR_PACKAGES_DIR}/celeritas-install/lib/cmake/Celeritas/CeleritasConfig.cmake")
            set(Celeritas_DIR "${NNBAR_PACKAGES_DIR}/celeritas-install/lib/cmake/Celeritas")
            set(NNBAR_CELERITAS_SOURCE "Bundled")
        endif()

        find_package(Celeritas QUIET)

        if(Celeritas_FOUND)
            if(NNBAR_CELERITAS_SOURCE STREQUAL "Bundled")
                set(NNBAR_CELERITAS_SOURCE "Bundled (v${Celeritas_VERSION})")
            else()
                set(NNBAR_CELERITAS_SOURCE "System (v${Celeritas_VERSION})")
            endif()
        else()
            message(STATUS "Celeritas not found. Install with: ./scripts/install_gpu_packages.sh --celeritas")
            set(WITH_CELERITAS OFF CACHE BOOL "" FORCE)
            set(NNBAR_CELERITAS_SOURCE "Not Found - install with scripts/install_gpu_packages.sh")
        endif()
    endif()
endif()

# ============================================================================
# 5. Opticks (GPU Optical Photons) - Optional
# ============================================================================
if(WITH_OPTICKS)
    check_cuda_availability()

    if(NOT CUDA_AVAILABLE)
        message(WARNING "Opticks requires CUDA. Disabling WITH_OPTICKS.")
        set(WITH_OPTICKS OFF CACHE BOOL "" FORCE)
    else()
        # Check bundled/environment
        if(EXISTS "${NNBAR_PACKAGES_DIR}/opticks/CMakeLists.txt")
            set(OPTICKS_HOME "${NNBAR_PACKAGES_DIR}/opticks")
        elseif(DEFINED ENV{OPTICKS_HOME})
            set(OPTICKS_HOME "$ENV{OPTICKS_HOME}")
        endif()

        if(OPTICKS_HOME)
            list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
            find_package(Opticks QUIET)

            if(Opticks_FOUND)
                set(NNBAR_OPTICKS_SOURCE "Found (${OPTICKS_HOME})")
            else()
                message(STATUS "Opticks found but not built. Run opticks build scripts.")
                set(WITH_OPTICKS OFF CACHE BOOL "" FORCE)
                set(NNBAR_OPTICKS_SOURCE "Not Built - follow Opticks install guide")
            endif()
        else()
            message(STATUS "Opticks not found. Install with: ./scripts/install_gpu_packages.sh --opticks")
            set(WITH_OPTICKS OFF CACHE BOOL "" FORCE)
            set(NNBAR_OPTICKS_SOURCE "Not Found")
        endif()
    endif()
endif()

# ============================================================================
# 6. Garfield++ (TPC Simulation) - Optional
# ============================================================================
if(WITH_GARFIELD)
    # Check bundled
    if(EXISTS "${NNBAR_PACKAGES_DIR}/garfield-install/lib/cmake/Garfield/GarfieldConfig.cmake")
        set(Garfield_DIR "${NNBAR_PACKAGES_DIR}/garfield-install/lib/cmake/Garfield")
        set(NNBAR_GARFIELD_SOURCE "Bundled")
    endif()

    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
    find_package(Garfield QUIET)

    if(Garfield_FOUND)
        if(NNBAR_GARFIELD_SOURCE STREQUAL "Bundled")
            set(NNBAR_GARFIELD_SOURCE "Bundled")
        else()
            set(NNBAR_GARFIELD_SOURCE "System")
        endif()
    else()
        message(STATUS "Garfield++ not found. Install with: ./scripts/install_gpu_packages.sh --garfield")
        set(WITH_GARFIELD OFF CACHE BOOL "" FORCE)
        set(NNBAR_GARFIELD_SOURCE "Not Found")
    endif()
endif()

# ============================================================================
# 7. Qt (Dashboard) - Auto-detect
# ============================================================================
# Try to find Qt even if not explicitly requested - enable dashboard if found
if(NOT DEFINED Qt_FOUND_CHECKED)
    set(Qt_FOUND_CHECKED TRUE)

    # Try Qt6 first
    find_package(Qt6 COMPONENTS Core Widgets QUIET)
    if(Qt6_FOUND)
        set(QT_AVAILABLE TRUE)
        set(QT_VERSION "6")
        set(QT_LIBRARIES Qt6::Core Qt6::Widgets)
    else()
        # Fall back to Qt5
        find_package(Qt5 COMPONENTS Core Widgets QUIET)
        if(Qt5_FOUND)
            set(QT_AVAILABLE TRUE)
            set(QT_VERSION "5")
            set(QT_LIBRARIES Qt5::Core Qt5::Widgets)
        else()
            set(QT_AVAILABLE FALSE)
        endif()
    endif()
endif()

# ============================================================================
# Export status variables
# ============================================================================
get_directory_property(NNBAR_HAS_PARENT_DIRECTORY PARENT_DIRECTORY)
if(NNBAR_HAS_PARENT_DIRECTORY)
    set(NNBAR_ARROW_SOURCE "${NNBAR_ARROW_SOURCE}" PARENT_SCOPE)
    set(NNBAR_JSON_SOURCE "${NNBAR_JSON_SOURCE}" PARENT_SCOPE)
    set(NNBAR_SPDLOG_SOURCE "${NNBAR_SPDLOG_SOURCE}" PARENT_SCOPE)
    set(NNBAR_CELERITAS_SOURCE "${NNBAR_CELERITAS_SOURCE}" PARENT_SCOPE)
    set(NNBAR_OPTICKS_SOURCE "${NNBAR_OPTICKS_SOURCE}" PARENT_SCOPE)
    set(NNBAR_GARFIELD_SOURCE "${NNBAR_GARFIELD_SOURCE}" PARENT_SCOPE)
    set(QT_AVAILABLE "${QT_AVAILABLE}" PARENT_SCOPE)
    set(QT_VERSION "${QT_VERSION}" PARENT_SCOPE)
    set(QT_LIBRARIES "${QT_LIBRARIES}" PARENT_SCOPE)
    set(CUDA_AVAILABLE "${CUDA_AVAILABLE}" PARENT_SCOPE)
    set(CUDA_VERSION "${CUDA_VERSION}" PARENT_SCOPE)
endif()
