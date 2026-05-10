# ============================================================================
# NNBAR Detector Simulation - CMake Build Configuration
# ============================================================================
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(nnbar-detector-simulation
    VERSION 1.0.0
    LANGUAGES CXX C
    DESCRIPTION "NNBAR/HIBEAM Detector Monte Carlo Simulation"
)

# ============================================================================
# C++ Standard Configuration
# ============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile commands for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add local cmake modules (includes FindCUDAToolkit for CMake < 3.17 compatibility)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

# ============================================================================
# Compiler Flags
# ============================================================================
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()
if(MSVC)
    add_compile_options(/W4)
endif()

# ============================================================================
# RPATH Configuration - Embed library paths in executable
# ============================================================================
# Don't skip RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# Use the install RPATH already during build (so it runs from build dir)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# Set the RPATH to use when installing. Additional dependency paths are appended
# after package discovery, once CMake knows where Geant4 and imported libraries
# came from.
set(CMAKE_INSTALL_RPATH "")
if(EXISTS "${CMAKE_SOURCE_DIR}/external/arrow-install/lib")
    list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/external/arrow-install/lib")
endif()
if(NOT APPLE)
    foreach(path
        "${CMAKE_SOURCE_DIR}/external/arrow-install-linux/lib"
        "/usr/local/cuda-12.4/lib64"
        "/usr/lib/x86_64-linux-gnu"
    )
        if(EXISTS "${path}")
            list(APPEND CMAKE_INSTALL_RPATH "${path}")
        endif()
    endforeach()
endif()

# Add the automatically determined parts of the RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Use RPATH instead of RUNPATH (RPATH has higher priority than LD_LIBRARY_PATH)
# --disable-new-dtags is Linux-only (ELF); skip on macOS (Mach-O uses LC_RPATH)
if(NOT APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags")
endif()

# ============================================================================
# Build Options
# ============================================================================
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
option(MCPL_BUILD "Enable MCPL particle source (0=ParticleGun, 1=MCPL)" OFF)
option(TARGET_BUILD "Build with carbon target foil" ON)
option(DEBUG_VERBOSE "Enable debug output" OFF)
# Physics options
option(WITH_SCINTILLATION "Enable optical photon generation (slower simulation)" OFF)
option(WITH_GARFIELD "Enable Garfield++ TPC simulation" OFF)
option(WITH_GARFIELD_GPU "Enable GPU-accelerated TPC drift simulation (custom CUDA)" OFF)

# GPU acceleration - explicit opt-in; default builds must stay portable.
option(WITH_CELERITAS "Enable Celeritas GPU-accelerated EM physics" OFF)
option(WITH_OPTICKS "Enable Opticks GPU optical photon propagation" OFF)

# Dashboard - AUTO-ENABLED when Qt is found
option(WITH_DASHBOARD "Enable Qt monitoring dashboard (auto-enabled if Qt found)" OFF)
option(WITH_CRY "Enable CRY cosmic ray generator" OFF)

# Convert options to 0/1 for config.h
if(MCPL_BUILD)
    set(MCPL_BUILD_VALUE 1)
else()
    set(MCPL_BUILD_VALUE 0)
endif()

if(TARGET_BUILD)
    set(TARGET_BUILD_VALUE 1)
else()
    set(TARGET_BUILD_VALUE 0)
endif()

if(DEBUG_VERBOSE)
    set(DEBUG_VERBOSE_VALUE 1)
else()
    set(DEBUG_VERBOSE_VALUE 0)
endif()

if(WITH_SCINTILLATION)
    set(WITH_SCINTILLATION_VALUE 1)
else()
    set(WITH_SCINTILLATION_VALUE 0)
endif()

if(WITH_GARFIELD)
    set(WITH_GARFIELD_VALUE 1)
else()
    set(WITH_GARFIELD_VALUE 0)
endif()

if(WITH_GARFIELD_GPU)
    set(WITH_GARFIELD_GPU_VALUE 1)
else()
    set(WITH_GARFIELD_GPU_VALUE 0)
endif()

if(WITH_CELERITAS)
    set(WITH_CELERITAS_VALUE 1)
else()
    set(WITH_CELERITAS_VALUE 0)
endif()

if(WITH_OPTICKS)
    set(WITH_OPTICKS_VALUE 1)
else()
    set(WITH_OPTICKS_VALUE 0)
endif()

if(WITH_DASHBOARD)
    set(WITH_DASHBOARD_VALUE 1)
else()
    set(WITH_DASHBOARD_VALUE 0)
endif()

if(WITH_CRY)
    set(WITH_CRY_VALUE 1)
else()
    set(WITH_CRY_VALUE 0)
endif()

# ============================================================================
# ANSI Color Codes for Professional Output
# ============================================================================
if(NOT WIN32)
    string(ASCII 27 ESC)
    set(RESET   "${ESC}[0m")
    set(BOLD    "${ESC}[1m")
    set(RED     "${ESC}[31m")
    set(GREEN   "${ESC}[32m")
    set(YELLOW  "${ESC}[33m")
    set(BLUE    "${ESC}[34m")
    set(MAGENTA "${ESC}[35m")
    set(CYAN    "${ESC}[36m")
    set(WHITE   "${ESC}[37m")
    set(CHECK   "${GREEN}[OK]${RESET}")
    set(CROSS   "${RED}[--]${RESET}")
    set(WARN    "${YELLOW}[!!]${RESET}")
else()
    set(RESET   "")
    set(BOLD    "")
    set(RED     "")
    set(GREEN   "")
    set(YELLOW  "")
    set(BLUE    "")
    set(MAGENTA "")
    set(CYAN    "")
    set(WHITE   "")
    set(CHECK   "[OK]")
    set(CROSS   "[--]")
    set(WARN    "[!!]")
endif()

# ============================================================================
# Initial Banner
# ============================================================================
message("")
message("${CYAN}╔══════════════════════════════════════════════════════════════════════════╗${RESET}")
message("${CYAN}║${RESET}${BOLD}          NNBAR/HIBEAM Detector Simulation - Build Configuration          ${RESET}${CYAN}║${RESET}")
message("${CYAN}║${RESET}          Geant4-based Monte Carlo for Neutron-Antineutron Oscillation     ${CYAN}║${RESET}")
message("${CYAN}╚══════════════════════════════════════════════════════════════════════════╝${RESET}")
message("")

# Note: config.h is generated AFTER all find_package calls (see below)
# to ensure correct feature detection values

# ============================================================================
# Find Geant4 (Required)
# ============================================================================
message("${BLUE}${BOLD}[1/4] Checking Required Dependencies${RESET}")
message("${BLUE}────────────────────────────────────────────────────────────────────────────${RESET}")

if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else()
    find_package(Geant4 REQUIRED)
endif()
include(${Geant4_USE_FILE})
message("  ${CHECK} Geant4 ${GREEN}v${Geant4_VERSION}${RESET} - Monte Carlo simulation toolkit")

set(NNBAR_GEANT4_INSTALL "")
if(Geant4_DIR)
    get_filename_component(NNBAR_GEANT4_INSTALL "${Geant4_DIR}/../../.." ABSOLUTE)
    if(EXISTS "${NNBAR_GEANT4_INSTALL}/lib")
        list(APPEND CMAKE_INSTALL_RPATH "${NNBAR_GEANT4_INSTALL}/lib")
    endif()
endif()

# ============================================================================
# Optional: Garfield++ for TPC simulation
# ============================================================================
message("")
message("${BLUE}${BOLD}[2/4] Checking Optional Physics Libraries${RESET}")
message("${BLUE}────────────────────────────────────────────────────────────────────────────${RESET}")

set(NNBAR_GARFIELD_STATUS "Disabled")
set(NNBAR_CELERITAS_STATUS "Disabled")
set(NNBAR_OPTICKS_STATUS "Disabled")
set(NNBAR_QT_STATUS "Disabled")
set(Celeritas_FOUND OFF)
set(Opticks_FOUND OFF)

if(WITH_GARFIELD)
    list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
    find_package(Garfield QUIET)
    if(Garfield_FOUND)
        message("  ${CHECK} Garfield++ ${GREEN}Found${RESET} - Realistic TPC electron drift simulation")
        set(NNBAR_GARFIELD_STATUS "Enabled")
        add_definitions(-DWITH_GARFIELD)
    else()
        message("  ${WARN} Garfield++ ${YELLOW}Not Found${RESET} - Set Garfield_DIR or GARFIELD_HOME")
        set(WITH_GARFIELD OFF)
        set(WITH_GARFIELD_VALUE 0)
    endif()
else()
    message("  ${CROSS} Garfield++ ${WHITE}Disabled${RESET} - Use -DWITH_GARFIELD=ON to enable")
endif()

# ============================================================================
# Optional: GarfieldGPU - Custom GPU-accelerated TPC drift (CUDA/OpenMP)
# ============================================================================
set(NNBAR_GARFIELD_GPU_STATUS "Disabled")
set(GARFIELD_GPU_USE_CUDA OFF)

if(WITH_GARFIELD_GPU)
    # Check for CUDA
    include(CheckLanguage)
    check_language(CUDA)

    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        # CMake 3.16 doesn't support CUDA_STANDARD 17, use nvcc flags directly
        if(CMAKE_VERSION VERSION_LESS "3.18")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++17")
        else()
            set(CMAKE_CUDA_STANDARD 17)
            set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        endif()

        # Set CUDA hints for CUDAToolkit detection
        if(NOT CUDAToolkit_ROOT AND EXISTS "/usr/local/cuda-12.4")
            set(CUDAToolkit_ROOT "/usr/local/cuda-12.4")
        elseif(NOT CUDAToolkit_ROOT AND EXISTS "/usr/local/cuda")
            set(CUDAToolkit_ROOT "/usr/local/cuda")
        endif()

        # Find CUDAToolkit for CUDA::cudart and CUDA::curand
        find_package(CUDAToolkit QUIET)
        if(CUDAToolkit_FOUND)
            set(GARFIELD_GPU_USE_CUDA ON)
            message("  ${CHECK} GarfieldGPU ${GREEN}CUDA ${CUDAToolkit_VERSION}${RESET} - GPU electron drift (100-1000x speedup)")
            set(NNBAR_GARFIELD_GPU_STATUS "Enabled (CUDA)")
            add_definitions(-DWITH_GARFIELD_GPU -DGARFIELD_GPU_CUDA)
        else()
            # Fall back to OpenMP if CUDAToolkit not found
            find_package(OpenMP QUIET)
            if(OpenMP_CXX_FOUND)
                message("  ${WARN} GarfieldGPU ${YELLOW}CUDAToolkit not found${RESET} - Using OpenMP CPU fallback")
                set(NNBAR_GARFIELD_GPU_STATUS "Enabled (OpenMP)")
                add_definitions(-DWITH_GARFIELD_GPU)
            else()
                message("  ${WARN} GarfieldGPU ${YELLOW}Neither CUDA nor OpenMP found${RESET} - Single-threaded fallback")
                set(NNBAR_GARFIELD_GPU_STATUS "Enabled (Single-threaded)")
                add_definitions(-DWITH_GARFIELD_GPU)
            endif()
        endif()
    else()
        # Fall back to OpenMP CPU parallelization
        find_package(OpenMP QUIET)
        if(OpenMP_CXX_FOUND)
            message("  ${WARN} GarfieldGPU ${YELLOW}CUDA not found${RESET} - Using OpenMP CPU fallback")
            set(NNBAR_GARFIELD_GPU_STATUS "Enabled (OpenMP)")
            add_definitions(-DWITH_GARFIELD_GPU)
        else()
            message("  ${WARN} GarfieldGPU ${YELLOW}Neither CUDA nor OpenMP found${RESET} - Single-threaded fallback")
            set(NNBAR_GARFIELD_GPU_STATUS "Enabled (Single-threaded)")
            add_definitions(-DWITH_GARFIELD_GPU)
        endif()
    endif()
else()
    message("  ${CROSS} GarfieldGPU ${WHITE}Disabled${RESET} - Use -DWITH_GARFIELD_GPU=ON to enable")
endif()

# ============================================================================
# Optional: Celeritas for GPU-accelerated EM physics
# ============================================================================
if(WITH_CELERITAS)
    # Check Geant4 version - Celeritas requires 11.0+
    if(Geant4_VERSION VERSION_LESS "11.0")
        message("  ${WARN} Celeritas ${YELLOW}Requires Geant4 ≥ 11.0${RESET} (you have ${Geant4_VERSION})")
        message("         Run: ${CYAN}./scripts/upgrade_geant4.sh${RESET}")
        set(WITH_CELERITAS OFF)
        set(WITH_CELERITAS_VALUE 0)
    else()
        find_package(Celeritas QUIET)
        if(Celeritas_FOUND)
            message("  ${CHECK} Celeritas ${GREEN}v${Celeritas_VERSION}${RESET} - GPU-accelerated EM physics (10-100x speedup)")
            set(NNBAR_CELERITAS_STATUS "Enabled")
            add_definitions(-DWITH_CELERITAS)
        else()
            message("  ${WARN} Celeritas ${YELLOW}Not Found${RESET} - Run: ./scripts/install_celeritas.sh")
            set(WITH_CELERITAS OFF)
            set(WITH_CELERITAS_VALUE 0)
        endif()
    endif()
else()
    message("  ${CROSS} Celeritas ${WHITE}Disabled${RESET} - Use -DWITH_CELERITAS=ON to enable")
endif()

# ============================================================================
# Optional: Opticks for GPU optical photon propagation
# ============================================================================
if(WITH_OPTICKS)
    list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
    find_package(Opticks QUIET)
    if(Opticks_FOUND)
        message("  ${CHECK} Opticks ${GREEN}Found${RESET} - GPU optical photon propagation (50-200x speedup)")
        set(NNBAR_OPTICKS_STATUS "Enabled (GPU)")
        set(OPTICKS_FULL_INSTALL ON)
        add_definitions(-DWITH_OPTICKS)
        # Disable assertions in Opticks headers (has strict binary tree assertions)
        add_definitions(-DNDEBUG)

        # Add Opticks include directories
        if(OPTICKS_HOME)
            set(OPTICKS_INCLUDE_DIRS
                ${OPTICKS_HOME}/include/G4CX
                ${OPTICKS_HOME}/include/U4
                ${OPTICKS_HOME}/include/CSGOptiX
                ${OPTICKS_HOME}/include/CSG
                ${OPTICKS_HOME}/include/SysRap
                ${OPTICKS_HOME}/include/GDXML
                ${OPTICKS_HOME}/include/QUDARap
                ${OPTICKS_HOME}/externals/plog/include
                ${OPTICKS_HOME}/externals/glm/glm-0.9.9.5/glm
            )
        endif()
    else()
        message("  ${WARN} Opticks ${YELLOW}Not Found${RESET} - using CPU optical simulation")
        set(WITH_OPTICKS OFF CACHE BOOL "" FORCE)
        set(WITH_OPTICKS_VALUE 0)
        set(NNBAR_OPTICKS_STATUS "Disabled")
        set(OPTICKS_FULL_INSTALL OFF)
    endif()
else()
    message("  ${CROSS} Opticks ${WHITE}Disabled${RESET} - Use -DWITH_OPTICKS=ON to enable")
endif()

# ============================================================================
# Optional: Qt Dashboard for monitoring
# ============================================================================
if(WITH_DASHBOARD)
    # Try Qt6 first, fall back to Qt5
    find_package(Qt6 COMPONENTS Core Widgets QUIET)
    if(Qt6_FOUND)
        set(QT_LIBRARIES Qt6::Core Qt6::Widgets)
        set(CMAKE_AUTOMOC ON)
        set(CMAKE_AUTORCC ON)
        set(CMAKE_AUTOUIC ON)
        # Pass WITH_DASHBOARD to MOC so it can see the Q_OBJECT macro
        set(CMAKE_AUTOMOC_MOC_OPTIONS "-DWITH_DASHBOARD=1")
        add_compile_definitions(WITH_DASHBOARD=1)
        message("  ${CHECK} Qt Dashboard ${GREEN}Qt6${RESET} - Real-time monitoring interface")
        set(NNBAR_QT_STATUS "Enabled (Qt6)")
    else()
        find_package(Qt5 COMPONENTS Core Widgets OpenGL QUIET)
        if(Qt5_FOUND)
            set(QT_LIBRARIES Qt5::Core Qt5::Widgets Qt5::OpenGL)
            # Phoenix display requires Qt5WebEngineWidgets which may not be installed
            # To enable, install: sudo apt install qtwebengine5-dev
            add_compile_definitions(WITH_PHOENIX=0)
            message("  ${CROSS} Qt WebEngine ${WHITE}Skipped${RESET} - Install qtwebengine5-dev to enable Phoenix")
            set(CMAKE_AUTOMOC ON)
            set(CMAKE_AUTORCC ON)
            set(CMAKE_AUTOUIC ON)
            # Pass WITH_DASHBOARD to MOC so it can see the Q_OBJECT macro
            set(CMAKE_AUTOMOC_MOC_OPTIONS "-DWITH_DASHBOARD=1")
            add_compile_definitions(WITH_DASHBOARD=1)
            message("  ${CHECK} Qt Dashboard ${GREEN}Qt5${RESET} - Real-time monitoring interface")
            set(NNBAR_QT_STATUS "Enabled (Qt5)")
        else()
            message("  ${WARN} Qt Dashboard ${YELLOW}Qt Not Found${RESET} - Install Qt5 or Qt6")
            set(WITH_DASHBOARD OFF)
            set(WITH_DASHBOARD_VALUE 0)
        endif()
    endif()
else()
    message("  ${CROSS} Qt Dashboard ${WHITE}Disabled${RESET} - Use -DWITH_DASHBOARD=ON to enable")
endif()


# ============================================================================
# Optional: CRY cosmic ray generator
# ============================================================================
set(NNBAR_CRY_STATUS "Disabled")
if(WITH_CRY)
    set(CRY_DIR "" CACHE PATH "CRY v1.7 installation prefix")
    if(NOT CRY_DIR AND DEFINED ENV{CRY_DIR})
        set(CRY_DIR "$ENV{CRY_DIR}" CACHE PATH "CRY v1.7 installation prefix" FORCE)
    endif()
    find_path(CRY_INCLUDE_DIR CRYSetup.h HINTS ${CRY_DIR}/src $ENV{CRY_DIR}/src REQUIRED)
    find_library(CRY_LIBRARY NAMES CRY HINTS ${CRY_DIR}/lib $ENV{CRY_DIR}/lib REQUIRED)
    add_compile_definitions(WITH_CRY=1)
    include_directories(${CRY_INCLUDE_DIR})
    set(NNBAR_CRY_STATUS "Enabled")
    message("  ${CHECK} CRY ${GREEN}Found${RESET} - Cosmic ray shower generator")
else()
    message("  ${CROSS} CRY ${WHITE}Disabled${RESET} - Use -DWITH_CRY=ON to enable")
endif()

# ============================================================================
# Generate config.h - AFTER all find_package calls for correct values
# ============================================================================
configure_file(
    ${PROJECT_SOURCE_DIR}/config.h.in
    ${CMAKE_BINARY_DIR}/config.h
    @ONLY
)

# ============================================================================
# Fetch External Dependencies (Arrow/Parquet, nlohmann_json, spdlog)
# ============================================================================
message("")
message("${BLUE}${BOLD}[3/4] Checking Data Processing Libraries${RESET}")
message("${BLUE}────────────────────────────────────────────────────────────────────────────${RESET}")

include(${PROJECT_SOURCE_DIR}/cmake/FetchDependencies.cmake)

message("  ${CHECK} Arrow/Parquet: ${GREEN}${NNBAR_ARROW_SOURCE}${RESET} - Columnar data output")
message("  ${CHECK} nlohmann_json: ${GREEN}${NNBAR_JSON_SOURCE}${RESET} - JSON configuration")
message("  ${CHECK} spdlog:        ${GREEN}${NNBAR_SPDLOG_SOURCE}${RESET} - Fast logging")

# ============================================================================
# MCPL Library (Monte Carlo Particle List) - Always Bundled
# ============================================================================
message("  ${CHECK} MCPL:          ${GREEN}Bundled${RESET} - Monte Carlo Particle List format")
message("")
add_library(mcpl STATIC
    ${PROJECT_SOURCE_DIR}/external/mcpl/mcpl.c
)
target_include_directories(mcpl PUBLIC
    ${PROJECT_SOURCE_DIR}/external/mcpl
)
set_target_properties(mcpl PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    C_STANDARD 99
)

# ============================================================================
# Parquet Writer Library
# ============================================================================
set(PARQUET_WRITER_DIR ${PROJECT_SOURCE_DIR}/external/parquet-writer/src/cpp)

add_library(parquet_writer STATIC
    ${PARQUET_WRITER_DIR}/parquet_writer.cpp
    ${PARQUET_WRITER_DIR}/parquet_writer_helpers.cpp
    ${PARQUET_WRITER_DIR}/parquet_writer_fill_types.cpp
    ${PARQUET_WRITER_DIR}/parquet_writer_exceptions.cpp
    ${PARQUET_WRITER_DIR}/logging.cpp
)

target_include_directories(parquet_writer PUBLIC
    ${PARQUET_WRITER_DIR}
)

# Suppress warnings from external parquet-writer code
target_compile_options(parquet_writer PRIVATE -Wno-pedantic -Wno-unused-parameter -Wno-deprecated-declarations)

# Link parquet_writer dependencies
if(TARGET Arrow::arrow_shared AND TARGET Parquet::parquet_shared)
    target_link_libraries(parquet_writer PUBLIC
        Arrow::arrow_shared
        Parquet::parquet_shared
        nlohmann_json::nlohmann_json
        spdlog::spdlog
    )
elseif(TARGET arrow_shared AND TARGET parquet_shared)
    target_link_libraries(parquet_writer PUBLIC
        arrow_shared
        parquet_shared
        nlohmann_json::nlohmann_json
        spdlog::spdlog
    )
else()
    message(FATAL_ERROR "Arrow/Parquet targets not found")
endif()

# ============================================================================
# Main Executable Sources
# ============================================================================
# Source files organized in subdirectories:
#   src/core/       - Core Geant4 actions and initialization
#   src/detector/   - Detector geometry construction
#   src/sensitive/  - Sensitive detector implementations
#   src/hits/       - Hit and run accumulation
#   src/generator/  - Primary particle generators (MCPL, ParticleGun)
#   src/output/     - Output management (Parquet)
#   src/physics/    - Physics lists and optional integrations
#   src/gui/        - Qt dashboard interface
#   src/util/       - Utility classes

file(GLOB SOURCES
    ${PROJECT_SOURCE_DIR}/src/core/*.cc
    ${PROJECT_SOURCE_DIR}/src/detector/*.cc
    ${PROJECT_SOURCE_DIR}/src/sensitive/*.cc
    ${PROJECT_SOURCE_DIR}/src/hits/*.cc
    ${PROJECT_SOURCE_DIR}/src/generator/*.cc
    ${PROJECT_SOURCE_DIR}/src/output/*.cc
    ${PROJECT_SOURCE_DIR}/src/util/*.cc
    ${PROJECT_SOURCE_DIR}/src/gpu/*.cc
)
list(FILTER SOURCES EXCLUDE REGEX "/\\._")

# Conditionally include physics sources based on feature flags
file(GLOB PHYSICS_SOURCES ${PROJECT_SOURCE_DIR}/src/physics/PhysicsList.cc)
list(APPEND SOURCES ${PHYSICS_SOURCES})

# Only include Garfield++ model if enabled and found
if(WITH_GARFIELD AND Garfield_FOUND)
    list(APPEND SOURCES ${PROJECT_SOURCE_DIR}/src/physics/GarfieldModel.cc)
endif()

# Only include Celeritas interface if enabled and found
if(WITH_CELERITAS AND Celeritas_FOUND)
    list(APPEND SOURCES ${PROJECT_SOURCE_DIR}/src/physics/CeleritasInterface.cc)
    list(APPEND SOURCES ${PROJECT_SOURCE_DIR}/src/physics/CeleritasCalorimeter.cc)
endif()

# Include Opticks interface if enabled (stub mode if libs not found)
if(WITH_OPTICKS)
    list(APPEND SOURCES ${PROJECT_SOURCE_DIR}/src/physics/OpticksInterface.cc)
endif()

# Include GarfieldGPU sources (CUDA or CPU fallback)
if(WITH_GARFIELD_GPU)
    if(GARFIELD_GPU_USE_CUDA)
        # CUDA implementation
        set_source_files_properties(
            ${PROJECT_SOURCE_DIR}/src/physics/GarfieldGPU.cu
            PROPERTIES LANGUAGE CUDA
        )
        list(APPEND SOURCES ${PROJECT_SOURCE_DIR}/src/physics/GarfieldGPU.cu)
    else()
        # CPU fallback (OpenMP or single-threaded)
        list(APPEND SOURCES ${PROJECT_SOURCE_DIR}/src/physics/GarfieldGPU_cpu.cc)
    endif()
    # TPC simulation components always included when GarfieldGPU is enabled
    list(APPEND SOURCES ${PROJECT_SOURCE_DIR}/src/physics/TPCDriftManager.cc)
    list(APPEND SOURCES ${PROJECT_SOURCE_DIR}/src/physics/TPCFieldMap.cc)
    list(APPEND SOURCES ${PROJECT_SOURCE_DIR}/src/physics/TPCPadReadout.cc)
endif()

# Conditionally include GUI sources
if(WITH_DASHBOARD)
    file(GLOB GUI_SOURCES ${PROJECT_SOURCE_DIR}/src/gui/*.cc)
    list(FILTER GUI_SOURCES EXCLUDE REGEX "/\\._")

    # Explicitly run MOC on the dashboard headers with Q_OBJECT
    if(Qt6_FOUND)
        qt6_wrap_cpp(MOC_SOURCES
            ${PROJECT_SOURCE_DIR}/include/gui/DashboardWindow.hh
            ${PROJECT_SOURCE_DIR}/include/gui/EventDisplay.hh
            ${PROJECT_SOURCE_DIR}/include/gui/MaterialBudgetPlot.hh
            ${PROJECT_SOURCE_DIR}/include/gui/G4CommandPanel.hh
            ${PROJECT_SOURCE_DIR}/include/gui/DashboardPanel.hh
            ${PROJECT_SOURCE_DIR}/include/gui/TPCPadDisplay.hh
            OPTIONS -DWITH_DASHBOARD=1
        )
    elseif(Qt5_FOUND)
        qt5_wrap_cpp(MOC_SOURCES
            ${PROJECT_SOURCE_DIR}/include/gui/DashboardWindow.hh
            ${PROJECT_SOURCE_DIR}/include/gui/EventDisplay.hh
            ${PROJECT_SOURCE_DIR}/include/gui/MaterialBudgetPlot.hh
            ${PROJECT_SOURCE_DIR}/include/gui/G4CommandPanel.hh
            ${PROJECT_SOURCE_DIR}/include/gui/DashboardPanel.hh
            ${PROJECT_SOURCE_DIR}/include/gui/TPCPadDisplay.hh
            OPTIONS -DWITH_DASHBOARD=1
        )
    endif()

    list(APPEND SOURCES ${GUI_SOURCES})
    list(APPEND SOURCES ${MOC_SOURCES})
endif()

file(GLOB_RECURSE HEADERS
    ${PROJECT_SOURCE_DIR}/include/*.hh
)
list(FILTER HEADERS EXCLUDE REGEX "/\\._")

# ============================================================================
# Build Main Executable
# ============================================================================
add_executable(nnbar-detector-simulation
    ${PROJECT_SOURCE_DIR}/src/main.cc
    ${SOURCES}
    ${HEADERS}
)

target_include_directories(nnbar-detector-simulation PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/core
    ${PROJECT_SOURCE_DIR}/include/detector
    ${PROJECT_SOURCE_DIR}/include/sensitive
    ${PROJECT_SOURCE_DIR}/include/hits
    ${PROJECT_SOURCE_DIR}/include/generator
    ${PROJECT_SOURCE_DIR}/include/output
    ${PROJECT_SOURCE_DIR}/include/physics
    ${PROJECT_SOURCE_DIR}/include/gui
    ${PROJECT_SOURCE_DIR}/include/util
    ${CMAKE_BINARY_DIR}
)

target_link_libraries(nnbar-detector-simulation PRIVATE
    ${Geant4_LIBRARIES}
    parquet_writer
    mcpl
)

# Optional GPU and advanced physics libraries
if(WITH_GARFIELD AND Garfield_FOUND)
    target_link_libraries(nnbar-detector-simulation PRIVATE ${Garfield_LIBRARIES})
    target_include_directories(nnbar-detector-simulation PRIVATE ${Garfield_INCLUDE_DIRS})
endif()

if(WITH_CELERITAS AND Celeritas_FOUND)
    # Need both celeritas (core physics) and accel (Geant4 integration)
    target_link_libraries(nnbar-detector-simulation PRIVATE
        Celeritas::accel
        Celeritas::celeritas)
    target_include_directories(nnbar-detector-simulation PRIVATE ${Celeritas_INCLUDE_DIRS})
endif()

# Only link Opticks libraries if full installation found
if(WITH_OPTICKS AND OPTICKS_FULL_INSTALL)
    target_link_libraries(nnbar-detector-simulation PRIVATE Opticks::Opticks)
    if(OPTICKS_INCLUDE_DIRS)
        target_include_directories(nnbar-detector-simulation PRIVATE ${OPTICKS_INCLUDE_DIRS})
    endif()
    # Add OptiX include directory if present
    if(EXISTS "/home/billy/nnbar/simulation/GEANT4_Packages/install/optix/include")
        target_include_directories(nnbar-detector-simulation PRIVATE
            /home/billy/nnbar/simulation/GEANT4_Packages/install/optix/include)
    endif()
endif()

if(WITH_CRY)
    target_include_directories(nnbar-detector-simulation PRIVATE ${CRY_INCLUDE_DIR})
    target_link_libraries(nnbar-detector-simulation PRIVATE ${CRY_LIBRARY})
endif()

# Qt Dashboard libraries
if(WITH_DASHBOARD)
    target_link_libraries(nnbar-detector-simulation PRIVATE ${QT_LIBRARIES})
    target_include_directories(nnbar-detector-simulation PRIVATE ${PROJECT_SOURCE_DIR}/include/gui)
endif()

# GarfieldGPU libraries (CUDA or OpenMP)
if(WITH_GARFIELD_GPU)
    if(GARFIELD_GPU_USE_CUDA)
        # Set CUDA architecture (auto-detect or specify)
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES "60;70;75;80;86")
        endif()
        set_target_properties(nnbar-detector-simulation PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON
        )
        # Link CUDA runtime and curand
        target_link_libraries(nnbar-detector-simulation PRIVATE
            CUDA::cudart
            CUDA::curand
        )
    else()
        # OpenMP for CPU fallback
        if(OpenMP_CXX_FOUND)
            target_link_libraries(nnbar-detector-simulation PRIVATE OpenMP::OpenMP_CXX)
        endif()
    endif()
endif()

# ============================================================================
# Copy Runtime Files to Build Directory
# ============================================================================
set(RUNTIME_DIRS
    batch
    config
    docs
    data
    macro
    profiling_tool
    scripts
)

foreach(dir ${RUNTIME_DIRS})
    if(EXISTS ${PROJECT_SOURCE_DIR}/${dir})
        file(COPY ${PROJECT_SOURCE_DIR}/${dir}
            DESTINATION ${CMAKE_BINARY_DIR}
            FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_EXECUTE
                            WORLD_READ WORLD_EXECUTE
        )
    endif()
endforeach()

# Copy individual files
set(RUNTIME_FILES
    gui.mac
    init_vis.mac
    vis.mac
    setup.file
)

foreach(file ${RUNTIME_FILES})
    if(EXISTS ${PROJECT_SOURCE_DIR}/${file})
        configure_file(
            ${PROJECT_SOURCE_DIR}/${file}
            ${CMAKE_BINARY_DIR}/${file}
            COPYONLY
        )
    endif()
endforeach()

# ============================================================================
# Create Output Directories
# ============================================================================
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/log)

# ============================================================================
# Installation
# ============================================================================
install(TARGETS nnbar-detector-simulation
    RUNTIME DESTINATION bin
)

install(DIRECTORY ${CMAKE_BINARY_DIR}/config
    DESTINATION share/nnbar-detector-simulation
)

install(DIRECTORY ${CMAKE_BINARY_DIR}/macro
    DESTINATION share/nnbar-detector-simulation
)

# ============================================================================
# Final Configuration Summary
# ============================================================================
message("")
message("${BLUE}${BOLD}[4/4] Configuration Summary${RESET}")
message("${BLUE}────────────────────────────────────────────────────────────────────────────${RESET}")
message("")

# Simulation Features Table
message("${CYAN}┌─────────────────────────────────────────────────────────────────────────────┐${RESET}")
message("${CYAN}│${RESET}${BOLD}                         Simulation Features                                ${RESET}${CYAN}│${RESET}")
message("${CYAN}├─────────────────────────────────────────────────────────────────────────────┤${RESET}")
if(WITH_SCINTILLATION)
    message("${CYAN}│${RESET}  Optical Photons (Scintillation/Cerenkov)  │  ${GREEN}ENABLED${RESET}  │ Realistic light  ${CYAN}│${RESET}")
else()
    message("${CYAN}│${RESET}  Optical Photons (Scintillation/Cerenkov)  │  ${WHITE}DISABLED${RESET} │ Fast mode        ${CYAN}│${RESET}")
endif()
if(WITH_GARFIELD AND Garfield_FOUND)
    message("${CYAN}│${RESET}  Garfield++ TPC Simulation                 │  ${GREEN}ENABLED${RESET}  │ Electron drift   ${CYAN}│${RESET}")
else()
    message("${CYAN}│${RESET}  Garfield++ TPC Simulation                 │  ${WHITE}DISABLED${RESET} │ Simple TPC       ${CYAN}│${RESET}")
endif()
if(WITH_GARFIELD_GPU)
    if(GARFIELD_GPU_USE_CUDA)
        message("${CYAN}│${RESET}  GarfieldGPU Custom Drift                  │  ${GREEN}CUDA${RESET}     │ 100-1000x speed  ${CYAN}│${RESET}")
    elseif(OpenMP_CXX_FOUND)
        message("${CYAN}│${RESET}  GarfieldGPU Custom Drift                  │  ${YELLOW}OpenMP${RESET}   │ Multi-threaded   ${CYAN}│${RESET}")
    else()
        message("${CYAN}│${RESET}  GarfieldGPU Custom Drift                  │  ${YELLOW}CPU${RESET}      │ Single-threaded  ${CYAN}│${RESET}")
    endif()
else()
    message("${CYAN}│${RESET}  GarfieldGPU Custom Drift                  │  ${WHITE}DISABLED${RESET} │ Use Garfield++   ${CYAN}│${RESET}")
endif()
if(WITH_CELERITAS AND Celeritas_FOUND)
    message("${CYAN}│${RESET}  Celeritas GPU EM Physics                  │  ${GREEN}ENABLED${RESET}  │ 10-100x speedup  ${CYAN}│${RESET}")
else()
    message("${CYAN}│${RESET}  Celeritas GPU EM Physics                  │  ${WHITE}DISABLED${RESET} │ CPU EM physics   ${CYAN}│${RESET}")
endif()
if(WITH_OPTICKS AND OPTICKS_FULL_INSTALL)
    message("${CYAN}│${RESET}  Opticks GPU Optical Photons               │  ${GREEN}ENABLED${RESET}  │ 50-200x speedup  ${CYAN}│${RESET}")
elseif(WITH_OPTICKS)
    message("${CYAN}│${RESET}  Opticks GPU Optical Photons               │  ${YELLOW}STUB${RESET}     │ CPU fallback     ${CYAN}│${RESET}")
else()
    message("${CYAN}│${RESET}  Opticks GPU Optical Photons               │  ${WHITE}DISABLED${RESET} │ CPU optical      ${CYAN}│${RESET}")
endif()
if(WITH_CRY)
    message("${CYAN}│${RESET}  CRY Cosmic Ray Generator                │  ${GREEN}ENABLED${RESET}  │ Cosmic shower   ${CYAN}│${RESET}")
else()
    message("${CYAN}│${RESET}  CRY Cosmic Ray Generator                │  ${WHITE}DISABLED${RESET} │ Particle gun    ${CYAN}│${RESET}")
endif()
if(WITH_DASHBOARD)
    message("${CYAN}│${RESET}  Qt Monitoring Dashboard                   │  ${GREEN}ENABLED${RESET}  │ Real-time GUI    ${CYAN}│${RESET}")
else()
    message("${CYAN}│${RESET}  Qt Monitoring Dashboard                   │  ${WHITE}DISABLED${RESET} │ Terminal only    ${CYAN}│${RESET}")
endif()
message("${CYAN}├─────────────────────────────────────────────────────────────────────────────┤${RESET}")
if(MCPL_BUILD)
    message("${CYAN}│${RESET}  Particle Source                           │  ${MAGENTA}MCPL File${RESET}          │             ${CYAN}│${RESET}")
else()
    message("${CYAN}│${RESET}  Particle Source                           │  ${MAGENTA}Particle Gun${RESET}       │             ${CYAN}│${RESET}")
endif()
if(TARGET_BUILD)
    message("${CYAN}│${RESET}  Carbon Target                             │  ${GREEN}Included${RESET}           │             ${CYAN}│${RESET}")
else()
    message("${CYAN}│${RESET}  Carbon Target                             │  ${WHITE}Not Included${RESET}       │             ${CYAN}│${RESET}")
endif()
message("${CYAN}└─────────────────────────────────────────────────────────────────────────────┘${RESET}")

# Build Info
message("")
message("${CYAN}┌─────────────────────────────────────────────────────────────────────────────┐${RESET}")
message("${CYAN}│${RESET}${BOLD}                            Build Information                               ${RESET}${CYAN}│${RESET}")
message("${CYAN}├─────────────────────────────────────────────────────────────────────────────┤${RESET}")
message("${CYAN}│${RESET}  Geant4 Version:     ${GREEN}${Geant4_VERSION}${RESET}")
message("${CYAN}│${RESET}  Build Type:         ${MAGENTA}${CMAKE_BUILD_TYPE}${RESET}")
message("${CYAN}│${RESET}  Install Prefix:     ${WHITE}${CMAKE_INSTALL_PREFIX}${RESET}")
message("${CYAN}│${RESET}  Output Directory:   ${WHITE}${CMAKE_BINARY_DIR}/output${RESET}")
message("${CYAN}└─────────────────────────────────────────────────────────────────────────────┘${RESET}")

message("")
message("${GREEN}╔══════════════════════════════════════════════════════════════════════════╗${RESET}")
message("${GREEN}║${RESET}${BOLD}                    Configuration Complete!                               ${RESET}${GREEN}║${RESET}")
message("${GREEN}╚══════════════════════════════════════════════════════════════════════════╝${RESET}")
message("")
message("  ${BOLD}Next Steps:${RESET}")
message("    1. Build:   ${CYAN}cmake --build . -- -j<N>${RESET}")
message("    2. Run:     ${CYAN}./nnbar-detector-simulation -m macro/signal/run_signal.mac${RESET}")
message("")

# ============================================================================
# Post-Build: Create Wrapper Script for Automatic Library Loading
# ============================================================================
# This creates a wrapper script that sets up LD_LIBRARY_PATH automatically
# so users can run ./nnbar-detector-simulation directly without manual setup

# Set paths for wrapper script generation from the packages CMake actually used.
set(GEANT4_INSTALL "${NNBAR_GEANT4_INSTALL}")
set(CELERITAS_INSTALL "")
if(Celeritas_DIR)
    get_filename_component(CELERITAS_INSTALL "${Celeritas_DIR}/../../.." ABSOLUTE)
elseif(EXISTS "${PROJECT_SOURCE_DIR}/../GEANT4_Packages/install/celeritas")
    set(CELERITAS_INSTALL "${PROJECT_SOURCE_DIR}/../GEANT4_Packages/install/celeritas")
endif()
set(OPTICKS_INSTALL "")
if(OPTICKS_HOME)
    set(OPTICKS_INSTALL "${OPTICKS_HOME}")
elseif(EXISTS "${PROJECT_SOURCE_DIR}/../GEANT4_Packages/install/opticks")
    set(OPTICKS_INSTALL "${PROJECT_SOURCE_DIR}/../GEANT4_Packages/install/opticks")
endif()

set(NNBAR_WRAPPER_LIBRARY_DIRS "")
foreach(path
    "${GEANT4_INSTALL}/lib"
    "${CELERITAS_INSTALL}/lib"
    "${OPTICKS_INSTALL}/lib"
)
    if(EXISTS "${path}")
        list(APPEND NNBAR_WRAPPER_LIBRARY_DIRS "${path}")
    endif()
endforeach()

foreach(package_dir "${Arrow_DIR}" "${Parquet_DIR}" "${spdlog_DIR}")
    if(package_dir)
        get_filename_component(package_prefix "${package_dir}/../../.." ABSOLUTE)
        if(EXISTS "${package_prefix}/lib")
            list(APPEND NNBAR_WRAPPER_LIBRARY_DIRS "${package_prefix}/lib")
        endif()
    endif()
endforeach()

foreach(path
    "${PROJECT_SOURCE_DIR}/external/arrow-install/lib"
    "${PROJECT_SOURCE_DIR}/external/spdlog-install/lib"
)
    if(EXISTS "${path}")
        list(APPEND NNBAR_WRAPPER_LIBRARY_DIRS "${path}")
    endif()
endforeach()

if(NOT APPLE)
    foreach(path
        "${PROJECT_SOURCE_DIR}/external/arrow-install-linux/lib"
        "${PROJECT_SOURCE_DIR}/external/spdlog-install-linux/lib"
        "/usr/local/cuda-12.4/lib64"
        "/usr/lib/x86_64-linux-gnu"
    )
        if(EXISTS "${path}")
            list(APPEND NNBAR_WRAPPER_LIBRARY_DIRS "${path}")
        endif()
    endforeach()
endif()

list(REMOVE_DUPLICATES NNBAR_WRAPPER_LIBRARY_DIRS)
if(APPLE)
    set(NNBAR_LIBRARY_PATH_VAR "DYLD_LIBRARY_PATH")
    set(NNBAR_WRAPPER_LIBRARY_PATH "")
else()
    set(NNBAR_LIBRARY_PATH_VAR "LD_LIBRARY_PATH")
    string(REPLACE ";" ":" NNBAR_WRAPPER_LIBRARY_PATH "${NNBAR_WRAPPER_LIBRARY_DIRS}")
endif()

# Generate wrapper script from template
configure_file(
    ${PROJECT_SOURCE_DIR}/scripts/wrapper_template.sh.in
    ${CMAKE_BINARY_DIR}/nnbar-detector-simulation.wrapper
    @ONLY
)

add_custom_command(TARGET nnbar-detector-simulation POST_BUILD
    # Rename the binary to .bin
    COMMAND ${CMAKE_COMMAND} -E rename
        $<TARGET_FILE:nnbar-detector-simulation>
        ${CMAKE_BINARY_DIR}/nnbar-detector-simulation.bin

    # Copy the generated wrapper script
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_BINARY_DIR}/nnbar-detector-simulation.wrapper
        ${CMAKE_BINARY_DIR}/nnbar-detector-simulation

    # Make the wrapper script executable
    COMMAND chmod +x ${CMAKE_BINARY_DIR}/nnbar-detector-simulation

    COMMENT "Creating wrapper script for automatic library loading"
)
