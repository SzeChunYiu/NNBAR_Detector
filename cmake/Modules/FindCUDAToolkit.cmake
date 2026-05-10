# FindCUDAToolkit.cmake - Simplified version for CMake 3.16 compatibility
# This finds CUDA toolkit libraries for projects built with older CMake

# Find CUDA compiler
find_program(CUDAToolkit_NVCC_EXECUTABLE
    NAMES nvcc
    HINTS
        ${CUDAToolkit_ROOT}
        $ENV{CUDA_HOME}
        $ENV{CUDA_PATH}
        /usr/local/cuda-12.4
        /usr/local/cuda
    PATH_SUFFIXES bin
)

if(CUDAToolkit_NVCC_EXECUTABLE)
    get_filename_component(CUDAToolkit_BIN_DIR "${CUDAToolkit_NVCC_EXECUTABLE}" DIRECTORY)
    get_filename_component(CUDAToolkit_ROOT_DIR "${CUDAToolkit_BIN_DIR}" DIRECTORY)

    # Get version from nvcc
    execute_process(
        COMMAND ${CUDAToolkit_NVCC_EXECUTABLE} --version
        OUTPUT_VARIABLE NVCC_VERSION_OUTPUT
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(NVCC_VERSION_OUTPUT MATCHES "release ([0-9]+)\\.([0-9]+)")
        set(CUDAToolkit_VERSION_MAJOR "${CMAKE_MATCH_1}")
        set(CUDAToolkit_VERSION_MINOR "${CMAKE_MATCH_2}")
        set(CUDAToolkit_VERSION "${CUDAToolkit_VERSION_MAJOR}.${CUDAToolkit_VERSION_MINOR}")
    endif()

    # Set include directory
    set(CUDAToolkit_INCLUDE_DIRS "${CUDAToolkit_ROOT_DIR}/include")

    # Set library directory
    if(EXISTS "${CUDAToolkit_ROOT_DIR}/lib64")
        set(CUDAToolkit_LIBRARY_DIR "${CUDAToolkit_ROOT_DIR}/lib64")
    else()
        set(CUDAToolkit_LIBRARY_DIR "${CUDAToolkit_ROOT_DIR}/lib")
    endif()

    # Find required libraries
    find_library(CUDA_cudart_LIBRARY
        NAMES cudart
        HINTS ${CUDAToolkit_LIBRARY_DIR}
    )

    find_library(CUDA_curand_LIBRARY
        NAMES curand
        HINTS ${CUDAToolkit_LIBRARY_DIR}
    )

    find_library(CUDA_cuda_driver_LIBRARY
        NAMES cuda
        HINTS ${CUDAToolkit_LIBRARY_DIR} /usr/lib/x86_64-linux-gnu
    )

    # Create imported targets
    if(CUDA_cudart_LIBRARY AND NOT TARGET CUDA::cudart)
        add_library(CUDA::cudart SHARED IMPORTED)
        set_target_properties(CUDA::cudart PROPERTIES
            IMPORTED_LOCATION "${CUDA_cudart_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${CUDAToolkit_INCLUDE_DIRS}"
        )
    endif()

    if(CUDA_curand_LIBRARY AND NOT TARGET CUDA::curand)
        add_library(CUDA::curand SHARED IMPORTED)
        set_target_properties(CUDA::curand PROPERTIES
            IMPORTED_LOCATION "${CUDA_curand_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${CUDAToolkit_INCLUDE_DIRS}"
        )
    endif()

    if(CUDA_cuda_driver_LIBRARY AND NOT TARGET CUDA::cuda_driver)
        add_library(CUDA::cuda_driver SHARED IMPORTED)
        set_target_properties(CUDA::cuda_driver PROPERTIES
            IMPORTED_LOCATION "${CUDA_cuda_driver_LIBRARY}"
        )
    endif()

    # Set found status
    if(CUDA_cudart_LIBRARY)
        set(CUDAToolkit_FOUND TRUE)
    endif()
endif()

# Handle standard find_package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CUDAToolkit
    REQUIRED_VARS CUDAToolkit_NVCC_EXECUTABLE CUDAToolkit_INCLUDE_DIRS CUDA_cudart_LIBRARY
    VERSION_VAR CUDAToolkit_VERSION
)

# Export variables
if(CUDAToolkit_FOUND)
    set(CUDAToolkit_INCLUDE_DIRS "${CUDAToolkit_INCLUDE_DIRS}" CACHE PATH "CUDA include dirs")
    set(CUDAToolkit_LIBRARY_DIR "${CUDAToolkit_LIBRARY_DIR}" CACHE PATH "CUDA library dir")
    set(CUDAToolkit_LIBRARIES "${CUDA_cudart_LIBRARY}" CACHE FILEPATH "CUDA libraries")
endif()

mark_as_advanced(
    CUDAToolkit_NVCC_EXECUTABLE
    CUDA_cudart_LIBRARY
    CUDA_curand_LIBRARY
    CUDA_cuda_driver_LIBRARY
)
