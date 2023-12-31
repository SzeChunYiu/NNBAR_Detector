cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(parquet-writer)

option(BUILD_TESTS "Whether or not to build the tests" OFF)
option(BUILD_DOCS "Whether or not to build the documentation" OFF)

set(CMAKE_CXX_STANDARD 17)

# need to handle this better
if(NOT DEFINED ENV{AtlasProject})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    execute_process(
        COMMAND  bash -c "awk -F= '/^ID=/{print $2}' /etc/os-release |tr -d '\n' | tr -d '\"'"
        OUTPUT_VARIABLE LINUX_DISTRO
    )
endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
endif()

# documentation
if(BUILD_DOCS)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindSphinx.cmake RESULT_VARIABLE FIND_SPHINX_FOUND)
    add_subdirectory(docs/)
endif()

find_package(PkgConfig REQUIRED)
find_package(Arrow REQUIRED)
find_package(Parquet REQUIRED)

add_subdirectory(exts/)

# don't propagate the tight compiler flags to the external packages
# so add them here
if(NOT DEFINED ENV{AtlasProject})
    #add_compile_options(-Wshadow)
    add_compile_options(-Werror)
    add_compile_options(-Wunused)
    add_compile_options(-Wall)
endif()

add_subdirectory(src/cpp/)
