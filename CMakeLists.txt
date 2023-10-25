cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(nnbar-calo-sim)
add_definitions(-w)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()
if ( MSVC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()


set(PROGRAM_PERMISSIONS_DEFAULT
    OWNER_WRITE OWNER_READ OWNER_EXECUTE
    GROUP_WRITE GROUP_READ GROUP_EXECUTE
    WORLD_WRITE WORLD_READ WORLD_EXECUTE)
    
option(WITH_GEANT4_UIVIS "Build nnbar-calo-sim with Geant4 UI and Vis drivers" ON)

if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
message(boost library, ${Boost_LIBRARY_DIRS})
message(${Boost_INCLUDE_DIR})


#Arrow/Parquet
find_package(Arrow REQUIRED)
message("------ ARROW Found : " ${Arrow_FOUND} , ${ARROW_FULL_SO_VERSION})
message("------ ARROW --------",arrow) 
#/sw/easybuild/software/Arrow/6.0.0-foss-2021b/lib/)


# NNBAR full detector configuration variables

if(NOT DEFINED MCPL_BUILD)
set (MCPL_BUILD 0)
message("--MCPL_BUILD not specified -- Set default value 0 (Particle gun mode)")
endif()

if(NOT DEFINED DEBUG_VERBOSE)
set (DEBUG_VERBOSE 0)
message("--DEBUG_VERBOSE not specified -- Set default value 0 (No Debug Output)")
endif()

if(NOT DEFINED TARGET_BUILD)
set (TARGET_BUILD 1)
message("--TARGET_BUILD not specified -- Set default value 1 (TARGET Foil will be built)")
endif()

message("= = = = = = = = = = = = = = = = = =")
message("Please check the following setting:")
message("MCPL_BUILD mode = ${MCPL_BUILD}")
message("TARGET build = ${TARGET_BUILD}")
message("DEBUG mode = ${DEBUG_VERBOSE}")
message("= = = = = = = = = = = = = = = = = =")

configure_file (config.h.in config.h) #exporting the config !! 

#pq writer
set(PQ_INCLUDE_DIRS "/home/scyiu/nnbar/parquet-writer-main/src/cpp")
set(PQ_LIBRARY_DIR "/home/scyiu/nnbar/parquet-writer-main/build/lib")


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#

include_directories(${PROJECT_SOURCE_DIR}/include 
        ${PQ_INCLUDE_DIRS}
		    ${CMAKE_CURRENT_BINARY_DIR}
  	    )

include_directories(/home/scyiu/nnbar/parquet-writer-main/src/cpp)
include_directories(/home/scyiu/nnbar/parquet-writer-main/exts/nlohmann/)
include_directories(/home/scyiu/nnbar/parquet-writer-main/exts/spdlog/)


file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*
                  ${PQ_INCLUDE_DIRS}/*.cpp
                  ${PROJECT_SOURCE_DIR}/src/Detector_Module/*
                  ${PROJECT_SOURCE_DIR}/src/Sensitive_Detector/*)
                  
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*
                  ${PQ_INCLUDE_DIRS}/*.h
                  ${PROJECT_SOURCE_DIR}/include/Detector_Module/*
                  ${PROJECT_SOURCE_DIR}/include/Sensitive_Detector/*)

add_executable(nnbar-calo-sim nnbar-calo-sim.cc ${sources} ${headers})
target_link_libraries(nnbar-calo-sim ${Geant4_LIBRARIES} arrow parquet ${PQ_LIBRARY_DIR} ${Boost_LIBRARY_DIRS} ${Boost_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

file(GLOB BATCH_FILES "batch_scripts/*")
file(GLOB CONFIG_FILES "config/*")
file(GLOB DOC_FILES "Documentations/*")
file(GLOB LEAD_GLASS_FILES "lead_glass_position/*")
file(GLOB MACRO_FILES "macro/*")
file(GLOB MCPL_FILES "mcpl_files/*")
file(GLOB PROFILE_FILES "profiling_tool/*")

file(COPY ${BATCH_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/batch_scripts
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${CONFIG_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/config
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${DOC_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Documentations
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${LEAD_GLASS_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/lead_glass_position
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${MACRO_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/macro
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${MCPL_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mcpl_files
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${PROFILE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/profiling_tool
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

set(SCRIPTS
  nnbar-calo-sim.in
  gui.mac
  init_vis.mac
  vis.mac
  setup.file
  )

foreach(_script ${SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS nnbar-calo-sim DESTINATION bin)
install(
	DIRECTORY
	DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/output";
  	DIRECTORY_PERMISSIONS ${PROGRAM_PERMISSIONS_DEFAULT}
)
install(
	DIRECTORY
	DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/log";
  	DIRECTORY_PERMISSIONS ${PROGRAM_PERMISSIONS_DEFAULT}
)
  	
