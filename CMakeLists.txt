cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(nnbar-calo-sim)

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

# NNBAR full detector configuration variables

if(NOT DEFINED VERSION_SHIELD)
set (VERSION_SHIELD 0)
message("--Shield build not specified -- Set default value 0 (No shield)")
endif()

if(NOT DEFINED VERSION_MCPL)
set (VERSION_MCPL 0)
message("--Option MCPL not specified -- Set default value 0 (No MCPL generator)")
endif()

if(NOT DEFINED CRY_BUILD)
set (CRY_BUILD 0)
message("--Option CRY not specified -- Set default value 0 (No CRY generator)")
endif()

message("= = = = = = = = = = = = = = = = = =")
message("Please check the following setting:")
message("Shield build = ${VERSION_SHIELD}")
message("MCPL generator = ${VERSION_MCPL}")
message("CRY generator = ${CRY_BUILD}")
message("= = = = = = = = = = = = = = = = = =")

configure_file (config.h.in config.h) #exporting the config !! 

# - - - - - - - - Setup CRY - - - - - - - - - -
if(CRY_BUILD)
MESSAGE ("Setting up Cry")
endif()

find_library(CRY NAMES CRY PATHS $ENV{CRYHOME}/lib NO_DEFAULT_PATH)

if(NOT DEFINED CRY_PATH)
message("WARNING:: Specify the path to folder -- example: CMake ... -DCRY_PATH=/home/billy/nnbar/cry_v1.7 ...")
message("CMake will use the example path since it is not specified")
set(CRY_PATH "/home/billy/nnbar/cry_v1.7")
endif()

set(CRY_INCLUDE_DIRS "${CRY_PATH}/src")
set(CRY_LIBRARY_DIR ${CRY_PATH}/lib)
set(CRY_LIBRARIES "${CRY_PATH}/lib/libCRY.a")

message(${CRY_INCLUDE_DIRS})
message(${CRY_LIBRARY_DIR})
message(${CRY_LIBRARIES})


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#

include_directories(${PROJECT_SOURCE_DIR}/include 
		    ${CRY_INCLUDE_DIRS}
		    ${CMAKE_CURRENT_BINARY_DIR}
		   )

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*)

add_executable(nnbar-calo-sim nnbar-calo-sim.cc ${sources} ${headers})

target_link_libraries(nnbar-calo-sim ${Geant4_LIBRARIES} ${CRY_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
file(GLOB MCPL_FILES "mcpl_files/*")
file(GLOB LEAD_GLASS_FILES "lead_glass_position/*")
file(COPY ${MCPL_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mcpl_files
	FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${LEAD_GLASS_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/lead_glass_position
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
	DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/output"
     	DIRECTORY_PERMISSIONS ${PROGRAM_PERMISSIONS_DEFAULT}
     )
