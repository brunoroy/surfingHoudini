
#
# This cmake script will build a basic Houdini DSO.
#
cmake_minimum_required(VERSION 2.6)

set(ProjectName SurfingSAT)

############## SET ############
set(SRC .)
set(ProjectSource ${SRC}/HoudiniSurfaceFromPoints)
set(Eigen ${SRC}/eigen)

project(${ProjectName})

# run our helper process to grab the compile and link flags from hcustom
execute_process(COMMAND hdk-config.py --notag --cmake=/usr/share/cmake-2.8/Modules/FindHDK.cmake)

############# FIND PACKAGE ###############
find_package(HDK)

############## HOUDINI ###################
# now get the HDK definitions from the FindHDK.cmake module that was just generated
set(CMAKE_MODULE_PATH "${PROJECT_BINARY_DIR}")



################# SOURCE ##################
file(GLOB source ${ProjectSource}/*.cpp ${ProjectSource}/*.C)
#set (source ${ProjectSource}/HoudiniUniformEdges.cpp)

set(CMAKE_INSTALL_PREFIX "${HDK_HIH_DIR}")


############# SUBDIRECTORY ###############
#ADD_SUBDIRECTORY(${Math} "${CMAKE_CURRENT_BINARY_DIR}/src/Math")

# setup to compile in the correct HDK environment
add_definitions(${HDK_DEFINITIONS})


########### INCLUDES #################

include_directories(${HDK_INCLUDE_DIRS})

include_directories(${SRC})
include_directories(${Eigen})
include_directories(${HDK_Include})
#include_directories(${Math})


# Prevents redefinition of MIN and MAX macros in windows.h (which causes conflicts with eigen library)
add_definitions(-DNOMINMAX)

########## LINKS ######################

link_directories(${HDK_LIBRARY_DIRS})

# this command calls hdk-config.py to build a special file that will be the repository of the taginfo
add_custom_command(OUTPUT sesitag.C
                   DEPENDS ${source}
                   COMMAND ../hdk-config.py --tagfile=${CMAKE_BINARY_DIR}/sesitag.C)

if (APPLE)
    set(LTYPE MODULE)
else (APPLE)
    set(LTYPE SHARED)
endif (APPLE)


#External Libraries
#add_library(fluidsCore libfluidsCore)


# this tells cmake to buid the DSO
add_library(${ProjectName} ${LTYPE} ${source}  sesitag.C)

# houdini doesn't want the DSO to have "lib" as it's prefix
set_target_properties(${ProjectName} PROPERTIES PREFIX "")

# unfortunately cmake's MODULE library type generates a ".so" file on OSX
if (APPLE)
    set_target_properties(${ProjectName} PROPERTIES SUFFIX ".dylib")
endif (APPLE)

#target_link_libraries (${ProjectName} ${FluidsCore})
target_link_libraries(${ProjectName} 
	${HDK_LIBRARIES}
)

# finally tell cmake where to put the results
install(TARGETS ${ProjectName} LIBRARY DESTINATION dso)
install(FILES ${ProjectName}.png DESTINATION config/Icons)


