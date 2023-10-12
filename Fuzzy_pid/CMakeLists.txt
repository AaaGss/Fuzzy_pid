cmake_minimum_required(VERSION 3.10)

# Setup the project.
# Its name is defined to be the controller directory name.
set(PROJECT "main")
project(${PROJECT})

set( CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set( CMAKE_CXX_STANDARD 14)
set( CMAKE_C_STANDARD 11)

# # Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message("Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for test_footStest_footPlacementtepPlan-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()


# include_directories(./inc)
# include_directories(./)
set(INCLUDE_DIRECTORIES ./inc ./)
include_directories(${INCLUDE_DIRECTORIES})

# Get C or C++ sources in the current directory (only).
file(GLOB C_SOURCES ./lib/*.c ./*.c)
file(GLOB CPP_SOURCES ./lib/*.cpp ./*.cpp)
set(SOURCES ${C_SOURCES} ${CPP_SOURCES})


# create and Link with the library.
#set( LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib) 
#add_library(circleArea SHARED ./lib/circleArea.c) #想得到动态库，参数就是SHARED
#link_directories(${PROJECT_SOURCE_DIR}/lib)
#set (LIBRARIES ${CMAKE_SHARED_LIBRARY_PREFIX}circleArea${CMAKE_SHARED_LIBRARY_SUFFIX} )


# Setup the target executable.
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}) 
add_executable(${PROJECT} ${SOURCES})
target_link_libraries(${PROJECT} ${LIBRARIES})


