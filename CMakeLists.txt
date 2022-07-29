cmake_minimum_required(VERSION 3.13.0 FATAL_ERROR) 

# Project name
project(SIGTestProject_allpccloud CXX)

message(STATUS PROJECT_NAME: ${PROJECT_NAME}) 

option(USE_STATIC "USE STATIC libary" ON)
message(STATUS "USE_STATIC: ${USE_STATIC}")

message(STATUS BUILD_PATCH: ${BUILD_PATCH})

message(STATUS WIN32: ${WIN32})
message(STATUS UNIX ${UNIX})
message(STATUS "OS: ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_SYSTEM_VERSION: ${CMAKE_SYSTEM_VERSION}")
message(STATUS "CMAKE_HOST_SYSTEM_NAME: ${CMAKE_HOST_SYSTEM_NAME}")
message(STATUS "CMAKE_SYSTEM_PROCESSOR: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "CMAKE_GENERATOR_PLATFORM: ${CMAKE_GENERATOR_PLATFORM}")

message(STATUS "CMAKE_SYSROOT: ${CMAKE_SYSROOT}")
message(STATUS "CMAKE_STAGING_PREFIX: ${CMAKE_STAGING_PREFIX}")

message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

if (WIN32)
    message(STATUS WIN32: ${WIN32})
endif()
# OpenCascade
find_package(OpenCASCADE)
# Configure C++ compiler's includes dir
include_directories(SYSTEM ${OpenCASCADE_INCLUDE_DIR})


# Add executable
add_executable(SIGTestProject_allpccloud main.cpp)

# set_target_properties(
#         SIGTestProject_allpccloud PROPERTIES
#         CXX_STANDARD 11
#         CXX_STANDARD_REQUIRED ON
# ) 

#for libdxf
include_directories(SYSTEM /home/oleksii/OpenCASCADE/libdxf/src/)

# Add linker options
foreach(LIB ${OpenCASCADE_LIBRARIES})
    target_link_libraries(SIGTestProject_allpccloud debug ${OpenCASCADE_LIBRARY_DIR}d/lib${LIB}.so)
    target_link_libraries(SIGTestProject_allpccloud optimized ${OpenCASCADE_LIBRARY_DIR}/lib${LIB}.so)
    #message(STATUS "OpenCASCADE_INCLUDE_DIR = ${OpenCASCADE_INCLUDE_DIR}")
    #message(STATUS "LIB = ${LIB}")
endforeach(LIB ${OpenCASCADE_LIBRARIES})

# Adjust runtime environment
set_property(TARGET SIGTestProject_allpccloud PROPERTY VS_DEBUGGER_ENVIRONMENT "PATH=$<$<CONFIG:DEBUG>:${OpenCASCADE_BINARY_DIR}d>$<$<NOT:$<CONFIG:DEBUG>>:${OpenCASCADE_BINARY_DIR}>;%PATH%")

set_target_properties(SIGTestProject_allpccloud PROPERTIES   
  #VS_DEBUGGER_COMMAND "debug_command"
    VS_DEBUGGER_COMMAND_ARGUMENTS "/home/oleksii/OpenCASCADE/SIGTestProject_allpccloud/RemfHead_Flat_D4p0_75329419.stp")