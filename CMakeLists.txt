cmake_minimum_required(VERSION 3.13.0 FATAL_ERROR) 

# Project name
project(SIGTestProject_allpccloud CXX)

# OpenCascade
find_package(OpenCASCADE)

# Configure C++ compiler's includes dir
include_directories(SYSTEM ${OpenCASCADE_INCLUDE_DIR})
# Add executable
add_executable(SIGTestProject_allpccloud main.cpp)

# Add linker options
foreach(LIB ${OpenCASCADE_LIBRARIES})
    target_link_libraries(SIGTestProject_allpccloud debug ${OpenCASCADE_LIBRARY_DIR}d/${LIB}.lib)
    target_link_libraries(SIGTestProject_allpccloud optimized ${OpenCASCADE_LIBRARY_DIR}/${LIB}.lib)
    #message(STATUS "OpenCASCADE_INCLUDE_DIR = ${OpenCASCADE_INCLUDE_DIR}")
    #message(STATUS "LIB = ${LIB}")
endforeach(LIB ${OpenCASCADE_LIBRARIES})

# Adjust runtime environment
set_property(TARGET SIGTestProject_allpccloud PROPERTY VS_DEBUGGER_ENVIRONMENT "PATH=$<$<CONFIG:DEBUG>:${OpenCASCADE_BINARY_DIR}d>$<$<NOT:$<CONFIG:DEBUG>>:${OpenCASCADE_BINARY_DIR}>;%PATH%")

set_target_properties(SIGTestProject_allpccloud PROPERTIES   
  #VS_DEBUGGER_COMMAND "debug_command"
    VS_DEBUGGER_COMMAND_ARGUMENTS "C:/OpenCASCADE/SIGTestProject_allpccloud/RemfHead_Flat_D4p0_75329419.stp")