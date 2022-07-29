cmake_minimum_required(VERSION 3.13.0 FATAL_ERROR)

project(libdxf VERSION 0.0.1)

configure_file(
        config.h.in 
        ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

set(SOURCES 
        main.cpp 
        version.cpp
)

set(HEADERS
        version.h
)

add_library(libdxf sumLib.cpp sumLib.h)
set_target_properties(
        libdxf PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED ON
)
target_include_directories(
        libdxf 
        PRIVATE 
        ${CMAKE_CURRENT_BINARY_DIR}
)
target_include_directories(
        libdxf 
        PRIVATE 
        ${CMAKE_CURRENT_BINARY_DIR}/src
)

install(TARGETS libdxf EXPORT libdxfTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include
)
install(FILES sumLib.h DESTINATION include)

install(EXPORT libdxfTargets
        FILE libdxfConfig.cmake
        NAMESPACE LibDXF::
        DESTINATION lib/cmake/LibDXF
)