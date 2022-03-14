cmake_minimum_required(VERSION 3.16)
project(flippy_test DESCRIPTION "flippy test framework")
set(CMAKE_CXX_STANDARD 20)
set(CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release)

if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    MESSAGE("building for MSVC")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
endif ()
message(${CMAKE_CXX_COMPILER_ID})
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")

include_directories(.)
include_directories(external)


if (${SINGLE_HEADER_TEST})
    message("testing SINGLE_HEADER_FLIPPY")
    include_directories(../single_header_flippy)
else()
    message("testing FLIPPY")
    include_directories(../flippy)
    include_directories(../flippy/external)
endif ()

set_source_files_properties(../external
        PROPERTIES
        COMPILE_FLAGS -W0
        )

add_executable(${PROJECT_NAME} main.cpp
        vec3_test.cpp
        Nodes_test.cpp
        local_geometry_test.cpp
        Triangulation_test.cpp
        Triangulator_test.cpp
        )

enable_testing()
add_test(${PROJECT_NAME} ${PROJECT_NAME})