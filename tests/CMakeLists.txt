cmake_minimum_required(VERSION 3.16)
project(flippy_test DESCRIPTION "flippy test framework")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")

include_directories(.)
include_directories(external)
include_directories(../flippy)
include_directories(../flippy/external)

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
        ../flippy/MonteCarloUpdater.hpp)

enable_testing()
add_test(${PROJECT_NAME} ${PROJECT_NAME})