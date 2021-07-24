cmake_minimum_required(VERSION 3.16)
project(flippy_test DESCRIPTION "flippy test framework")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE Debug)
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
        )

set_target_properties( ${PROJECT_NAME} PROPERTIES DESCRIPTION ${PROJECT_DESCRIPTION})
