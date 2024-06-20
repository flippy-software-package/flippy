cmake_minimum_required(VERSION 3.10)

project(planar_membrane_sheet_fluctuation_mc_demo)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_definitions(-DDEBUG)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 ")

include_directories(../../flippy)
add_executable(${PROJECT_NAME} main.cpp)
