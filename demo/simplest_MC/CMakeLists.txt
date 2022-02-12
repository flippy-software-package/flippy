cmake_minimum_required(VERSION 3.10)

project(flippy_simplest_vesicle_mc_demo DESCRIPTION "flippy demo create a simple triangulated sphere")
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 ")

add_executable(${PROJECT_NAME} main.cpp)
