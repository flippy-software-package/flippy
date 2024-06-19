cmake_minimum_required(VERSION 3.16)

file(READ "VERSION.json" ver)
string(JSON ver GET ${ver} VERSION)
string(REGEX MATCHALL "[0-9]+.[0-9]+.[0-9]+" version ${ver})
string(REGEX MATCHALL "-.+" pre_release_flag ${ver})

if("${pre_release_flag}" STREQUAL "")
    message("no prerelease flags to set")
    project(flippy VERSION ${version})
else()
    project(flippy VERSION ${version} DESCRIPTION ${pre_release_flag})
    message("prerelease flag: ${pre_release_flag}")
endif()

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")

add_subdirectory(tests)
add_subdirectory(demo)
add_subdirectory(benchmarks)