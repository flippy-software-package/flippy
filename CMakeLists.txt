cmake_minimum_required(VERSION 3.16)
project(flippy)
file(READ "VERSION.json" ver)
string(JSON ver GET ${ver} VERSION)
string(REGEX MATCHALL "[0-9]+.[0-9]+.[0-9]+" version ${ver})
string(REGEX MATCHALL "-.+" pre_release_flag ${ver})
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if("${pre_release_flag}" STREQUAL "")
  message("no prerelease flags to set")
  project(flippy VERSION ${version})
else()
  project(
    flippy
    VERSION ${version}
    DESCRIPTION ${pre_release_flag})
  message("prerelease flag: ${pre_release_flag}")
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  message("building for MSVC")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -fno-omit-frame-pointer")

set(FLIPPY_COMPILE_OPTIONS
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
    -Wswitch-enum
    -Werror
    -ffast-math
    # gToDo: remove this after fixing nlohmann dependency
    -Wno-nan-infinity-disabled
    # gToDo: remove this after fixing catch2 dependency
    -Wno-c2y-extensions)

set(CMAKE_CXX_STANDARD 26)
set(CXX_STANDARD_REQUIRED ON)

enable_testing()

add_subdirectory(third_party)
add_subdirectory(tests)
add_subdirectory(demo)
add_subdirectory(benchmarks)
