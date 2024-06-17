cmake_minimum_required(VERSION 3.18)
project(vec3_microbenchmark)

set(CMAKE_CXX_STANDARD 20)
set(CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release)

add_compile_definitions(TEST_ASSET_PATH="${CMAKE_CURRENT_SOURCE_DIR}/assets/")
if (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    MESSAGE("building for MSVC")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
endif ()
message(${CMAKE_CXX_COMPILER_ID})
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")

include_directories(.)
include_directories(external)

Include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.3.2 # or a later release
)

FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)

# Create an interface library to suppress warnings for Catch2

if (${SINGLE_HEADER_TEST})
    message("testing SINGLE_HEADER_FLIPPY")
    include_directories(../single_header_flippy)
else()
    message("testing FLIPPY")
    include_directories(../flippy)
    include_directories(../flippy/external)
endif ()



message(${catch2_SOURCE_DIR})
add_executable(${PROJECT_NAME}
        vec3_bench.cpp
        node_neighbour_insertion_delition_bench.cpp
)

catch_discover_tests(${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} PRIVATE Catch2::Catch2WithMain)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Werror)
enable_testing()
add_test(${PROJECT_NAME} ${PROJECT_NAME})

