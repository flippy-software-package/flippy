set(CMAKE_BUILD_TYPE Release)

set(demo_name biconcave_shapes_mc_demo)

include_directories(../../flippy)

add_executable(${demo_name} main.cpp)
target_compile_options(${demo_name} PRIVATE ${FLIPPY_COMPILE_OPTIONS})
