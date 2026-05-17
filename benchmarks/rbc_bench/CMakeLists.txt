set(CMAKE_BUILD_TYPE Release)

set(benchmark_name flippy_rbc_benchmark)

add_definitions(-DLOGGING_ON=1)
add_definitions(-DPRINTING_ON=1)
add_definitions(-DALLOCATIONTRACKING=1)

include_directories(../benchmark_logger)
include_directories(../external)
include_directories(../../flippy)

add_executable(${benchmark_name} main.cpp)

target_compile_options(${benchmark_name} PRIVATE ${FLIPPY_COMPILE_OPTIONS})
