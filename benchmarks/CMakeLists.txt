cmake_minimum_required(VERSION 3.16)

# include_directories(micro_benchmarks)
include_directories(external/google_benchmark)
include_directories(external)
include_directories(benchmark_logger)
include_directories(../flippy)
include_directories(../flippy/utilities)
include_directories(.)

# add_subdirectory(micro_benchmarks/vec3_bench)
add_subdirectory(rbc_bench)
