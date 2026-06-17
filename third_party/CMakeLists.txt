include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.7.1)
FetchContent_MakeAvailable(Catch2)

target_compile_options(Catch2 PRIVATE -w)
target_compile_options(Catch2WithMain PRIVATE -w)
