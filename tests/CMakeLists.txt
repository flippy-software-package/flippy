add_compile_definitions(TEST_ASSET_PATH="${CMAKE_CURRENT_SOURCE_DIR}/assets/")

include_directories(.)
include_directories(external)

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(Catch)

message("testing FLIPPY")
include_directories(../flippy)
include_directories(../flippy/external)

option(FLIPPY_ENABLE_ASAN "Enable AddressSanitizer for test targets" OFF)
option(FLIPPY_ENABLE_TSAN "Enable ThreadSanitizer for test targets" OFF)
option(FLIPPY_ENABLE_UBSAN "Enable MemorySanitizer for test targets" OFF)

set(FP_TEST_COMPILE_OPTIONS ${FLIPPY_COMPILE_OPTIONS})

if(FLIPPY_ENABLE_ASAN)
  set(FP_TEST_COMPILE_OPTIONS ${FP_TEST_COMPILE_OPTIONS} -fsanitize=address)
  set(FP_TEST_LINK_OPTIONS ${FP_TEST_LINK_OPTIONS} -fsanitize=address)
elseif(FLIPPY_ENABLE_TSAN)
  set(FP_TEST_COMPILE_OPTIONS ${FP_TEST_COMPILE_OPTIONS} -fsanitize=thread)
  set(FP_TEST_LINK_OPTIONS ${FP_TEST_LINK_OPTIONS} -fsanitize=thread)
elseif(FLIPPY_ENABLE_UBSAN)
  set(FP_TEST_COMPILE_OPTIONS ${FP_TEST_COMPILE_OPTIONS} -fsanitize=undefined)
  set(FP_TEST_LINK_OPTIONS ${FP_TEST_LINK_OPTIONS} -fsanitize=undefined)
endif()

set(FLIPPY_TEST_SOURCES
    Nodes_test.cpp
    vec3_test.cpp
    local_geometry_test.cpp
    Triangulation_test.cpp
    Triangulator_test.cpp
    Triangulation_update_test.cpp
    stlSerializer_test.cpp
    sim_utils_test.cpp
    dynamic_move_size_test.cpp
    xyzDataIO_test.cpp
    local_update_energy_test.cpp)

foreach(test_source ${FLIPPY_TEST_SOURCES})
  get_filename_component(test_name ${test_source} NAME_WE)

  add_executable(${test_name} ${test_source})

  target_compile_options(${test_name} PRIVATE ${FP_TEST_COMPILE_OPTIONS})
  target_link_options(${test_name} PRIVATE ${FP_TEST_LINK_OPTIONS})
  target_include_directories(${test_name} SYSTEM
                             PRIVATE ${Catch2_SOURCE_DIR}/src)
  target_include_directories(
    ${test_name} SYSTEM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../flippy/external)
  target_link_libraries(${test_name} PRIVATE Catch2::Catch2WithMain)
  catch_discover_tests(${test_name})
endforeach()
