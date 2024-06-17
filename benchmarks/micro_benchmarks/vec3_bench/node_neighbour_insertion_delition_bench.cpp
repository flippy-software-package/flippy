#include "flippy.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/generators/catch_generators_all.hpp>

using REAL = double;
using INDEX = unsigned int;

#define SETUP \
    fp::vec3<REAL> v1{0., 1., 1.}; \
    fp::vec3<REAL> v2{1., 0., 7.}; \
    fp::vec3<REAL> v3{1., 0., 0.}; \
    fp::vec3<REAL> v4{0., 2., 3.}; \
    fp::vec3<REAL> v5{0., 0., 9.}; \
    fp::vec3<REAL> v6{1., 1., 3.}; \
    fp::vec3<REAL> v7{0., 3., 3.}; \
    fp::vec3<REAL> v8{1., 0., 3.}; \
    Node node0; \
    node0.id = 0;                  \
    node0.nn_ids = {1, 2, 3 ,4 ,5 ,6 ,7 ,8}; \
    node0.nn_distances = {v1, v2, v3, v4, v5, v6, v7, v8};

#define BASE_EXPERIMENTAL_COMPARE(bench)             \
    BENCHMARK_ADVANCED("base")(Catch::Benchmark::Chronometer meter) {       \
        typedef fp::Node<REAL,INDEX> Node; \
        SETUP                    \
        bench                     \
    };
//    BENCHMARK_ADVANCED("experimental")(Catch::Benchmark::Chronometer meter) { \
//        typedef fp::experimental::Node<REAL, INDEX> Node; \
//        SETUP \
//        bench \
//    };


TEST_CASE("Benchmark node", "[!benchmark]")
{


    INDEX new_id{11};
    fp::vec3<REAL> new_pos{0., 0., 0};


    SECTION("emplace at 0") {
#define emplace_at_0 meter.measure([&] { return node0.emplace_nn_id(new_id, new_pos, 0); });
        BASE_EXPERIMENTAL_COMPARE(emplace_at_0)
    };

    SECTION("emplace at the end") {
#define emplace_at_the_end \
    auto end = static_cast<unsigned int>(node0.nn_ids.size() - 1); \
    meter.measure([&] { return node0.emplace_nn_id(new_id, new_pos, end); });
        BASE_EXPERIMENTAL_COMPARE(emplace_at_the_end)
    };

    SECTION("emplace in the middle") {
#define emplace_in_the_middle \
    auto end = static_cast<unsigned int>(node0.nn_ids.size()/2); \
    meter.measure([&] { return node0.emplace_nn_id  (new_id, new_pos, end); });
        BASE_EXPERIMENTAL_COMPARE(emplace_in_the_middle)
    };

}