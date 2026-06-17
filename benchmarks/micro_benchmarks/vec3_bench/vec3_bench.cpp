#include "../experimental/vec3_valarray.h"
#include "experimental/vec3_valarray.h"
#include "vec3.hpp"
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

using REAL               = double;
constexpr int  numTrials = 1;
constexpr REAL min = -1e6, max = 1e6;

#define TWO_VEC_SETUP                                                                                 \
    const std::vector<REAL> rand = GENERATE(take(numTrials, chunk<REAL>(6, random<REAL>(min, max)))); \
    [[maybe_unused]] vec3   v0{rand[0], rand[1], rand[2]};                                            \
    [[maybe_unused]] vec3   v1{rand[3], rand[4], rand[5]};

#define BASE_EXPERIMENTAL_COMPARE(bench)                                      \
    BENCHMARK_ADVANCED("base")(Catch::Benchmark::Chronometer meter) {         \
        typedef fp::vec3<REAL> vec3;                                          \
        TWO_VEC_SETUP                                                         \
        bench                                                                 \
    };                                                                        \
    BENCHMARK_ADVANCED("experimental")(Catch::Benchmark::Chronometer meter) { \
        typedef fp::experimental::vec3<REAL> vec3;                            \
        TWO_VEC_SETUP                                                         \
        bench                                                                 \
    };

TEST_CASE("Benchmark vec3", "[!benchmark]") {

    SECTION("named creation") {
        const std::vector<REAL> rand = GENERATE(take(numTrials, chunk<REAL>(6, random<REAL>(min, max))));
        BENCHMARK_ADVANCED("base")(Catch::Benchmark::Chronometer meter) {
            using fp::vec3;
            meter.measure([&rand] {
                vec3<REAL> v1{rand[3], rand[4], rand[5]};
                return v1;
            });
        };

        BENCHMARK_ADVANCED("experimental")(Catch::Benchmark::Chronometer meter) {
            using fp::experimental::vec3;
            meter.measure([&rand] {
                vec3<REAL> v1{rand[3], rand[4], rand[5]};
                return v1;
            });
        };
    }

    SECTION("unnamed creation") {
        const std::vector<REAL> rand = GENERATE(take(numTrials, chunk<REAL>(6, random<REAL>(min, max))));
        BENCHMARK_ADVANCED("base") { return fp::vec3<REAL>{rand[3], rand[4], rand[5]}; };
        BENCHMARK_ADVANCED("experimental") { return fp::experimental::vec3<REAL>{rand[3], rand[4], rand[5]}; };
    }

    SECTION("dot") {
#define dot_bench          \
    meter.measure([&] {    \
        return v0.dot(v1); \
    });
        BASE_EXPERIMENTAL_COMPARE(dot_bench)
    }

    SECTION("norm") {
#define norm_bench        \
    meter.measure([&] {   \
        return v0.norm(); \
    });
        BASE_EXPERIMENTAL_COMPARE(norm_bench)
    }

    SECTION("norm_square") {
#define norm_square_bench        \
    meter.measure([&] {          \
        return v0.norm_square(); \
    });
        BASE_EXPERIMENTAL_COMPARE(norm_square_bench)
    }

    SECTION("two vector subtract") {
#define subtraction_bench \
    meter.measure([&] {   \
        return v0 - v1;   \
    });
        BASE_EXPERIMENTAL_COMPARE(subtraction_bench)
    };

    SECTION("vector minus") {
#define vector_minus_bench \
    meter.measure([&] {    \
        return -v1;        \
    });
        BASE_EXPERIMENTAL_COMPARE(vector_minus_bench)
    };

    SECTION("two vector add") {
#define add_bench       \
    meter.measure([&] { \
        return v0 + v1; \
    });
        BASE_EXPERIMENTAL_COMPARE(add_bench)
    };

    SECTION("multiply left") {
#define multiply_left_bench \
    auto r = rand[4];       \
    meter.measure([&] {     \
        return r * v1;      \
    });
        BASE_EXPERIMENTAL_COMPARE(multiply_left_bench)
    };

    SECTION("multiply right") {
#define multiply_right_bench \
    auto r = rand[4];        \
    meter.measure([&] {      \
        return v1 * r;       \
    });
        BASE_EXPERIMENTAL_COMPARE(multiply_left_bench)
    };

    SECTION("complex expression") {
#define complex_expression_bench                          \
    auto r = rand[4];                                     \
    meter.measure([&] {                                   \
        return (v1 * r + v0 * v1.dot(v0)).cross(v0) - v0; \
    });
        BASE_EXPERIMENTAL_COMPARE(complex_expression_bench)
    };
}
