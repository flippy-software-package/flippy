#ifndef NODES_HPP_TEST_HELPER_H
#define NODES_HPP_TEST_HELPER_H
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

using Approx = Catch::Approx;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

[[maybe_unused]] static float      SMALL_FRACTION = 1e-5f;
[[maybe_unused]] static const auto ZERO_ISH       = WithinAbs(0, 0.00001);
[[maybe_unused]] static auto       ish(auto num) {
    return (WithinRel(static_cast<float>(num), static_cast<float>(SMALL_FRACTION)));
}

[[maybe_unused]] static auto operator""_ish(long double value) { return ish(value); }
[[maybe_unused]] static auto operator""_ish(unsigned long long value) { return ish(value); }

#endif // NODES_HPP_TEST_HELPER_H
