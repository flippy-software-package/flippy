/*
 *```txt
 *
 *  .d888 888 d8b
 * d88P"  888 Y8P
 * 888    888
 * 888888 888 888 88888b.  88888b.  888  888
 * 888    888 888 888 "88b 888 "88b 888  888     simulating package for
 * 888    888 888 888  888 888  888 888  888     dynamically triangulated
 * 888    888 888 888 d88P 888 d88P Y88b 888     surfaces
 * 888    888 888 88888P"  88888P"   "Y88888
 *                888      888           888
 *                888      888      Y8b d88P
 *                888      888       "Y88P"
 *
 * https://github.com/flippy-software-package/flippy
 *
 *
 * MIT License
 *
 * Copyright (c) 2021 George Dadunashvili
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *```
 */

#ifndef FLIPPY_CUSTOM_CONCEPTS_HPP
#define FLIPPY_CUSTOM_CONCEPTS_HPP
#include <concepts>
#include <iterator>
#include <numbers>
/**
 * @file
 * @brief This file contains the concepts that are costomly defined for the flippy class templates.
 */
namespace fp {
/**
 * @defgroup myConcepts Special concepts
 * @brief definitions of concepts that restrict template parameters in  flippy's classes.
 * @note For more information on what `c++` concepts are in general see
 * [this](https://en.cppreference.com/w/cpp/language/constraints) cpp reference page.
 *@{
 */

/**
 * @brief Here we implement the concepts of a floating point number.
 *
 * This concept is equivalent to [std::floating_point](https://en.cppreference.com/w/cpp/concepts/floating_point).
 * However, not all `c++20` compilers implement these features yet.
 * In particular some apple clang versions that are technically `c++20` capable.
 * @tparam T This concept requires the type T to be a floating point number.
 */
template <class T>
concept floating_point_number = std::is_floating_point_v<T>;

/**
 * @brief Here we implement the concepts of a positive integer number that is used throughout the code for indexing.
 *
 * @tparam T This concept requires the type T to be unsigned and an integral type.
 */
template <class T>
concept IndexingNumber = std::is_unsigned_v<T> && std::is_integral_v<T>;
/**@}*/

using Index = std::size_t;
using Real  = double;

namespace implementation {
template <typename C>
concept Beginable = requires(C c) { std::begin(c); };
template <typename C>
concept Endable = requires(C c) { std::end(c); };
template <typename C>
concept NeqableBeginAndEnd = requires(C c) {
    { std::begin(c) != std::end(c) } -> std::same_as<bool>;
};
template <typename C>
concept BeginIncrementable = requires(C c) { std::begin(c)++; };
template <typename C>
concept BeginDerefable = requires(C c) { *std::begin(c); };
template <typename C>
concept BeginDerefToVoid = requires(C c) {
    { *std::begin(c) } -> std::same_as<void>;
};
template <typename C>
concept BeginAndEndCopyConstructibleAndDestructible = requires(C c) {
    requires std::destructible<decltype(std::begin(c))> && std::destructible<decltype(std::end(c))> &&
                 std::copy_constructible<decltype(std::begin(c))> && std::copy_constructible<decltype(std::end(c))>;
};

} // namespace implementation
template <typename C>
concept Container =
    implementation::Beginable<C> && implementation::Endable<C> && implementation::NeqableBeginAndEnd<C> &&
    implementation::BeginIncrementable<C> && implementation::BeginDerefable<C> &&
    !implementation::BeginDerefToVoid<C> && implementation::BeginAndEndCopyConstructibleAndDestructible<C>;

[[maybe_unused]] static constexpr Real operator""_r(long double value) { return static_cast<Real>(value); }
[[maybe_unused]] static constexpr Real operator""_r(unsigned long long value) { return static_cast<Real>(value); }

static constexpr auto PI = static_cast<Real>(std::numbers::pi);
} // namespace fp

#endif // FLIPPY_CUSTOM_CONCEPTS_HPP
