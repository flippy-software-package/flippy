#ifndef FLIPPY_CUSTOM_CONCEPTS_HPP
#define FLIPPY_CUSTOM_CONCEPTS_HPP
#include <concepts>

namespace fp{
/**
 * Here we implement the concepts of a floating point number and integral number.
 * Since apple clang does not implement concepts fully we need to do this here.
 */
template<class T> concept floating_point_number = std::is_floating_point_v<T>;
template<class T> concept integer_number = std::is_integral_v<T>;
}


#endif //FLIPPY_CUSTOM_CONCEPTS_HPP
