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

#ifndef FLIPPY_UTILITIES_UTILS_HPP
#define FLIPPY_UTILITIES_UTILS_HPP
/** @file
 *  @brief This file contains helper functions that are used throughout flippy, but are not specific to any given class.
 */

#include <external/json.hpp>
#include <fstream>
#include <utility>

namespace fp {
/**
 * @GlobalsStub
 * @{
 */
//! shortening of the nlohmann::json namespace, which is an [external open source
//! library](https://github.com/nlohmann/json) bundled by flippy.
using Json = nlohmann::json;

/**
 * @brief Simple wrapper function around Json objects built in dump() method.
 * @param file_name @FileNameOrPathFileNameStub
 * @param data json data object that is supposed to be stored.
 */
static inline void json_dump(const std::string &file_name, const Json &data) {
    std::ofstream o(file_name + ".json");
    o << data.dump();
    o.close();
}

/**
 * @brief Simple wrapper function  that reads the content of a text file into a json object.
 *
 * The file name onb the disk needs to end in '.json' for this function to work.
 * @param file_name @FileNameOrPathFileNameStub
 * @return Json object that was parsed from the provided file.
 *
 * @warning This function will stream any file into the json object.
 * If the provided file is not a valid json file this will cause runtime errors.
 */
static Json inline json_read(std::string file_name) {
    auto pos_json = file_name.find_last_of(".json");
    auto not_json = (file_name.size() - 1 != pos_json);
    if (not_json) { file_name = file_name + ".json"; }
    std::ifstream o(file_name);
    Json          data;
    o >> data;
    o.close();
    return data;
}

/**
 * @brief Convenient wrapper around std::find, which only works for std::vectors.
 *
 * @tparam T type of the vector elements.
 * @param v std::vector in which we want to search.
 * @param el the value of the element that we want to check for.
 * @return The function returns `true` if `el` is contained in vector `v` (at least once), otherwise it returns `false`.
 */
template <typename T> [[maybe_unused]] static bool is_member(const std::vector<T> &v, const T &el) {
    return (std::find(v.begin(), v.end(), el) != v.end());
}
/**@}*/
} // namespace fp
#endif // FLIPPY_UTILITIES_UTILS_HPP
