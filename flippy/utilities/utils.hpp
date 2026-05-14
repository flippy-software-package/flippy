#ifndef FLIPPY_UTILITIES_UTILS_HPP
#define FLIPPY_UTILITIES_UTILS_HPP
/** @file
 *  @brief This file contains helper functions that are used throughout flippy, but are not specific to any given class.
 */

#include <external/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <type_traits>
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
