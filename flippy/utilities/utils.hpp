#ifndef FLIPPY_UTILS_H
#define FLIPPY_UTILS_H
/*
 *  helper functions for flippy
 */

#include <iostream>
#include <fstream>
#include <utility>
#include <filesystem>
#include <type_traits>

namespace fp {

using Json = nlohmann::json;

[[maybe_unused]] static inline void json_dump(std::string const& file_name, const Json& data)
{
    std::ofstream o(file_name + ".json");
    o << data.dump();
    o.close();
}
[[maybe_unused]] static Json inline json_read(std::string file_name)
{
    auto pos_json = file_name.find_last_of(".json");
    auto not_json = (file_name.size() - 1!=pos_json);
    if (not_json) { file_name = file_name + ".json"; }
    std::ifstream o(file_name);
    Json data;
    o >> data;
    o.close();
    return data;
}

template<typename T>
[[maybe_unused]] static bool is_member(std::vector<T> const& v, T const& el){
    /***
     * if the function returns true, then el is contained in v (at least once).
     */
    return (std::find(v.begin(),v.end(), el) != v.end());
}

[[maybe_unused]] static bool is_dir(std::string const& path){
    std::filesystem::path dir(path);
    return std::filesystem::is_directory(dir);
}

[[maybe_unused]] static void make_dir(std::string const& path){
    std::filesystem::path dir(path);
    std::filesystem::create_directories(dir);
}

}
#endif