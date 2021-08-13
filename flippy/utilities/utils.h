#ifndef FLIPPY_UTILS_H
#define FLIPPY_UTILS_H
#include <iostream>
#include <chrono>
#include <ctime>
#include <fstream>

#if LOGGING_ON
#define LOGN(x) std::cout<< #x <<": "; print(x)
#else
#define LOGN(x)
#endif

#if LOGGING_ON
#define LOG(...) print(__VA_ARGS__)
#else
#define LOG(...);
#endif

#if PRINTING_ON
#define PRINT(...) print(__VA_ARGS__)
#else
#define PRINT(...)
#endif

namespace fp {

using Json = nlohmann::json;

static inline Json json_read(std::string file_name);
static inline void json_dump(std::string const& file_name, Json data);

static void print();

template<template<typename, typename> typename C, typename T>
static void print(const C<T, std::allocator<T>>& output);

template<typename Tfirst, typename... Trest>
static void print(const Tfirst& first, const Trest& ... rest);


// ------------------------ IMPLEMENTATIONS ---------------------------- //

static inline void json_dump(std::string const& file_name, Json data)
{
    std::ofstream o(file_name + ".json");
    o << data.dump();
    o.close();
}
Json inline json_read(std::string file_name)
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

[[maybe_unused]] void print() { std::cout << '\n'; }

template<template<typename, typename> typename C, typename T>
[[maybe_unused]] void print(const C<T, std::allocator<T>>& output)
{
    std::cout << "{";
    auto last = output.size() - 1;
    for (unsigned long i = 0; i<last; ++i) {
        std::cout << output[i] << ", ";
    }
    std::cout << output[last] << "}" << std::endl;
}

template<typename Tfirst, typename... Trest>
[[maybe_unused]] void print(const Tfirst& first, const Trest& ... rest)
{
    std::cout << first << ' ';
    print(rest...);
}

//template<typename T>
//[[maybe_unused]] inline bool is_member(std::vector<T> const& v, T const& el){
//  /***
//   * if the function returns true than el is contained in v (at least once).
//   */
//  return (std::find(v.begin(),v.end(), el) != v.end());
//}


class Timer
{
/**
 * this class keeps time form its creation to destrtuction. I.e. it can
 * time the duration of a scope.
 *
 * */
private:
    using Clock =
    std::conditional_t<std::chrono::high_resolution_clock::is_steady,
                       std::chrono::high_resolution_clock,
                       std::chrono::steady_clock>;
    const Clock::time_point Start = Clock::now();
    Clock::time_point Now;
public:
    std::string unit;
    std::string unit_fine;
    std::string scope_name = "GLOBAL (guess)";
    [[maybe_unused]] Timer() = default;
    [[maybe_unused]] explicit Timer(std::string const& _scope_name)
            :scope_name(_scope_name) { };
    void stop()
    {
        Now = Clock::now();
        auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>
                (Now - Start).count();
        auto diff_fine = diff;

        if (diff/(long double)1e3<1ll) { unit = " ns"; }
        else if (diff/((long double)1e6)<1) {
            unit = " µs";
            diff /= 1e3;
        }
        else if (diff/((long double)1e9)<1) {
            unit = " ms";
            diff /= 1e6;
            unit_fine = " µs";
            diff_fine /= 1e3;
        }
        else if (diff/((long double)1e9*60)<1) {
            unit = " s";
            diff /= 1e9;
            unit_fine = " ms";
            diff_fine /= 1e6;
        }
        else if (diff/((long double)1e9*3600)<1) {
            unit = " m";
            diff /= 6e10;
            unit_fine = " s";
            diff_fine /= 1e9;
        }
        else {
            unit = " h";
            diff /= 36e11;
            unit_fine = " m";
            diff_fine /= 6e10;
        }

        auto end = std::chrono::system_clock::now();
        auto end_time = std::chrono::system_clock::to_time_t(end);
        std::cout << "this timer ran in the following scope: " << scope_name << '\n';
        std::cout << "finished computation at " << std::ctime(&end_time)
                  << "elapsed time: " << diff << unit << " (" << diff_fine << unit_fine << ")" << '\n';
    }

    ~Timer() { stop(); }
};
}

#endif