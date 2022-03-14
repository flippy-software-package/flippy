#ifndef FLIPPY_DEBUG_UTILS_H
#define FLIPPY_DEBUG_UTILS_H
#include <iostream>
#include <chrono>
#include <ctime>
#include <utility>

namespace cutils {

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

static void print();

template<template<typename, typename> typename C, typename T>
static void print(const C<T, std::allocator<T>>& output);

template<typename Tfirst, typename... Trest>
static void print(const Tfirst& first, const Trest& ... rest);

struct HumanReadableTime{
  std::string unit;
  std::string unit_fine;
  unsigned long long diff;
  unsigned long long diff_fine;
  unsigned long long diff_ns;
};

// ------------------------ IMPLEMENTATIONS ---------------------------- //

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

static HumanReadableTime human_readable_time(unsigned long long diff){
    std::string unit;
    std::string unit_fine;
    unsigned long long diff_ns = diff;
    unsigned long long diff_fine = diff;
    if (diff/(1.e3l)<1.l) { unit = " ns"; }
    else if (diff/1.e6l<1.l) {
        unit = " µs";
        diff /= 1.e3l;
    }
    else if (diff/1.e9l<1.l) {
        unit = " ms";
        diff /= 1.e6l;
        unit_fine = " µs";
        diff_fine /= 1.e3l;
    }
    else if (diff/(1.e9l*60.l)<1.l) {
        unit = " s";
        diff /= 1e9l;
        unit_fine = " ms";
        diff_fine /= 1e6l;
    }
    else if (diff/(1.e9l*3600.l)<1.l) {
        unit = " m";
        diff /= 6.e10l;
        unit_fine = " s";
        diff_fine /= 1.e9l;
    }
    else {
        unit = " h";
        diff /= 36.e11l;
        unit_fine = " m";
        diff_fine /= 6.e10l;
    }
    return {.unit=unit, .unit_fine=unit_fine, .diff=diff,
            .diff_fine=diff_fine, .diff_ns=diff_ns};
}

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
    Clock::time_point Start = Clock::now();
    Clock::time_point Now;
    bool stopped=false;
public:
    [[maybe_unused]] Timer() = default;
    void restart(){Start = Clock::now();}
    HumanReadableTime stop()
    {
        Now = Clock::now();
        unsigned long long diff_ = std::chrono::duration_cast<std::chrono::nanoseconds>
                (Now - Start).count();

        HumanReadableTime hrt = human_readable_time(diff_);

        auto end = std::chrono::system_clock::now();
        auto end_time = std::chrono::system_clock::to_time_t(end);
        std::cout << "finished computation at " << std::ctime(&end_time)
        << "elapsed time: " << hrt.diff << hrt.unit << " (" << hrt.diff_fine << hrt.unit_fine << ")" << '\n';
        stopped = true;
        return hrt;
    }

    ~Timer() { if(not stopped){stop();} }
};

}
#endif
