#pragma once
#include <iostream>
#include <chrono>
#include <ctime>
#include <type_traits>
#include <concepts>

namespace cutils {

#if LOGGING_ON
#define LOGN(x) std::cout<< #x <<": "; cutils::print(x)
#else
#define LOGN(x)
#endif

#if LOGGING_ON
#define LOG(...) cutils::print(__VA_ARGS__)
#else
#define LOG(...);
#endif

#if PRINTING_ON
#define PRINT(...) cutils::print(__VA_ARGS__)
#else
#define PRINT(...)
#endif


template<typename C> concept Beginable= requires(C c) {std::begin(c);};
template<typename C> concept Endable= requires(C c) {std::end(c);};
template<typename C> concept NeqableBeginAndEnd = requires(C c) {{std::begin(c) != std::end(c)} -> std::same_as<bool>;};
template<typename C> concept BeginIncrementable = requires(C c) {std::begin(c)++;};
template<typename C> concept BeginDerefable = requires(C c) {*std::begin(c);};
template<typename C> concept BeginDerefToVoid = requires(C c) {{*std::begin(c)} -> std::same_as<void>;};
template<typename C> concept BeginAndEndCopyConstructibleAndDestructible = requires(C c) {
  std::destructible<decltype(std::begin(c))> &&
      std::destructible<decltype(std::end(c))> &&
      std::copy_constructible<decltype(std::begin(c))> &&
      std::copy_constructible<decltype(std::end(c))>;
};


template<typename C> concept Container =
Beginable<C> &&
    Endable<C> &&
    NeqableBeginAndEnd<C> &&
    BeginIncrementable<C> &&
    BeginDerefable<C> &&
    !BeginDerefToVoid<C> &&
    BeginAndEndCopyConstructibleAndDestructible<C>;

template<typename T> concept StdOutStreamable = requires(T t) {std::cout << t;};
template<typename T> concept OutStreamable = requires(std::ostream os, T t) {os << t;};
template<typename T> concept InStreamable = requires(std::istream os, T t) {os >> t;};

template<typename T> concept Printable = StdOutStreamable<T> || Container<T>;


template<char sep, char end> static void print();
template<char sep, char end> static void print(Container auto const & output);
template<char sep, char end> static void print(Printable auto const & first, Printable auto const & ... rest);


struct HumanReadableTime{
  std::string unit;
  std::string unit_fine;
  long long diff;
  long long diff_fine;
  long long diff_ns;
};

// ------------------------ IMPLEMENTATIONS ---------------------------- //



template<char sep=' ', char end='\n'>
[[maybe_unused]] void print() { std::cout << end; }

template<char sep=' ', char end='\n'>
[[maybe_unused]] void print(Container auto const & output)
{
  print<sep, ' '>('{');
  auto const & last_elem = *(output.end()-1);
  for (const auto& elem : output) {
    if (elem != last_elem){
      print<sep, ' '>(elem, ',');
    }
    else{
      print<sep, end>(elem, '}');
    }
  }
}


/**
 * Print function that can take an arbitrary number of arguments. It was implemented to work similarly to Pythons print function.
 * @tparam Tfirst Type of the first argument (automatically deduced)
 * @tparam Trest Type of the remaining part of the variadic number of arguments (automatically deduced)
 * @param first First argument
 * @param rest rest of the variadic arguments
 * @returns nothing.
 * @see print(const C<T, std::allocator<T>>& output)
 * @attention can not handle containers of containers.
 * @warning if a container of zero length is provided the code will crush.
 */
//template<typename Tfirst, typename... Trest>
//[[maybe_unused]] void print(const Tfirst& first, const Trest& ... rest)
template<char sep=' ', char end='\n'>
[[maybe_unused]] void print(Printable auto const & first, Printable auto const & ... rest)
{
  if constexpr (StdOutStreamable<decltype(first)>){
      std::cout<<first<<sep;
  }
  else{
      print<sep, ' '>(first);
  }
  print<sep, end>(rest...);
}

/**
 * given a time difference in the highest clock resolution this struct constructs a human readable string of elapsed time.
 */
static HumanReadableTime human_readable_time(long long diff){
    std::string unit;
    std::string unit_fine;
    long long diff_ns = diff;
    long long diff_fine = diff;
    double diff_d = static_cast<double>(diff);
    if (diff_d/(1.e3l)<1.l) { unit = " ns"; }
    else if (diff_d/1.e6<1.) {
        unit = " µs";
        diff /= 1000;
    }
    else if (diff_d/1.e9<1.) {
        unit = " ms";
        diff /= 1000'000ll;
        unit_fine = " µs";
        diff_fine /= 1000ll;
    }
    else if (diff_d/(1.e9l*60.l)<1.l) {
        unit = " s";
        diff /= 1000'000'000ll;
        unit_fine = " ms";
        diff_fine /= 1000'000ll;
    }
    else if (diff_d/(1.e9*3600.)<1.) {
        unit = " m";
        diff /= 60'000'000'000ll;
        unit_fine = " s";
        diff_fine /= 1000'000'000ll;
    }
    else {
        unit = " h";
        diff /= 3600'000'000'000ll;
        unit_fine = " m";
        diff_fine /= 60'000'000ll;
    }
    return {.unit=unit, .unit_fine=unit_fine, .diff=diff,
            .diff_fine=diff_fine, .diff_ns=diff_ns};
}

class Timer
{
/**
 * this class keeps time form its creation to destruction. I.e. it can
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
        long long diff_ = std::chrono::duration_cast<std::chrono::nanoseconds>(Now - Start).count();

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


template<typename G>
concept UniformRandomBitGenerator = requires(G g) {
  typename G::result_type; // Must have result_type
  { G::min() } -> std::same_as<typename G::result_type>; // min function
  { G::max() } -> std::same_as<typename G::result_type>; // max function
  { g() } -> std::same_as<typename G::result_type>; // operator()
  requires std::is_unsigned_v<typename G::result_type>; // result_type is unsigned integer type
  requires G::min() < G::max(); // min is strictly less than max
};


template <typename E>
concept RandomNumberEngine = requires(E e, const E x, const E y, typename E::result_type s, unsigned long long z, std::ostream& os, std::istream& is) {
  typename E::result_type;
  { E() } -> std::same_as<E>;
  { E(x) } -> std::same_as<E>;
  { E(s) } -> std::same_as<E>;
  { e.seed() } -> std::same_as<void>;
  { e.seed(s) } -> std::same_as<void>;
  { e() } -> std::same_as<typename E::result_type>;
  { e.discard(z) } -> std::same_as<void>;
  { x == y } -> std::same_as<bool>;
  { x != y } -> std::same_as<bool>;
  { os << x } -> std::same_as<std::ostream&>;
  { is >> e } -> std::same_as<std::istream&>;
  { E::min() } -> std::same_as<typename E::result_type>;
  { E::max() } -> std::same_as<typename E::result_type>;
};


// source https://en.wikipedia.org/wiki/Xorshift#xoshiro_and_xoroshiro
//template<typename ResultType=uint32_t>
//requires std::is_integral<ResultType>::value
//requires std::is_same_v<ResultType, uint32_t>
class xorshift32
{
  using ResultType = uint32_t;
  ResultType seed_;
public:
  typedef ResultType result_type;

  xorshift32(): seed_(12) { }

  explicit xorshift32(ResultType seed_): seed_(seed_) { }

  /* The state word must be initialized to non-zero */
  ResultType operator()()
  {
    /* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */
    ResultType x = seed_;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return seed_ = x;
  }
  void discard(long long z)
  {
    for (long long i = 0; i<z; ++i) {
      (*this)();
    }
  }
  constexpr static ResultType min() { return std::numeric_limits<ResultType>::min(); }
  constexpr static ResultType max() { return std::numeric_limits<ResultType>::max(); }

  void seed(){*this = xorshift32();}
  void seed(ResultType seed_inp){*this = xorshift32(seed_inp);}

  bool operator==(xorshift32 const& rhs) const{
    return seed_==rhs.seed_;
  }
  bool operator!=(xorshift32 const& rhs) const{
    return seed_!=rhs.seed_;
  }


  friend auto operator<<( std::ostream& os, xorshift32 const& rng) -> std::ostream& {
    os.flags(std::ostream::dec | std::ostream::skipws);
    os << rng.seed_;
    return os;
  }

  friend auto operator>>(std::istream& is, xorshift32 & rng) -> std::istream & {
    is.flags(std::istream::dec | std::istream::skipws);
    is >> rng.seed_;
    return is;
  }

};


}
