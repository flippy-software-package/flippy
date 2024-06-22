#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "test_helper.h"
#include "utilities/DataIO.hpp"

using REAL = float;
using INDEX = unsigned int;
fp::DataIOInterface<REAL, INDEX> auto reset_memory_and_give_back_copy(fp::DataIOInterface<REAL, INDEX> auto dio){
    dio.reset_memory();
    return dio;
}
TEST_CASE("construct xyzDataIO"){
    fp::Triangulation<REAL, INDEX> trg;

    auto dio = fp::xyzDataIO("./");
    SECTION("xyzDataIO matches the interface"){
        dio.reset_memory();
        auto dio_copy = reset_memory_and_give_back_copy(dio);
        STATIC(dio_copy==dio);
    }
}