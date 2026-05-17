#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "Triangulation.hpp"
#include "test_helper.h"
#include "utilities/DataIO.hpp"

using REAL  = float;
using INDEX = unsigned int;
fp::experimental::DataSaver reset_memory_and_give_back_copy(fp::DataIOInterface auto dio) {
    dio.reset_memory();
    return dio;
}
TEST_CASE("construct xyzDataIO") {

    fp::Triangulation::experimental_load_sphere_from_stl(const std::filesystem::path &stl_file_path,
                                                         Real verlet_radius_inp) fp::Triangulation trg(1);

    auto dio = fp::experimental::xyzDataSaver{};
    SECTION("xyzDataIO matches the interface") {
        dio.reset_memory();
        auto dio_copy = reset_memory_and_give_back_copy(dio);
        STATIC(dio_copy == dio);
    }
}
