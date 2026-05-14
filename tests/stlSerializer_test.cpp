#include "flippy.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

const std::string ASSETS_PATH{TEST_ASSET_PATH};

TEST_CASE("testing stlSerializer") {

    std::filesystem::path stl_solids{ASSETS_PATH};

    SECTION("basic serialization, test 1") {
        std::filesystem::path stl_solid_path{stl_solids.append("icosasphere.stl")};

        auto trg = fp::Triangulation::experimental_load_sphere_from_stl(stl_solid_path, 1.);

        CHECK(trg.size() == 42);
    }

    SECTION("basic serialization, test 2") {
        std::filesystem::path stl_solid_path{stl_solids.append("icosasphere_1082.stl")};

        auto trg = fp::Triangulation::experimental_load_sphere_from_stl(stl_solid_path, 1.);

        CHECK(trg.size() == 1082);
        double V  = trg.global_geometry().volume;
        double A  = trg.global_geometry().area;
        double Rv = std::pow(3. * V / (4. * M_PI), 1. / 3.);
        double Ra = std::sqrt(A / (4. * M_PI));
        CHECK(Rv == Catch::Approx(Ra).epsilon(0.01));
    }
}

