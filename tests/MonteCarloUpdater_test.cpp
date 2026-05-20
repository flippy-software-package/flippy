#include "MonteCarloUpdater.hpp"
#include "flippy.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

namespace fp::test {

namespace {

struct FakeEnergyParameters {
    Real p;
};

auto surface_energy(const fp::Node &,
                    const fp::Triangulation &,
                    const FakeEnergyParameters &,
                    const std::vector<fp::Index> &) {
    return 0_r;
};

} // namespace

TEST_CASE("Test Monte Carlo Updater") {
    Index nIter = 4;

    auto                 trg = Triangulation::make_spherical_triangulation(nIter, 1., 0.);
    FakeEnergyParameters prms{.p = 1};

    std::mt19937 rng(0);

    SECTION("Monte Carlo Updater is initialized correctly") {

        auto l_min = 2_r;
        auto l_max = 2_r * l_min;
        auto mc_updater =
            MonteCarloUpdater<FakeEnergyParameters, std::mt19937>(trg, prms, surface_energy, rng, l_min, l_max);
        CHECK(mc_updater.kBT() == 1);
        CHECK(mc_updater.move_attempt_count() == 0);
        CHECK(mc_updater.bond_length_move_rejection_count() == 0);
        CHECK(mc_updater.move_back_count() == 0);
        CHECK(mc_updater.flip_attempt_count() == 0);
        CHECK(mc_updater.bond_length_flip_rejection_count() == 0);
        CHECK(mc_updater.flip_back_count() == 0);
    }
}

} // namespace fp::test
